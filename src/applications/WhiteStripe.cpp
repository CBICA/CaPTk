#include "WhiteStripe.h"
#include "ws_utils.h"
WhiteStripe::WhiteStripe()
{
  m_wsWidth = 0.05;
  m_sliceStartZ = 80;
  m_sliceStopZ = 120;
  m_tissuesMax = 5;
  m_smoothMax = 10.0;
  m_smoothDelta = 0.5;
  m_histSize = 2000;
  m_bT1 = true;
}

void WhiteStripe::getHisInfo(std::vector<float>& mids, std::vector<float>& origHist, std::vector<float>& smoothHist, std::vector<int>& peakIds, int& modeId)
{
  m_hstObj.getHisInfo(mids, origHist, smoothHist, peakIds);
  if (mids.empty()) return;
  float mode;
  if (m_bT1)
  {
    mode = m_hstObj.getLastMode(false);
  }
  else
  {
    mode = m_hstObj.getLargestMode();
  }
  modeId = std::find(mids.begin(), mids.end(), mode) - mids.begin();
  return;
}
void WhiteStripe::setParams(float wsWidth, int sliceStartZ, int sliceStopZ, int tissuesMax, float smoothMax, float smoothDelta, int histSize, bool bT1)
{
  m_wsWidth = wsWidth;
  m_sliceStartZ = sliceStartZ;
  m_sliceStopZ = sliceStopZ;
  m_tissuesMax = tissuesMax;
  m_smoothMax = smoothMax;
  m_smoothDelta = smoothDelta;
  m_histSize = histSize;
  m_bT1 = bT1;
}
ImageType::Pointer WhiteStripe::process(ImageType::Pointer img, ImageType::Pointer& maskOut)
{
  ImageType::Pointer normalizedImg;
  try
  {

    vector<float> voi = makeImageVoi(img);
    if (voi.empty())
    {
      return normalizedImg;
    }
    m_hstObj.setParams(m_tissuesMax, m_smoothMax, m_smoothDelta, m_histSize);
    m_hstObj.compute(voi);
    float mode = 0.0;
    if (m_bT1)
    {
      mode = m_hstObj.getLastMode();
    }
    else
    {
      mode = m_hstObj.getLargestMode();
    }

    //mode = 84.22038;
    float modeQ = threshMean(voi, mode);//Find Fraction of pixels in voi less than threshold( mode) 
    //cout << "mode =" << mode << endl;
    //cout << "modeQ =" << modeQ << endl;
    std::vector<float>  whiteStripe, probs;
    float probsMin = max(modeQ - m_wsWidth, 0.0f);
    float probsMAx = min(modeQ + m_wsWidth, 1.0f);
    probs.push_back(probsMin);
    probs.push_back(probsMAx);
    whiteStripe = Utils::quantile(voi, probs);
    std::vector<ImageType::IndexType> whiteStripeInd = Utils::getPixelIndicesThresholded(img, whiteStripe[0], whiteStripe[1]);
    if (whiteStripeInd.empty())
    {
      //using whole brain normalization
      float mean = Utils::getMeanValue(img);
      whiteStripeInd = Utils::getPixelIndicesThresholded(img, mean);
    }
    //float whiteStripeSigma =calcStdDev(img, whiteStripeInd);//TBD used nowhere
    normalizedImg = whitestripeNorm(img, whiteStripeInd);
    maskOut = whitestripeIndToMask(img, whiteStripeInd);
  }
  catch (...)
  {
    cout << "Exception in process" << endl;
    return NULL;
  }
  return normalizedImg;
}

float WhiteStripe::threshMean(const std::vector<float>& voi, float thresh)
{
  double sum = 0;
  for (size_t i = 0; i < voi.size(); i++)
  {
    if (voi[i] < thresh)
    {
      sum = sum + 1;
    }
  }
  return float(sum / voi.size());
}
std::vector<float> WhiteStripe::makeImageVoi(const ImageType::Pointer img)
{

  double mean = Utils::getMeanValue(img);
  if (std::isnan(mean))
  {
    cout << "Error: Invalid image with mean  NAN \n ";
    return vector<float>();
  }
  // TBD something wrong some how img now giving different mean when used without duplicate
  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(img);
  duplicator->Update();
  ImageType::Pointer imgDup = duplicator->GetOutput();

  typedef itk::RegionOfInterestImageFilter< ImageType, ImageType> ExtractFilterType;
  ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
  ImageType::SizeType size3D = img->GetLargestPossibleRegion().GetSize();
  if (m_sliceStartZ < 0)// if unspecified, use Full ROI
  {
    m_sliceStartZ = 1;//Only for R compatability 
  }
  if (m_sliceStopZ < 0)// if unspecified, use Full ROI
  {
    m_sliceStopZ = size3D[2] - 1;
  }

  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = m_sliceStartZ - 1;// -1 as the index In R

  ImageType::SizeType size;
  size[0] = size3D[0];
  size[1] = size3D[1];
  size[2] = (m_sliceStopZ - start[2]);
  if (start[2]<0 || (start[2] + size[2])>size3D[2])
  {
    cout << "Error: Z slice out of range. sliceStartZ: " << m_sliceStartZ << ", sliceStopZ: " << m_sliceStopZ;
    cout << ", Image Z size: " << size3D[2] << endl;
    return vector<float>();
  }

  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize(size);
  desiredRegion.SetIndex(start);
  extractFilter->SetRegionOfInterest(desiredRegion);
  extractFilter->SetInput(imgDup);
  extractFilter->Update();
  ImageType::Pointer imgVoi = ImageType::New();
  imgVoi->Graft(extractFilter->GetOutput());
  //cout << imgVoi->GetLargestPossibleRegion().GetSize() << ", mean=" << getMeanValue(img) << ", MeanVoi=" << getMeanValue(imgVoi) << endl;

  std::vector<float> voi;
  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  ConstIteratorType in(imgVoi, imgVoi->GetRequestedRegion());
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    float val = in.Get();
    if (val > mean)
    {
      voi.push_back(val);
    }
  }
  return voi;
}




double WhiteStripe::calcStdDev(const ImageType::Pointer& img, const std::vector<ImageType::IndexType>& indices)
{
  std::vector<float> pixelVals;
  for (size_t i = 0; i < indices.size(); i++)
  {
    pixelVals.push_back(img->GetPixel(indices[i]));
  }
  double sum = std::accumulate(pixelVals.begin(), pixelVals.end(), 0.0);
  double mean = sum / pixelVals.size();
  double sqSum = std::inner_product(pixelVals.begin(), pixelVals.end(), pixelVals.begin(), 0.0);
  double stdev = std::sqrt(sqSum / pixelVals.size() - mean * mean);
  return stdev;
}
void WhiteStripe::zeroOutImg(ImageType::Pointer inputImage)
{
  typedef itk::ImageRegionIterator< ImageType > IteratorType;
  IteratorType itr(inputImage, inputImage->GetRequestedRegion());
  for (itr.GoToBegin(); !itr.IsAtEnd(); ++itr)
  {
    itr.Set(0);
  }
}
ImageType::Pointer WhiteStripe::whitestripeIndToMask(ImageType::Pointer img, std::vector<ImageType::IndexType> indices)
{
  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(img);
  duplicator->Update();
  ImageType::Pointer mask = duplicator->GetOutput();
  zeroOutImg(mask);
  for (size_t i = 0; i < indices.size(); i++)
  {
    mask->SetPixel(indices[i], 1);
  }
  return mask;
}
std::vector<ImageType::IndexType> WhiteStripe::whitestripeMasktoInd(ImageType::Pointer mask)
{
  std::vector<ImageType::IndexType> indices;
  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  ConstIteratorType in(mask, mask->GetRequestedRegion());
  for (in.GoToBegin(); !in.IsAtEnd(); ++in)
  {
    if (in.Get() > 0)
    {
      indices.push_back(in.GetIndex());
    }
  }
  return indices;
}
ImageType::Pointer WhiteStripe::whitestripeNorm(ImageType::Pointer img, std::vector<ImageType::IndexType> indices) {
  std::vector<float> pixelVals;
  for (size_t i = 0; i < indices.size(); i++)
  {
    pixelVals.push_back(img->GetPixel(indices[i]));
  }
  double sum = std::accumulate(pixelVals.begin(), pixelVals.end(), 0.0);
  double mean = sum / pixelVals.size();
  double sq_sum = std::inner_product(pixelVals.begin(), pixelVals.end(), pixelVals.begin(), 0.0);
  double stdev = std::sqrt(sq_sum / pixelVals.size() - mean * mean);

  //img = (img - mu) / sig  :that is:   img = (img - mean) / stdev
  typedef itk::ImageRegionIterator< ImageType > IteratorType;
  IteratorType out(img, img->GetRequestedRegion());
  for (out.GoToBegin(); !out.IsAtEnd(); ++out)
  {
    float val = (out.Get() - mean) / stdev;
    out.Set(val);
  }
  //TBD img = cal_img(img); img = zero_trans(img);
  return img;
}