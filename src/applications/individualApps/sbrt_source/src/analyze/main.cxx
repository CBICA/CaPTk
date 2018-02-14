#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <exception>
#include <typeinfo>
#include <stdexcept>
#include <stdio.h>
#include <numeric>
#include <functional>
#include <vector>
#include "generateTextureFeatures.h"
#include "LBPFeatures.h"
#include "cbicaCmdParser.h"
#include "cbicaLogging.h"
#include "FeatureScalingClass.h"

#include "cbicaITKSafeImageIO.h"
#include "cbicaUtilities.h"
//#include "cbicaPreProcessImage.h"
#include "itkImageFileReader.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include <itkVectorImage.h>
#include <itkDenseFrequencyContainer2.h>
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorContainer.h"
#include "itkStatisticsImageFilter.h"
//#include "opencv/cv.h"
#include "opencv2/core/core.hpp"
#include "itkShiftScaleImageFilter.hxx"
//#include "itkOpenCVImageBridge.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryMorphologicalOpeningImageFilter.h"
#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkOtsuMultipleThresholdsImageFilter.h"
#include "itkPasteImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkBinaryContourImageFilter.h"
#include "itkNormalizeImageFilter.h"

#include "opencv2/core/core.hpp"
#include "opencv2/ml.hpp"
#include <math.h>
#include <cmath>


typedef  float           PixelType;
typedef itk::Image< PixelType, 3 > ImageType;
typedef ImageType::Pointer InputImagePointerType;
typedef itk::ImageFileReader< ImageType> ImageReaderType;
typedef ImageReaderType::Pointer InputImageFileReaderPointerType;

typedef itk::ShiftScaleImageFilter<
  ImageType, ImageType >  ShiftScaleFilterType;
typedef itk::Neighborhood<float, 3> NeighborhoodType;
typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<ImageType>Image2CoOccuranceType;
typedef Image2CoOccuranceType::HistogramType HistogramType;
typedef itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType> Hist2FeaturesType1;
typedef ImageType::OffsetType OffsetType;
typedef itk::VectorContainer< unsigned char, OffsetType > OffsetVector;
typedef OffsetVector::Pointer  OffsetVectorPointer;
typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> roiType;
typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;
typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;
typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter
<ImageType, HistogramFrequencyContainerType> RunLengthFilterType;

typedef short                                    RunLengthFeatureName;
typedef itk::VectorContainer<unsigned char, RunLengthFeatureName>  FeatureNameVector;
typedef FeatureNameVector::Pointer      FeatureNameVectorPointer;
typedef FeatureNameVector::ConstPointer FeatureNameVectorConstPointer;
typedef itk::VectorContainer< unsigned char, double > FeatureValueVector;
typedef FeatureValueVector::Pointer     FeatureValueVectorPointer;
std::fstream FileToWorkWith;
typedef itk::Image< itk::Vector<PixelType, 2>, Dimension > MaskImageType;
typedef itk::ImageRegionIteratorWithIndex< MaskImageType > IteratorType_withmask;

typedef itk::VariableSizeMatrix< double > MatrixType;
MatrixType Featurematrix;


VariableSizeMatrixType ScaleGivenTrainingFeatures(VariableSizeMatrixType &inputdata)
{
  unsigned int NumberOfSamples = inputdata.Rows();
  unsigned int NumberOfFeatures = inputdata.Cols() - 1;
  std::vector<float> feature;
  double minValue = inputdata[0][0], maxValue = inputdata[0][0];
  for (unsigned int y = 0; y < NumberOfFeatures; y++) // for each row
    for (unsigned int x = 0; x <NumberOfSamples; x++) { // for each element of row
      minValue = min(minValue, inputdata[x][y]);
      maxValue = max(maxValue, inputdata[x][y]);
    }
  //change the matrix values 
  for (unsigned int y = 0; y < NumberOfFeatures; y++)
    for (unsigned int x = 0; x <NumberOfSamples; x++)
      inputdata[x][y] = (inputdata[x][y] - minValue) / (maxValue - minValue);

  return inputdata;
}


int read_data_from_csv(string filename, cv::Mat data, int n_samples)
{

  std::string token;
  std::string line1;

  std::ifstream f(filename);
  if (!f)
  {
    printf("ERROR: cannot read file %s\n", filename.c_str());
    return 0; // all not OK
  }
  int attribute = 0;
  int line = 0;

  while (std::getline(f, line1)){
    std::stringstream convertor(line1);
    while (std::getline(convertor, token, ','))
    {
      //  std::cout << token;
      data.at<float>(line, attribute) = atof(token.c_str());
      attribute = attribute + 1;
    }
    line = line + 1;
    attribute = 0;
  }
  return 1; // all OK
}


template <class TImageType = ImageTypeFloat3D>
void Intensity_featurs(typename TImageType::Pointer image, typename TImageType::Pointer mask)
{

}


template<typename TPixel, unsigned int VDimension>
void getFeature(itk::Image<TPixel, VDimension> *imagein,
  itk::Image<TPixel, VDimension>* PET_imgin, itk::Image<TPixel, VDimension>* maskin, std::vector<std::string> *featurename,
  std::vector<float> *featurevec)
{
  try
  {
    if (maskin->GetLargestPossibleRegion() != PET_imgin->GetLargestPossibleRegion())
    {
      return;
    }

    typedef itk::NormalizeImageFilter < ImageType, ImageType >
      NormalizeFilterType;
    NormalizeFilterType::Pointer normalizeFilter = NormalizeFilterType::New();
    normalizeFilter->SetInput(PET_imgin);
    normalizeFilter->Update();
    ImageType::Pointer PETimg_windowed = normalizeFilter->GetOutput();

    typedef itk::MedianImageFilter<ImageType, ImageType > FilterType;
    FilterType::Pointer medianFilter = FilterType::New();
    typedef  itk::ImageFileWriter< ImageType  > WriterType;
    WriterType::Pointer writer = WriterType::New();
    //FilterType::InputSizeType radius1;
    typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
    MaskFilterType::Pointer maskFilter = MaskFilterType::New();


    /*Extract intnsity based features like mean , variance , skewness , kurtosis*/
    typedef MaskImageType::Pointer NewMask;
    NewMask newmask = MaskImageType::New();
    MaskImageType::RegionType region1;

    region1.SetSize(PETimg_windowed->GetLargestPossibleRegion().GetSize());
    newmask->SetRegions(region1);
    newmask->Allocate();
    MaskImageType::PixelType  initialValue;
    initialValue.Fill(0.0);

    newmask->FillBuffer(initialValue);

    IteratorType inputIt(maskin, maskin->GetLargestPossibleRegion());
    std::vector< ImageType::IndexType> index_vec;
    std::vector< ImageType::PixelType > nonzero_pix;
    for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    {

      ImageType::IndexType idx = inputIt.GetIndex();
      if (idx[0] << maskin->GetLargestPossibleRegion().GetSize(0) && idx[1] << maskin->GetLargestPossibleRegion().GetSize(1) &&
        idx[2] << maskin->GetLargestPossibleRegion().GetSize(2))
      {
        MaskImageType::PixelType   pixelValue;
        pixelValue[0] = PETimg_windowed->GetPixel(idx);   // x component
        pixelValue[1] = maskin->GetPixel(inputIt.GetIndex());
        newmask->SetPixel(inputIt.GetIndex(), pixelValue);
      }
      //  std::cout << newmask->GetPixel(inputIt.GetIndex());

      if (maskin->GetPixel(inputIt.GetIndex()) != 0)
      {
        index_vec.push_back(inputIt.GetIndex());
        nonzero_pix.push_back(PETimg_windowed->GetPixel(inputIt.GetIndex()));
      }
    }
    double min = *min_element(nonzero_pix.begin(), nonzero_pix.end());
    //std::cout << "Min value: " << min;
    double max = *max_element(nonzero_pix.begin(), nonzero_pix.end());
    //std::cout << "Max value: " << max;

    double sum = std::accumulate(nonzero_pix.begin(), nonzero_pix.end(), 0.0);
    double mean = sum / nonzero_pix.size();

    std::vector<double> diff(nonzero_pix.size());
    std::transform(nonzero_pix.begin(), nonzero_pix.end(), diff.begin(),
      std::bind2nd(std::minus<double>(), mean));
    //  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    //  double stdev = std::sqrt(sq_sum / nonzero_pix.size());
    //std::cout << "Mean value: " << mean;
    //std::cout << "std value: " << stdev;
    auto tempSpacing = PETimg_windowed->GetSpacing();
    double voxvol = tempSpacing[0] * tempSpacing[1] * tempSpacing[2];
    double volume = nonzero_pix.size() * voxvol;
    //std::cout << "volume: " << volume;
    const int numlevels = 32; std::vector<double> gray_lim = { min, max };
    float slope = (numlevels - 1) / (gray_lim[1] - gray_lim[0]);
    float intercept = 1 - (slope*gray_lim[0]);
    std::transform(nonzero_pix.begin(), nonzero_pix.end(), nonzero_pix.begin(),
      std::bind1st(std::multiplies<float>(), slope));
    transform(nonzero_pix.begin(), nonzero_pix.end(), nonzero_pix.begin(), bind2nd(std::plus<float>(), intercept));
    std::vector< float > hist(numlevels, 0);
    //std::cout << ceil(nonzero_pix.at(1));

    for (unsigned int j = 0; j < nonzero_pix.size(); j++)
    {
      int dummy = std::round(nonzero_pix[j]);
      hist[dummy - 1]++;
    }

    std::transform(hist.begin(), hist.end(), hist.begin(),
      std::bind2nd(std::divides<float>(), nonzero_pix.size()));
    std::vector< float > hist_grey(numlevels, 0);
    std::vector< float > grey(numlevels, 0);
    std::vector< float > hist_var(numlevels, 0);
    std::vector< float > hist_sk(numlevels, 0);
    std::vector< float > hist_kur(numlevels, 0);
    for (int j = 0; j < numlevels; j++)
    {
      hist_grey.at(j) = hist.at(j)*(j + 1);
      grey.at(j) = j + 1;
    }
    double histmean = std::accumulate(hist_grey.begin(), hist_grey.end(), 0.0);

    for (unsigned int i = 0; i < hist.size(); i++)
    {
      grey.at(i) = grey.at(i) - histmean;
      hist_var.at(i) = hist.at(i)* pow(grey.at(i), 2);
      hist_sk.at(i) = hist.at(i)* pow(grey.at(i), 3);
      hist_kur.at(i) = hist.at(i)* pow(grey.at(i), 4);
    }
    double variance = std::accumulate(hist_var.begin(), hist_var.end(), 0.0);
    //std::cout << std::accumulate(hist_sk.begin(), hist_sk.end(), 0.0) << ",,,," << std::accumulate(hist_kur.begin(), hist_kur.end(), 0.0);
    //std::cout << pow(sqrt(variance), 3) << ",,,," << pow(sqrt(variance), 4);
    double skewness = std::accumulate(hist_sk.begin(), hist_sk.end(), 0.0) / pow(sqrt(variance), 3);
    double kurtosis = std::accumulate(hist_kur.begin(), hist_kur.end(), 0.0) / pow(sqrt(variance), 4);
    std::set<int> x_value, y_value, z_value;
    for (unsigned int i = 0; i < index_vec.size(); i++)
    {
      x_value.insert(index_vec[i][0]);
      y_value.insert(index_vec[i][1]);
      z_value.insert(index_vec[i][2]);
    }


    std::vector < double> intensity_feature;
    std::vector <std::string> intensity_featurename;
    intensity_featurename.push_back("volume"); intensity_feature.push_back(volume);
    intensity_featurename.push_back("min_int"); intensity_feature.push_back(min);
    intensity_featurename.push_back("max_int"); intensity_feature.push_back(max);
    intensity_featurename.push_back("mean_int"); intensity_feature.push_back(mean);
    intensity_featurename.push_back("variance"); intensity_feature.push_back(variance);
    intensity_featurename.push_back("skewness"); intensity_feature.push_back(skewness);
    intensity_featurename.push_back("kurtosis"); intensity_feature.push_back(kurtosis);

    featurename->insert(featurename->end(), intensity_featurename.begin(), intensity_featurename.end());
    featurevec->insert(featurevec->end(), intensity_feature.begin(), intensity_feature.end());


    //set region of interest in the input image where feature extraction will take place

    ImageType::Pointer  newroi = ImageType::New();
    ImageType::RegionType region;

    const ImageType::SizeType  size = { { x_value.size(), y_value.size(), z_value.size() } };
    const ImageType::IndexType start = { { *x_value.begin(), *y_value.begin(), *z_value.begin() } };

    region.SetSize(size);
    //  region.SetIndex(start);
    newroi->SetRegions(region);
    newroi->Allocate();
    newroi->FillBuffer(0.0);
    //intensity features
    // typedef itk::ImageRegionIteratorWithIndex< MaskImageType > IteratorType1;
    const ImageType::SizeType image_size = newmask->GetLargestPossibleRegion().GetSize();

    int x1 = 0; int y1 = 0; int z1 = 0;
    for (unsigned int z = start[2]; z < start[2] + size[2]; z1++, z++)
    {
      y1 = 0;
      for (unsigned int y = start[1]; y < start[1] + size[1]; y1++, y++)
      {
        x1 = 0;
        for (unsigned int x = start[0]; x < start[0] + size[0]; x1++, x++)
        {
          ImageType::IndexType ind1 = { x, y, z };
          ImageType::IndexType ind2 = { x1, y1, z1 };
          if (x < image_size[0] && y < image_size[1] && z < image_size[2])
          {
            newroi->SetPixel(ind2, newmask->GetPixel(ind1)[0]);
          }
          else newroi->SetPixel(ind2, 0);
        }
      }
    }


    NeighborhoodType neighborhood;
    neighborhood.SetRadius(1);
    unsigned int centerIndex = neighborhood.GetCenterNeighborhoodIndex();
    OffsetType offset;
    OffsetVector::Pointer offsets = OffsetVector::New();
    for (unsigned int d = 0; d < centerIndex; d++)
    {
      offset = neighborhood.GetOffset(d);
      offsets->push_back(offset);
    }
    generateTextureFeatures texturefeature;
    std::vector<std::string> glcm_featurename;
    std::vector<double> glcm_feature;

    texturefeature.calculateTextureFeatures<ImageType, OffsetVector::Pointer>(PETimg_windowed, newroi, min, max, offsets, glcm_featurename, glcm_feature);
    featurename->insert(featurename->end(), glcm_featurename.begin(), glcm_featurename.end());
    featurevec->insert(featurevec->end(), glcm_feature.begin(), glcm_feature.end());
    glcm_featurename.clear();
    glcm_feature.clear();
    texturefeature.calculateRunLength<ImageType, OffsetVector::Pointer>(PETimg_windowed, newroi, min, max, offsets, glcm_featurename, glcm_feature);
    featurename->insert(featurename->end(), glcm_featurename.begin(), glcm_featurename.end());
    featurevec->insert(featurevec->end(), glcm_feature.begin(), glcm_feature.end());
    glcm_featurename.clear();
    glcm_feature.clear();
    texturefeature.ShapeFeatures<ImageType>(newroi, glcm_featurename, glcm_feature);
    featurename->insert(featurename->end(), glcm_featurename.begin(), glcm_featurename.end());
    featurevec->insert(featurevec->end(), glcm_feature.begin(), glcm_feature.end());
    glcm_featurename.clear();
    glcm_feature.clear();


    int radius = 1; int neigh = 8;
    //Pad image before calculating LBP
    typedef itk::ConstantPadImageFilter < ImageType, MaskImageType >
      ConstantPadImageFilterType;

    ImageType::SizeType lowerExtendRegion;
    lowerExtendRegion[0] = radius;
    lowerExtendRegion[1] = radius;
    lowerExtendRegion[2] = radius;

    ImageType::SizeType upperExtendRegion;
    upperExtendRegion[0] = radius;
    upperExtendRegion[1] = radius;
    upperExtendRegion[2] = radius;

    ImageType::PixelType constantPixel = 0;

    ConstantPadImageFilterType::Pointer padFilter
      = ConstantPadImageFilterType::New();
    padFilter->SetInput(newroi);
    //padFilter->SetPadBound(outputRegion); // Calls SetPadLowerBound(region) and SetPadUpperBound(region)
    padFilter->SetPadLowerBound(lowerExtendRegion);
    padFilter->SetPadUpperBound(upperExtendRegion);
    padFilter->SetConstant(constantPixel);
    padFilter->Update();
    MaskImageType::Pointer  lbproi = padFilter->GetOutput();
    x1 = -1;  y1 = -1; z1 = -1;
    for (unsigned int z = start[2] - radius; z < start[2] + size[2] + radius; z1++, z++)
    {
      y1 = -1;
      for (unsigned int y = start[1] - radius; y < start[1] + radius + size[1]; y1++, y++)
      {
        x1 = -1;
        for (unsigned int x = start[0] - radius; x < start[0] + radius + size[0]; x1++, x++)
        {
          ImageType::IndexType ind1 = { x, y, z };
          ImageType::IndexType ind2 = { x1, y1, z1 };
          MaskImageType::PixelType   pixelValue;
          if (x < image_size[0] && y < image_size[1] && z < image_size[2])
          {
            pixelValue[0] = newmask->GetPixel(ind1)[0];   // x component
            pixelValue[1] = newmask->GetPixel(ind1)[1];
            lbproi->SetPixel(ind2, pixelValue);
          }
          else
          {
            pixelValue[0] = 0;  // x component
            pixelValue[1] = 0;
            lbproi->SetPixel(ind2, pixelValue);
          }
        }
      }
    }
    LBPFeatures lbpfeatures;

    lbpfeatures.calculateLBP<float, 3>(newmask, maskin, lbproi, radius, neigh, &glcm_featurename, &glcm_feature);
    //  Lbpfeatures(newmask, maskin, lbproi, radius, neigh);
    featurename->insert(featurename->end(), glcm_featurename.begin(), glcm_featurename.end());
    featurevec->insert(featurevec->end(), glcm_feature.begin(), glcm_feature.end());
  }
  catch (std::exception ex)
  {
    cout << std::string("Error in getFeature:") + ex.what() << std::endl;
  }
  catch (...)
  {
    cout << "Unknown Error in getFeature:" << std::endl;
  }

}

std::vector<float> getLables(const string filename)
{
  std::vector<float> labels;

  ifstream file(filename); // declare file stream: 
  std::string value;
  while (file.good())
  {
    getline(file, value);
    labels.push_back(atof(value.c_str()));
  }
  return labels;
}


int main(int argc, char** argv)
{
  try
  {
    cbica::CmdParser parser(argc, argv);


    string CTimage, PETimage, mask, modelfilename, outputfile;
    InputImagePointerType imgin = ImageType::New();
    InputImagePointerType maskin = ImageType::New();
    InputImagePointerType PET_imgin = ImageType::New();
    string filename1;
    int Number_features;

    parser.addRequiredParameter("c", "CTimage", cbica::Parameter::FILE, "NIfTI", "Input CT Image");
    parser.addRequiredParameter("p", "PETimage", cbica::Parameter::FILE, "NIfTI", "Input PET Image");

    parser.addRequiredParameter("m", "mask", cbica::Parameter::FILE, "NIfTI", "Input mask image");
    parser.addRequiredParameter("t", "model", cbica::Parameter::FILE, ".xml", "Trained SVM model");
    parser.addRequiredParameter("o", "outputfile", cbica::Parameter::FILE, ".txt", "output file");
    parser.addOptionalParameter("L", "Logger", cbica::Parameter::FILE, "log file which user has write access to", "Full path to log file to store console outputs", "By default, only console output is generated");


    std::string logFile, inputPatterns, labelfile, dataDir, flag;
    cbica::Logging logger;

    if (parser.isPresent("L"))
    {
      parser.getParameterValue("L", logFile);
      logger.UseNewFile(logFile);
    }

    if (parser.isPresent("c"))
    {
      parser.getParameterValue("c", CTimage);

    }
    if (parser.isPresent("p"))
    {
      parser.getParameterValue("p", PETimage);
    }
    if (parser.isPresent("m"))
    {
      parser.getParameterValue("m", mask);
    }
    if (parser.isPresent("t"))
    {
      parser.getParameterValue("t", modelfilename);
    }
    if (parser.isPresent("o"))
    {
      parser.getParameterValue("o", outputfile);
    }

    if (!cbica::fileExists(CTimage))
    {
      logger.WriteError("Aborting:   Could not read " + CTimage);
      return EXIT_FAILURE;
    }
    else
    {
      imgin = cbica::ReadImage(CTimage);
    }

    //latest working version
    if (!cbica::fileExists(PETimage))
    {
      logger.WriteError("Aborting:   Could not read " + PETimage);
      return EXIT_FAILURE;
    }
    else
    {
      PET_imgin = cbica::ReadImage(PETimage);
    }

    if (!cbica::fileExists(mask))
    {
      logger.WriteError("Aborting:   Could not read " + mask);
      return EXIT_FAILURE;
    }
    else
    {
      maskin = cbica::ReadImage(mask);

    }

    std::vector<std::string> featurename;
    std::vector<float> featurevec;

    getFeature<float, 3>(imgin, PET_imgin, maskin, &featurename, &featurevec);

    Number_features = featurename.size();
    cv::Mat Feature_matrix;

    Feature_matrix.create(1, Number_features, CV_32F);

    for (int j = 0; j < Number_features; j++)
    {

      Feature_matrix.at<float>(0, j) = featurevec.at(j);
    }

    auto svm = cv::ml::SVM::create();
    cv::Mat data = cv::Mat(1, Number_features + 1, CV_32FC1);
    ofstream myfile;
    myfile.open(outputfile, std::fstream::app);

    //  read_data_from_csv(dataDir + "scaled_training_features.csv", data, 1);
   // svm = cv::ml::SVM::load(modelfilename);


    int result = svm->predict(Feature_matrix);
    float confidence = 1.0 / (1.0 + exp(-result));//TBD

    if (myfile.is_open())
    {
      if (confidence > 0.6)
        myfile << confidence << endl;
      else
        myfile << confidence << endl;
      myfile.close();
    }

    return EXIT_SUCCESS;

  }
  catch (std::exception ex)
  {
    cout << std::string("Error in getFeature:") + ex.what() << std::endl;
  }
  catch (...)
  {
    cout << "Unknown Error in getFeature:" << std::endl;
  }
}



