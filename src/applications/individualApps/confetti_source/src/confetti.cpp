#include "confetti.h"
//Genrate connectivity
int Confetti::genConnectivitySig( string tdiCsvFile, string signatureFileName, string fiberFile)
{
  Mat mP2V;
  Point3f origin;
  Utils::resetTimer();

  Utils::logText("Image Loading : " + tdiCsvFile);
  vector<ImageType::Pointer> Images = readImages(tdiCsvFile, mP2V, origin);
  if (Images.empty())
  {
    Utils::logText("Error. Cannot read images from:" + tdiCsvFile);
    return EXIT_FAILURE;
  }
  //Check image for size integrity
  Point3i imageSize;
  imageSize.x = Images[0]->GetLargestPossibleRegion().GetSize()[0] - 1;
  imageSize.y = Images[0]->GetLargestPossibleRegion().GetSize()[1] - 1;
  imageSize.z = Images[0]->GetLargestPossibleRegion().GetSize()[2] - 1;
  for (size_t i = 1; i < Images.size(); i++)
  {
    Point3i curSize;
    curSize.x = Images[i]->GetLargestPossibleRegion().GetSize()[0] - 1;
    curSize.y = Images[i]->GetLargestPossibleRegion().GetSize()[1] - 1;
    curSize.z = Images[i]->GetLargestPossibleRegion().GetSize()[2] - 1;
    if (curSize.x != imageSize.x || curSize.y != imageSize.y || curSize.z != imageSize.z)
    {
      stringstream ss;
      ss << imageSize << " vs" << curSize << "At Image" << i;
      Utils::logText("Error.Exiting as images size does not match:" + ss.str());
      return EXIT_FAILURE;
    }
  }
  //For limit checking
  imageSize.x = imageSize.x - 1;
  imageSize.y = imageSize.y - 1;
  imageSize.z = imageSize.z - 1;

  Utils::logText("Fiber Loading: " + fiberFile);
  vector<long> fiberIndexs;
  vector<Point3f> fiberXYZSet = Utils::readFiberBfloat(fiberFile, fiberIndexs);
  if (fiberXYZSet.empty() || fiberIndexs.empty())
  {
    Utils::logText("Error. Cannot read:" + fiberFile);
    return EXIT_FAILURE;
  }

  Utils::logText("Affine Transform..");


  vector<Point3i> fiberIJK = affineTransform(fiberXYZSet, mP2V, origin, imageSize);
  Point3i* currentPtr = &fiberIJK.front();
  Utils::logText("Generating Connectivity..");
  vector<vector< float > >  sigs;
  for (size_t idx = 0; idx < fiberIndexs.size(); idx++)
  {
    vector<float> sig = connectivityCoreFast(currentPtr, fiberIndexs[idx], Images);
    sigs.push_back(sig);
    currentPtr = currentPtr + fiberIndexs[idx];
    Utils::progressUpdate(idx, fiberIndexs.size(), "Connectivity");
  }
  Utils::logText("Saving output file " + signatureFileName);
  Utils::writeCsv(sigs, signatureFileName);
  Utils::logText("Finished Connectivity Signature generation!");
  return EXIT_SUCCESS;
}
//Extract Trcks
int Confetti::extractTract(string outDIR, string tempRoot, string inputIDFile, string inputNamesFile, string inputLabelFile, string inputBfloatFile)
{
  Utils::resetTimer();
  Utils::logText("Extracting tracks for:" + inputIDFile);
  vector< vector<int> > bundles = Utils::readCsv<int>(inputIDFile);
  if (bundles.empty())
  {
    Utils::logText("Error reading:" + inputIDFile);
    return EXIT_FAILURE;
  }

  Utils::logText("Reading:" + inputNamesFile);
  vector<string> bundleNames;
  ifstream fin(inputNamesFile.c_str());
  if (fin.is_open())
  {
    for (string line; getline(fin, line);)
    {
      bundleNames.push_back(line);
      line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    }
    fin.close();
  }
  if (bundleNames.empty())
  {
    Utils::logText("Error reading:" + inputNamesFile);
    return EXIT_FAILURE;
  }

  Utils::logText("Reading:" + inputLabelFile);
  vector< vector<int> > vecClusterIDs = Utils::readCsv<int>(inputLabelFile);
  if (vecClusterIDs.empty())
  {
    Utils::logText("Error reading:" + inputLabelFile);
    return EXIT_FAILURE;
  }

  vector<int> clusterIDs = Utils::flatten(vecClusterIDs);
  vector<string> outputFileNames;
  for (size_t bunId = 0; bunId < bundleNames.size(); bunId++)
  {
    outputFileNames.push_back(outDIR + bundleNames[bunId] + ".Bfloat");
  }

  extractTractCore(inputBfloatFile, clusterIDs, bundles, outputFileNames);
  Utils::logText("Finished running extraction");
  return EXIT_SUCCESS;
}
vector<ImageType::Pointer> Confetti::readImages( const string& tdiCsvFile, Mat&mP2V, Point3f& origin)
{
  vector<ImageType::Pointer> images;
  vector<string> tdiFileNames = Utils::flatten<string>(Utils::readCsv<string>(tdiCsvFile));
  if (tdiFileNames.empty())
  {
    Utils::logText("Error: TDI file list is empty :"+ tdiCsvFile);
    return images;
  }
  string basePath = string();
  if (!Utils::fileExist(tdiFileNames[0]))
  {
    basePath + Utils::getBasePath(tdiCsvFile);
  }
  Utils::logText("Reading TDI files (" + Utils::toString(tdiFileNames.size()) + ") ...");
  for (size_t i = 0; i < tdiFileNames.size(); i++)
  {
    string niftiImageFile = basePath+tdiFileNames[i];
    if (!Utils::fileExist(niftiImageFile))
    {
      return images;
    }
    if (i == 0)
    {
      getAffine(niftiImageFile, mP2V, origin);
    }
    images.push_back(Utils::readNifti(niftiImageFile));
    Utils::logText("Loading: " + niftiImageFile);
    Utils::progressUpdate(i, tdiFileNames.size(), "Reading Images");
  }
  return images;
}

vector<Point3i> Confetti::affineTransform(vector<Point3f> &fiberXYZ, const Mat &mP2V, const Point3f &origin, const Point3i &maxSize)
{
  vector<Point3i> fiberIJK(fiberXYZ.size());
  for (size_t i = 0; i < fiberXYZ.size(); i++)
  {
    //Dot product and ciel
    fiberIJK[i] = Point3f(Mat(mP2V*Mat(fiberXYZ[i] - origin))) + Point3f(0.5, 0.5, 0.5);
    fiberIJK[i].x = min(fiberIJK[i].x, maxSize.x);
    fiberIJK[i].y = min(fiberIJK[i].y, maxSize.y);
    fiberIJK[i].z = min(fiberIJK[i].z, maxSize.z);
    fiberIJK[i].x = max(fiberIJK[i].x, 0);
    fiberIJK[i].y = max(fiberIJK[i].y, 0);
    fiberIJK[i].z = max(fiberIJK[i].z, 0);
    Utils::progressUpdate(i, fiberXYZ.size(), "Affine Transform");
  }
  return fiberIJK;
}
vector<float> Confetti::connectivityCoreFast(const Point3i*  _fiberIJK, const size_t _size, const vector<ImageType::Pointer> &Images)
{
  vector<Point3i> fiberIJK;
  fiberIJK.assign(_fiberIJK, _fiberIJK + _size);
  //Remove duplicates near by ( only works if sorted otherwise )
  fiberIJK.erase(unique(fiberIJK.begin(), fiberIJK.end()), fiberIJK.end());

  Mat Pf = Mat::zeros(fiberIJK.size(), Images.size(), CV_32F);
  for (size_t v = 0; v < fiberIJK.size(); v++)
  {
    for (size_t r = 0; r < Images.size(); r++)
    {
      ImageType::IndexType index;
      index[0] = fiberIJK[v].x;
      index[1] = fiberIJK[v].y;
      index[2] = fiberIJK[v].z;
      Pf.at<float>(v, r) = Images[r]->GetPixel(index);
    }
  }
  const double sigma = 5.0;
  //apply smoothing
  Mat Lf;
  cv::GaussianBlur(Pf, Lf, Size((8 * sigma + 1), 1), sigma);
  const int lw = int(4.0 * sigma + 0.5);
  const int kernalSize = (2 * lw + 1);
  cv::GaussianBlur(Pf, Lf, Size(1, kernalSize), sigma);

  //weighted average of voxels
  vector<float> ww = createBellCurve(fiberIJK.size());
  vector<float> F = applyBellCurve(Lf, ww);
  for (size_t i = 0; i < F.size(); i++)
  {
    F[i] = roundf(F[i]);
  }

  return F;
}
vector<float> Confetti::getAffine(const string& imgFileName, Mat&mP2V, Point3f& origin)
{
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(imgFileName.c_str(), itk::ImageIOFactory::ReadMode);
  imageIO->SetFileName(imgFileName);
  imageIO->ReadImageInformation();

  vector<float> vecOrigin, spacing, dims;
  unsigned int numDim = imageIO->GetNumberOfDimensions();
  for (unsigned i = 0; i <numDim; i++)
  {
    vecOrigin.push_back(imageIO->GetOrigin(i));
    spacing.push_back(imageIO->GetSpacing(i));
    dims.push_back(imageIO->GetDimensions(i));
  }

  //Reading http://public.kitware.com/pipermail/insight-users/2010-June/037400.html
  //LPS TO RAS https://itk.org/Wiki/ITK_Release_4.0_Orientation_by_Torsten
  //"LPS" to "RAS" means x becomes-x, y becomes - y, and z stays z.
  //"ARS" to "RAS" means that x and y elements are swapped in all direction vectors and also in the space origin.
  Point3f RasMask = Point3f(-1.0, -1.0, 1.0);//TBD

  Mat S = Mat::zeros(3, 3, CV_32FC1);
  S.at<float>(0, 0) = spacing[0] * RasMask.x;
  S.at<float>(1, 1) = spacing[1] * RasMask.y;
  S.at<float>(2, 2) = spacing[2] * RasMask.z;

  mP2V = 1.0 / S;
  origin = Point3f(vecOrigin[0] * RasMask.x, vecOrigin[1] * RasMask.y, vecOrigin[2] * RasMask.z);
  return dims;
}
vector<float> Confetti::createBellCurve(size_t size)
{
  vector<float> ss;
  float dvSr = float(size - 1.);
  for (size_t i = 0; i < size; ++i)
  {
    ss.push_back(i / dvSr);
  }
  ss[0] = 0.01f;
  ss[size - 1] = 1.f - 0.01f;
  for (size_t i = 0; i < ss.size(); ++i)
  {
    ss[i] = pow(ss[i], -0.5)*pow(1 - ss[i], -0.5) + 0.01f;
  }
  float sumF = cv::sum(ss).val[0];
  for (size_t i = 0; i < ss.size(); ++i)
  {
    ss[i] = ss[i] / sumF;
  }
  return ss;
}
vector<float> Confetti::applyBellCurve(Mat& data, const vector<float>& curve)
{
  vector<float> res;
  assert(uint(data.rows) == curve.size());
  for (int r = 0; r < data.rows; r++)//rows=28
  {
    for (int c = 0; c < data.cols; c++)//cols =87
    {
      data.at<float>(r, c) = data.at<float>(r, c)*curve[r];
    }
  }
  Mat dataReduced;
  cv::reduce(data, dataReduced, 0, CV_REDUCE_SUM, CV_32F);
  dataReduced.row(0).copyTo(res);
  return res;
}
void Confetti::extractTractCore(string bfloatFile, vector<int> clusterIDs, vector< vector<int> > selectedTracts, vector<string> outputFileNames)
{
  ifstream f(bfloatFile.c_str(), ios::binary);
  if (!f.is_open())
  {
    Utils::logText("Error. Cannot Read: " + bfloatFile);
    return;
  }
  //output files
  vector<ofstream*> tractFiles;
  for (size_t fileId = 0; fileId < outputFileNames.size(); fileId++)
  {
    tractFiles.push_back(new ofstream(outputFileNames[fileId].c_str(), ios::binary));
    if (!(*tractFiles.back()).is_open())
    {
      Utils::logText("Error. Cannot create: " + outputFileNames[fileId]);
      return;
    }
  }

  Utils::logText("Extracting clusters.. ");
  for (size_t fib = 0; fib < clusterIDs.size(); fib++)
  {
    Utils::progressUpdate(fib, clusterIDs.size(), "Extract clusters");
    int numPoints = Utils::readVal(f);
    int seedPoint = Utils::readVal(f);
    //read physical coordinates
    vector<Point3f> coords(numPoints);
    for (size_t i = 0; i < coords.size(); i++)
    {
      float x = Utils::readVal(f);
      float y = Utils::readVal(f);
      float z = Utils::readVal(f);

      coords[i] = Point3f(x, y, z);
    }

    //check if we want this fiber and then write it to the right file
    for (size_t tr = 0; tr< selectedTracts.size(); tr++)
    {
      if (find(selectedTracts[tr].begin(), selectedTracts[tr].end(), clusterIDs[fib]) != selectedTracts[tr].end())
      {
        Utils::writeVal(*tractFiles[tr], numPoints);
        Utils::writeVal(*tractFiles[tr], seedPoint);
        for (size_t p = 0; p<(uint)numPoints; p++)
        {
          Utils::writeVal(*tractFiles[tr], coords[p].x);
          Utils::writeVal(*tractFiles[tr], coords[p].y);
          Utils::writeVal(*tractFiles[tr], coords[p].z);
        }
        break;
      }
    }
  }
  f.close();
  for (size_t fileId = 0; fileId < tractFiles.size(); fileId++)
  {
    tractFiles[fileId]->close();
  }

  return;
}
