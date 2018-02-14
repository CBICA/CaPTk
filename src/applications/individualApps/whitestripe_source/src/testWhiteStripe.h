#ifndef TEST_WHITE_STRIPE_H
#define TEST_WHITE_STRIPE_H
#include "whiteStripe.h"
#include "utils.h"
#include <ctime>
class TestWhiteStripe
{

  static void getMetaInfo(string filename, vector<float> &mids, vector<float> &counts, vector<float> &fit, float& mode)
  {
    vector< vector<float> > values;
    ifstream fin(filename.c_str());
    for (string line; getline(fin, line);)
    {
      while (line.find("x") == string::npos) {
        getline(fin, line);
      }// skip 
      while (true)
      {
        getline(fin, line);
        if (line.find("modeQ") != string::npos)
        {
          break;
        }
        replace(line.begin(), line.end(), ',', ' ');
        stringstream ss(line);
        float val;
        ss >> val;
        mids.push_back(val);
        ss >> val;
        counts.push_back(val);
        ss >> val;
        fit.push_back(val);
      }
      getline(fin, line);//Mode q val
      getline(fin, line);//Mode heaser
      getline(fin, line);//actual val;
      stringstream ss(line);
      ss >> mode;
      break;
    }
    return;
  }

public:
  static int testImgMean()
  {
    string fileName = "../../data/image_1_T1.nii.gz";
    itk::SmartPointer<ImageType> img = Utils::readNifti(fileName);
    if (img.IsNull())
    {
      cout << "Cannot read - " << fileName << endl;
      return false;
    }
    double meanGT = 15.384206743161940;
    double mean = Utils::getMeanValue(img);
    if (mean == meanGT)
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }  
  static int testQuantile()
  {
    vector<float> testVals = { 16, 17, 18, 19, 20, 21, 22, 23, 24 };
    vector<float> probs = { 0.5f, 0.6f, 0.1f };
    vector<float> res = Utils::quantileFast(testVals, probs);
    vector<float> res2 = Utils::quantile(testVals, probs);
    if (res[0] == 20.f && res[1] == 20.f && res[2] == 16.f)
    {
      if (res2[0] == 20.f && res2[1] == 20.f && res2[2] == 16.f)
      {
        return EXIT_SUCCESS;
      }
    }
    return EXIT_FAILURE;
  }
 
  static int testWhiteStripe()
  {
    string fileName = "../../data/image_1_T1.nii.gz";
    string fileNameRes = "../../data/image_1_T1_Res.nii.gz";
    string fileNameMask = "../../data/image_1_T1_mask.nii.gz";

    itk::SmartPointer<ImageType> img = Utils::readNifti(fileName);
    itk::SmartPointer<ImageType> imgResGT = Utils::readNifti(fileNameRes);
    itk::SmartPointer<ImageType> imgMaskGT = Utils::readNifti(fileNameMask);
    if (img.IsNull() || imgResGT.IsNull() || imgMaskGT.IsNull())
    {
      cout << "Cannot read - " << fileName << endl;
      return EXIT_FAILURE;
    }
    itk::SmartPointer<ImageType> imgRes, imgMask;
    WhiteStripe obj;
    imgRes = obj.process(img, imgMask);
    bool bRes = Utils::pixelMatchImages(imgRes, imgResGT);
    bool bMask = Utils::pixelMatchImages(imgMask, imgMaskGT);
    bool bBs = Utils::pixelMatchImages(img, imgMaskGT);
    if (bRes&&bMask&&!bBs)
    {
      return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
  }
  static void test(int(*testFun)(), string name = "")
  {
    clock_t  start=std::clock();
    string resString = "\t Failed";
    if (testFun() == EXIT_SUCCESS)
    {
      resString = "\t Passed";
    }
    double millSec = (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
    cout << "Test:" << name + resString  << "\t Time:" << millSec << "ms" << std::endl;
  }
  static void runAllTests()
  {

     test(&testQuantile, "Quantile");
     test(&testImgMean, "ImgMean");
     test(&testWhiteStripe, "WhiteStripe");
     
  }
};
#endif
