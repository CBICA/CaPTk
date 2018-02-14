#ifndef _CONFETTI_H
#define _CONFETTI_H

/**
\file  confetti.h

\brief Implementation of the Confetti Algorithm

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2017 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html
Author: Ratheesh Kalarot
*/
#include "dtiUtils.h"
using namespace cv;
/**
\class Confetti
\brief Implementation of  Confetti algorithm
Reference:
Tunç B, Smith AR, Wasserman D, et al.
Multinomial Probabilistic Fiber Representation for Connectivity Driven Clustering. 
Information processing in medical imaging : proceedings of the . conference.
2013;23:730-741.
*/
class Confetti {

public:

	//--------------------------Core----------------------------
  //Listner for log and progress updates
	void setListner(TXT_CALLBACK_TYPE listner)
	{
		Utils::setListner(listner);
		return;
	}
  //Genrate connectivity
  int genConnectivitySig(string tdiCsvFile, string signatureFileName, string fiberFile);
  //Extract Trcks
  int extractTract(string outDIR, string tempRoot, string inputIDFile, string inputNamesFile, string inputLabelFile, string inputBfloatFile);
	
  int createTDI(string fiberFile, string niftiImageFile)
  {
    Utils::logText("Fiber Loading: " + fiberFile);
    vector<long> fiberIndexs;
    vector<Point3f> fiberXYZSet = Utils::readFiberBfloat(fiberFile, fiberIndexs);
    if (fiberXYZSet.empty() || fiberIndexs.empty())
    {
      Utils::logText("Error. Cannot read:" + fiberFile);
      return EXIT_FAILURE;
    }

    Utils::logText("Affine Transform..");
    Mat mP2V;
    Point3f origin;
    vector<float> imgSize= getAffine(niftiImageFile, mP2V, origin);
    if (imgSize.size()<3)
    {
      Utils::logText("Error. Cannot read size of the file:" + niftiImageFile);
      return EXIT_FAILURE;
    }
    Point3i imageSize(imgSize[0] - 1, imgSize[1] - 1, imgSize[2] - 1);
    vector<Point3i> fiberIJK = affineTransform(fiberXYZSet, mP2V, origin, imageSize);
    ImageType::IndexType idx;
    ImageType::PixelType pixVal;
    ImageType::Pointer img = Utils::readNifti(niftiImageFile);
    itk::ImageRegionIterator<ImageType> imageIterator(img, img->GetLargestPossibleRegion());
    imageIterator.GoToBegin();
    while (!imageIterator.IsAtEnd())
    {
      imageIterator.Set(0);
      ++imageIterator;
    }
    for (size_t i = 0; i < fiberIJK.size(); i++)
    {
      idx[0] = fiberIJK[i].x;
      idx[1] = fiberIJK[i].y;
      idx[2] = fiberIJK[i].z;
      pixVal = img->GetPixel(idx) + 1;
      img->SetPixel(idx, pixVal);
    }
  }
private :
  // Read probtrax images from data directory
  vector<ImageType::Pointer> readImages(const string& tdiCsvFile, Mat&mP2V, Point3f& origin);
  //Affine Transform
  vector<Point3i> affineTransform(vector<Point3f> &fiberXYZ, const Mat &mP2V, const Point3f &origin, const Point3i &maxSize);
  //Core of connectiviy 
  vector<float> connectivityCoreFast(const Point3i*  _fiberIJK, const size_t _size, const vector<ImageType::Pointer> &Images);
  vector<float> getAffine(const string& imgFileName, Mat&mP2V, Point3f& origin);
  //Create Bellcurve ofspecified lenth
  vector<float> createBellCurve(size_t size);
  vector<float> applyBellCurve(Mat& data, const vector<float>& curve);
  //Core of extracting tracks
  void extractTractCore(string bfloatFile, vector<int> clusterIDs, vector< vector<int> > selectedTracts, vector<string> outputFileNames);
};

#endif
