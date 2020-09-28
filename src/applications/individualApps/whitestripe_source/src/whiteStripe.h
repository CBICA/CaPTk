/**
\file  whiteStripe.h

\brief Implementation of the WhiteStripe Algorithm

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2017 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html
Author: Ratheesh Kalarot
*/

#ifndef WHITE_STRIPE_H
#define WHITE_STRIPE_H

#define WHITESTRIPE_VERSION "1.0.0"

#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <iomanip>

#include "itkNiftiImageIO.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkExtractImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkContinuousIndex.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageDuplicator.h"
#include "KernelFit.h"
#include "utils.h"

using namespace std;
using namespace cv;

typedef float PixelType;
typedef itk::Image<PixelType, 3> ImageType;
/**
\Struct HistSruct
\brief structure for basic histogram operations
*/
struct HistSruct
{
  vector<float> mids;
  vector<float> counts;
  vector<float> countsOrig;
  int m_tissuesMax;
  float m_smoothMax;
  float m_smoothDelta;
  int m_histSize;
  HistSruct()
  {
    m_tissuesMax = 5;
    m_smoothMax = 10.0;
    m_smoothDelta = 0.5;
    m_histSize = 2000;
    mids.clear();
    counts.clear();
  }
  void setParams(int tissuesMax, float smoothMax, float smoothDelta, int histSize)
  {
    m_tissuesMax = tissuesMax;
    m_smoothMax = smoothMax;
    m_smoothDelta = smoothDelta;
    m_histSize = histSize;
  }
  void setVals(vector<float> _mids, vector<float> _counts)
  {
    mids = _mids;
    counts = _counts;
    countsOrig = counts;
  }
  void compute(vector<float> imagePts)
  {
    mids.clear();
    counts.clear();
    countsOrig.clear();


    float max = *max_element(imagePts.begin(), imagePts.end());
    float min = *min_element(imagePts.begin(), imagePts.end());

    float range[] = { min, max };
    const float* histRange = { range };
    bool uniform = true; bool accumulate = false;
    Mat histCounts;
    vector<Mat> images;
    images.push_back(Mat(imagePts));

    calcHist(&images[0], images.size(), 0, Mat(), histCounts, 1, &m_histSize, &histRange, uniform, accumulate);
    counts = histCounts;
    mids.resize(counts.size());
    float delMid = (max - min) / float(counts.size());
    for (size_t i = 0; i < counts.size(); i++)
    {
      if (i == 0)
      {
        mids[i] = min + delMid / 2;
      }
      else
      {
        mids[i] = mids[i - 1] + delMid;
      }
    }
    countsOrig = counts;
  }

  vector<int> getPeaks(vector<float>& peakMidsOut, float delMin = 0.01, float valMin = 0.25)
  {
    vector<int> peaks;
    peakMidsOut.clear();
    for (size_t i = 1; i < counts.size() - 1; i++)
    {
      if ((counts[i - 1] < counts[i]) && (counts[i + 1] < counts[i]))
      {
        float delta = fabsf(counts[i] - counts[i - 1]) + fabsf(counts[i] - counts[i + 1]);
        if (delta>delMin)
        {
          //TBD calc Max only once
          float max = *max_element(counts.begin(), counts.end());
          if (counts[i] / max > valMin)
          {
            peaks.push_back(i);
            peakMidsOut.push_back(mids[i]);
          }
        }
      }
    }
    return peaks;
  }
  void smooth(float bandWidth = 3.0)
  {
    if (counts.empty()) return;
    KernelFit1D<float> kernelFun(mids, counts, bandWidth);
    counts = kernelFun.Solve(mids);
  }
  Mat getHisDrawing(Mat canvas, Scalar color = Scalar(255, 255, 0), int yOffSet = 0)
  {
    if (counts.empty() || canvas.empty()) return Mat();

    int width = canvas.cols;
    int height = canvas.rows;
    Mat hist = Mat(counts).clone();
    double binWidth = (width *1.0) / hist.rows;
    normalize(hist, hist, 0, canvas.rows, NORM_MINMAX, -1, Mat());
    for (size_t i = 1; i < (uint)width; i++)
    {
      int idx = min(int(i / binWidth + 0.5), hist.rows - 1);
      int idxLast = min(int((i - 1) / binWidth + 0.5), hist.rows - 1);
      Point firstPt(i, height + yOffSet - cvRound(hist.at<float>(idx)));
      Point lastPt(i, height + yOffSet - cvRound(hist.at<float>(idxLast)));
      line(canvas, firstPt, lastPt, color);
    }
    return canvas;
  }
  void drawPeaks(Mat histImage, vector<int>& peaks, Scalar color = Scalar(0, 0, 255), int thickness = 1)
  {
    if (counts.empty()) return;

    double biWidth = (histImage.cols*1.0) / counts.size();
    for (size_t i = 0; i < peaks.size(); i++)
    {
      Point p1(biWidth * peaks[i], histImage.rows);
      Point p2(biWidth * peaks[i], 0);
      line(histImage, p1, p2, color, thickness);
    }
    return;
  }
  float  drawStdevPeak(Mat histImage, vector<float> histOrig, int peak, Scalar color = Scalar(0, 255, 255 / 2), int thickness = 2)
  {
    if (histOrig.empty()) return 0.0;
    int width = 0.05* histOrig.size();
    vector<float> histGaauss;

    for (size_t i = 0; i < (uint)width; i++)
    {
      if ((peak - width) < 0 || (peak + width) >= (int)histOrig.size())
      {
        width = i - 1;
        break;
      }
      histGaauss.push_back(histOrig[peak - i]);
      histGaauss.push_back(histOrig[peak + i]);
    }
    double sum = std::accumulate(histGaauss.begin(), histGaauss.end(), 0.0);
    double mean = sum / histGaauss.size();
    double sqSum = std::inner_product(histGaauss.begin(), histGaauss.end(), histGaauss.begin(), 0.0);
    double stdev = std::sqrt(sqSum / histGaauss.size() - mean * mean);


    vector<float> peaks;
    peaks.push_back(peak - width);
    peaks.push_back(peak + width);

    double biWidth = (histImage.cols*1.0) / counts.size();
    for (size_t i = 0; i < peaks.size(); i++)
    {
      Point p1(biWidth * peaks[i], histImage.rows);
      Point p2(biWidth * peaks[i], 0);
      line(histImage, p1, p2, color, thickness);
    }
    std::ostringstream stdevStr;
    stdevStr << std::setprecision(3) << stdev;
    string text = "STD=" + stdevStr.str();
    putText(histImage, text, Point(biWidth * peaks[0], histImage.rows / 2), FONT_HERSHEY_PLAIN, 2, color);

    return stdev;
  }
  float  getLastMode(bool removeTail = true, float rareProp = 0.2, float* smoothParam = NULL, int* peakSize = NULL)
  {
    float mode = 0;
    if (removeTail) {
      float maxHistCount = *max_element(counts.begin(), counts.end());
      float thresh = rareProp*maxHistCount;
      vector<float> histMidsNew, histCountsNew, countsOrigNew;
      for (size_t i = 0; i < mids.size(); i++)
      {
        if (counts[i] > thresh)
        {
          countsOrigNew.push_back(countsOrig[i]);
          histCountsNew.push_back(counts[i]);
          histMidsNew.push_back(mids[i]);
        }
      }
      mids = histMidsNew;
      counts = histCountsNew;
      countsOrig = countsOrigNew;
    }
    vector<float> peakMids;
    vector<int> peaks;

    float smtParam;
    for (smtParam = m_smoothDelta; smtParam < m_smoothMax; smtParam = smtParam + m_smoothDelta)
    {
      smooth(smtParam);
      peaks = getPeaks(peakMids, 0.001);
      if (peaks.size() < (uint)m_tissuesMax) break;
    }
    int idx = 0;
    if (peaks.size() > 0)
    {
      idx = peaks.back();
      mode = mids[idx];
    }
    if (smoothParam != NULL) *smoothParam = smtParam;
    if (peakSize != NULL) *peakSize = peaks.size();
    return mode;
  }
  float getLargestMode(float* smoothParam = NULL, int* peakSize = NULL)
  {
    float mode = 0;
    vector<float> peakMids;
    vector<int> peaks;
    float smtParam;
    for (smtParam = m_smoothDelta; smtParam < m_smoothMax; smtParam = smtParam + m_smoothDelta)
    {
      smooth(smtParam);
      peaks = getPeaks(peakMids, 0.001);

      if (peaks.size() < (uint)m_tissuesMax)
      {
        if (peaks.size() == 0)//TBD need proper fix
        {
          peaks = getPeaks(peakMids, 0.0001);
        }
        break;
      }
    }
    if (peaks.size() > 0)
    {


      float largestNode = counts[peaks[0]];
      mode = mids[peaks[0]];
      for (size_t i = 0; i < peaks.size(); i++)
      {
        if (largestNode < counts[peaks[i]])
        {
          largestNode = counts[peaks[i]];
          mode = mids[peaks[i]];
        }
      }
    }
    else
    {
      mode = 0;
    }
    if (smoothParam != NULL) *smoothParam = smtParam;
    if (peakSize != NULL) *peakSize = peaks.size();
    return mode;
  }
  void getHisInfo(vector<float>& _mids, vector<float>& _origHist, vector<float>& _smoothHist, vector<int>& peakIds)
  {
    if (mids.empty()) return;

    _mids = mids;
    _origHist = countsOrig;
    _smoothHist = counts;
    vector<float> peakMids;
    peakIds = getPeaks(peakMids, 0.001);
  }

};

/**
\class WhiteStripeGUI
\brief Implementation of  WhiteStripe algorithm
Reference:
@article{Shinohara20149,
title = "Statistical normalization techniques for magnetic resonance imaging ",
journal = "NeuroImage: Clinical ",
volume = "6",
pages = "9 - 19",
year = "2014",
issn = "2213-1582",
}
*/
class WhiteStripe {

public:
  WhiteStripe();
  //Set WhiteStripe parameters 
  void setParams(float wsWidth, int sliceStartZ, int sliceStopZ, int tissuesMax, float smoothMax, float smoothDelta, int histSize, bool bT1);
  //Process an mage with optional return of computed mask(maskOut)
  ImageType::Pointer process(ImageType::Pointer img, ImageType::Pointer& maskOut);

  //Return Hist info for display purpose
  void getHisInfo(vector<float>& mids, vector<float>& origHist, vector<float>& smoothHist, vector<int>& peakIds, int& modeId);


private:
  //Private member functions
  ImageType::Pointer whitestripeNorm(ImageType::Pointer img, vector<ImageType::IndexType> indices);
  float threshMean(const vector<float>& voi, float thresh);
  vector<float> makeImageVoi(const ImageType::Pointer img);
  template <typename T>
  double getMeanValue(vector<T> data)
  {

    double mean = 0.0;
    for (size_t i = 0; i < data.size(); i++)
    {
      mean += data[i];
    }
    mean = mean / data.size();
    return mean;
  }
  double calcStdDev(const ImageType::Pointer& img, const vector<ImageType::IndexType>& indices);
  void zeroOutImg(ImageType::Pointer inputImage);
  ImageType::Pointer whitestripeIndToMask(ImageType::Pointer img, vector<ImageType::IndexType> indices);
  vector<ImageType::IndexType> whitestripeMasktoInd(ImageType::Pointer mask);

  //Private Member variables
  HistSruct m_hstObj;
  float m_wsWidth;
  int m_sliceStartZ;
  int m_sliceStopZ;
  int m_tissuesMax;
  float m_smoothMax;
  float m_smoothDelta;
  int m_histSize;
  bool m_bT1;

};

#endif
