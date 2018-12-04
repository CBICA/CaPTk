#ifndef _UTILS_H
#define _UTILS_H

#include <fstream>
#include <iostream>
#include <iterator>
#include <ctime>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include <itkAbsoluteValueDifferenceImageFilter.h>
#include "itkStatisticsImageFilter.h"

#include <itkNiftiImageIO.h>
#include <itkImageFileWriter.h>
#include "itkImageFileReader.h"
#include <inttypes.h>
#include <sstream>

using namespace std;
typedef int(*TXT_CALLBACK_TYPE)(const char*);

typedef float PixelType;
typedef itk::Image<PixelType, 3> ImageType;

namespace Utils
{
  //Wh reusable code
  inline vector<ImageType::IndexType> getPixelIndicesThresholded(const ImageType::Pointer inputImage, const float lowerThreshold)
  {
    typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
    ConstIteratorType in(inputImage, inputImage->GetRequestedRegion());
    vector<ImageType::IndexType> res;
    for (in.GoToBegin(); !in.IsAtEnd(); ++in)
    {
      if (in.Get() > lowerThreshold)
      {
        res.push_back(in.GetIndex());
      }
    }
    return res;
  }
  inline vector<ImageType::IndexType> getPixelIndicesThresholded(const ImageType::Pointer inputImage, const float minVal, const float maxVal)
  {
    typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
    ConstIteratorType in(inputImage, inputImage->GetRequestedRegion());
    vector<ImageType::IndexType> res;
    for (in.GoToBegin(); !in.IsAtEnd(); ++in)
    {
      if (in.Get()>minVal && in.Get()<maxVal)
      {
        res.push_back(in.GetIndex());
      }
    }
    return res;
  }
  inline  vector<float> quantileFast(vector<float> data, vector<float> probs)
  {
    vector<float> qVals;
    for (size_t i = 0; i < probs.size(); i++)
    {
      int idx = int(data.size()*probs[i] + 0.5) - 1;
      nth_element(data.begin(), data.begin() + idx, data.end());
      qVals.push_back(data[idx]);
    }
    return qVals;
  }
  inline vector<float> quantile(vector<float> data, vector<float> probs)
  {
    vector<float> qVals;
    sort(data.begin(), data.end());
    for (size_t i = 0; i < probs.size(); i++)
    {
      uint idx = uint(data.size()*probs[i] + 0.5) - 1;
      if (idx < 0)  idx = 0;
      if (idx >= data.size())  idx = data.size() - 1;
      qVals.push_back(data[idx]);
    }
    return qVals;
  }
  inline double getMeanValue(const ImageType::Pointer inputImage)
  {
    typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
    ConstIteratorType in(inputImage, inputImage->GetRequestedRegion());
    double mean = 0.0;
    long long length = 0;
    for (in.GoToBegin(); !in.IsAtEnd(); ++in)
    {
      mean += in.Get();
      length++;
    }
    mean = mean / length;
    return mean;
  }


	//Utility 
  inline string getCmdOption(char ** begin, char ** end, const std::string & option)
	{
		char ** itr = std::find(begin, end, option);
		if (itr != end && ++itr != end)
		{
			return string(*itr);
		}
		return string();
	}
  inline  bool cmdOptionExists(char** begin, char** end, const std::string& option)
	{
		return std::find(begin, end, option) != end;
	}

	//Strings functions
  inline void replaceSubStr(string &s, const string &search, const string &replace)
	{
		for (size_t pos = 0;; pos += replace.length()) 
		{
			pos = s.find(search, pos);
			if (pos == string::npos) break;
			s.erase(pos, search.length());
			s.insert(pos, replace);
		}
	}

	//File finctions
	inline bool fileExist(const std::string& name) {
		ifstream f(name.c_str());
		return f.good();
	}
  inline string toString(const int& value)
	{
		std::ostringstream ostr;
		ostr << value;
		return  ostr.str(); 
	}
  inline bool pixelMatchImages(itk::SmartPointer<ImageType> baseImg, itk::SmartPointer<ImageType> testImg, double intensityTolerance = 2.0, int radiusTolerance = 0)
  {
    if (baseImg.IsNull() || testImg.IsNull())
    {
      return false;
    }
    if (intensityTolerance < 0.0)
    {
      intensityTolerance = 0.0;
    }
    if (radiusTolerance < 0)
    {
      radiusTolerance = 0;
    }

    ImageType::SizeType baselineSize = baseImg->GetLargestPossibleRegion().GetSize();
    ImageType::SizeType testSize = testImg->GetLargestPossibleRegion().GetSize();
    if (baselineSize != testSize)
    {
      return false;
    }

    // Now compare the two images
    typedef itk::AbsoluteValueDifferenceImageFilter<ImageType, ImageType, ImageType> DiffType;
    DiffType::Pointer diff = DiffType::New();
    diff->SetInput1(baseImg);
    diff->SetInput1(testImg);
    diff->Update();


    typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
    StatisticsImageFilterType::Pointer statisticsImageFilter
      = StatisticsImageFilterType::New();
    statisticsImageFilter->SetInput(diff->GetOutput());
    statisticsImageFilter->Update();

    float meanError = statisticsImageFilter->GetMean();
    if (meanError < 0.001)
    {
      return true;
    }
    return false;
  }
  inline itk::SmartPointer<ImageType> readNifti(const string& fileName)
	{
		typedef itk::ImageFileReader<ImageType> ReaderType;
    itk::SmartPointer<ImageType>  image = ImageType::New();
		try
		{
			ReaderType::Pointer reader = ReaderType::New();
			reader->SetFileName(fileName);
			reader->Update();
			image->Graft(reader->GetOutput());
		}
		catch (...)
		{
			cout<<"Exception in reading Nifti:" + fileName;
      return NULL;
		}
		return image;
	}
  inline void writeNifti(const string& inputFile, const string& outFile, const ImageType::Pointer image)
	{
		itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(inputFile.c_str(), itk::ImageIOFactory::ReadMode);
		imageIO->SetFileName(inputFile);
		imageIO->ReadImageInformation();
		itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();
		nifti_io->SetPixelType(imageIO->GetPixelType());

		itk::ImageFileWriter<ImageType>::Pointer writer = itk::ImageFileWriter<ImageType>::New();
		writer->SetFileName(outFile);
		writer->SetInput(image);
		writer->SetImageIO(nifti_io);
		writer->Update();
	}
}
#endif
