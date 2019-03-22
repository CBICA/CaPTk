#pragma once

#include "itkImage.h"
#include "itkHistogram.h"
#include "itkNumericTraits.h"
#include "itkVectorContainer.h"
#include "itkNeighborhoodIterator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMaskImageFilter.h"

#include <map>
#include <string>
#include <numeric>
#include <vector>

#include "FeatureBase.h"

#include "cbicaStatistics.h"

template< typename TImageType >
class HistogramFeatures : public FeatureBase < TImageType >
{
public:
  //! Default constructor
  HistogramFeatures() { };

  //! Default destructor
  ~HistogramFeatures() { };
  
  void SetMinimum(typename TImageType::PixelType min) { m_minimum = min; };

  void SetMaximum(typename TImageType::PixelType max) { m_maximum = max; };
  /**
  \brief Set the number of bins for quantization of the input image; defaults to 10 unless overridden by user

  \param numBinValue integer value for the number of bins you want for image quantization
  **/
  void SetNumBins(unsigned int numBinValue)
  {
    m_Bins = numBinValue;
  }

  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      using HistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter< TImageType, TImageType >;
      using HistogramMeasurementType = typename HistogramFilterType::HistogramType::MeasurementVectorType;

      HistogramMeasurementType lowerBound(1), upperBound(1);
      lowerBound.Fill(m_minimum);
      upperBound.Fill(m_maximum);

      typename HistogramFilterType::HistogramType::SizeType size(1); // this is always grayscale
      size.Fill(m_Bins);

      auto histogramCalculator = typename HistogramFilterType::New();
      histogramCalculator->SetInput(this->m_inputImage);
      histogramCalculator->SetMaskImage(this->m_Mask);
      histogramCalculator->SetMaskValue(1);
      histogramCalculator->SetHistogramBinMinimum(lowerBound);
      histogramCalculator->SetHistogramBinMaximum(upperBound);
      histogramCalculator->SetHistogramSize(size);

      try
      {
        histogramCalculator->Update();
      }
      catch (itk::ExceptionObject & error)
      {
        std::cerr << "Error: " << error << "\n";
      }

      auto histogram = histogramCalculator->GetOutput();

      auto totalFrequency = histogram->GetTotalFrequency();
      std::vector< typename TImageType::PixelType > histogramVector/*(histogram->GetTotalFrequency())*/;
      double entropy = 0, uniformity = 0;

      for (size_t bin = 0; bin < histogram->Size(); ++bin)
      {
        auto currentFrequency = histogram->GetFrequency(bin, 0);
        for (size_t f = 0; f < currentFrequency; ++f)
        {
          histogramVector.push_back(bin);
        }
        auto prob = static_cast<double>(currentFrequency) / static_cast<double>(totalFrequency);
        this->m_features["Bin-" + std::to_string(bin) + "_Probability"] = prob;
        this->m_features["Bin-" + std::to_string(bin) + "_Frequency"] = currentFrequency;
        this->m_features["Bin-" + std::to_string(bin) + "_Max"] = histogram->GetBinMax(0, bin);
        this->m_features["Bin-" + std::to_string(bin) + "_Min"] = histogram->GetBinMin(0, bin);

        entropy += prob * std::log2(prob);
        uniformity += std::pow(prob, 2); // IBSI 3.4.19
      }

      entropy = -entropy; // IBSI 3.4.18

      cbica::Statistics< typename TImageType::PixelType > histogramStatsCalculator;
      histogramStatsCalculator.SetInput(histogramVector);

      typename FeatureBase< TImageType >::TConstIteratorType imageIterator(this->m_inputImage, this->m_inputImage->GetBufferedRegion()),
        maskIterator(this->m_Mask, this->m_Mask->GetBufferedRegion());
      auto fifthPercentileValue = histogramStatsCalculator.GetNthPercentileElement(5);
      auto ninetyFifthPercentileValue = histogramStatsCalculator.GetNthPercentileElement(95);

      auto fifthPercentileMean = 0.0;
      auto ninetyFifthPercentileMean = 0.0;
      double fifthN = 0.0;
      double ninetyFifthN = 0.0;
      double N = 0.0;

      for (maskIterator.GoToBegin(); !maskIterator.IsAtEnd(); ++maskIterator)
      {
        if (maskIterator.Get() > 0)
        {
          imageIterator.SetIndex(maskIterator.GetIndex());
          auto currentValue = imageIterator.Get();

          N++;

          if (currentValue <= fifthPercentileValue)
          {
            fifthPercentileMean += currentValue;
            fifthN++;
          }
          else if (currentValue >= ninetyFifthPercentileValue)
          {
            ninetyFifthPercentileMean += currentValue;
            ninetyFifthN++;
          }
        }
      }

      // calculate final statistics
      fifthPercentileMean /= fifthN;
      ninetyFifthPercentileMean /= ninetyFifthN;

      this->m_features["Minimum"] = histogramStatsCalculator.GetMinimum();
      this->m_features["Maximum"] = histogramStatsCalculator.GetMaximum();
      this->m_features["Mean"] = histogramStatsCalculator.GetMean();
      this->m_features["Sum"] = histogramStatsCalculator.GetSum();
      this->m_features["Mode"] = histogramStatsCalculator.GetMode();
      this->m_features["Median"] = histogramStatsCalculator.GetMedian();
      this->m_features["Variance"] = histogramStatsCalculator.GetVariance();
      this->m_features["StandardDeviation"] = histogramStatsCalculator.GetStandardDeviation();
      this->m_features["Skewness"] = histogramStatsCalculator.GetSkewness();
      this->m_features["Kurtosis"] = histogramStatsCalculator.GetKurtosis();
      this->m_features["Range"] = histogramStatsCalculator.GetRange();
      this->m_features["Energy"] = histogramStatsCalculator.GetEnergy();
      this->m_features["Entropy"] = entropy;
      this->m_features["Uniformity"] = uniformity;
      this->m_features["RootMeanSquare"] = histogramStatsCalculator.GetRootMeanSquare();
      this->m_features["TenthPercentile"] = histogramStatsCalculator.GetNthPercentileElement(10);
      this->m_features["NinetiethPercentile"] = histogramStatsCalculator.GetNthPercentileElement(90);
      this->m_features["FifthPercentile"] = fifthPercentileValue;
      this->m_features["NinetyFifthPercentile"] = ninetyFifthPercentileValue;
      this->m_features["FifthPercentileMean"] = fifthPercentileMean;
      this->m_features["NinetyFifthPercentileMean"] = ninetyFifthPercentileMean;
      this->m_features["TwentyFifthPercentile"] = histogramStatsCalculator.GetNthPercentileElement(25);
      this->m_features["SeventyFifthPercentile"] = histogramStatsCalculator.GetNthPercentileElement(75);
      this->m_features["InterQuartileRange"] = histogramStatsCalculator.GetInterQuartileRange();
      this->m_features["MeanAbsoluteDeviation"] = histogramStatsCalculator.GetMeanAbsoluteDeviation();
      this->m_features["RobustMeanAbsoluteDeviation1090"] = histogramStatsCalculator.GetRobustMeanAbsoluteDeviation(10, 90);
      this->m_features["MedianAbsoluteDeviation"] = histogramStatsCalculator.GetMedianAbsoluteDeviation();
      this->m_features["CoefficientOfVariation"] = histogramStatsCalculator.GetCoefficientOfVariation();
      this->m_features["QuartileCoefficientOfVariation"] = histogramStatsCalculator.GetQuartileCoefficientOfDispersion();
      /// histogram calculation from ITK -- for texture feature pipeline ends

      for (size_t i = 0; i < histogram->GetSize()[0]; i++)
      {
        this->m_features["Bin_" + std::to_string(i) + "_Frequency"] = histogram->GetFrequency(i);
        this->m_features["Bin_" + std::to_string(i) + "_Max"] = histogram->GetBinMax(0, i);
        this->m_features["Bin_" + std::to_string(i) + "_Min"] = histogram->GetBinMin(0, i);
      }

      const float maxRescaleVal = 1000;
      auto rescaleFilter = itk::RescaleIntensityImageFilter< TImageType, TImageType >::New();
      rescaleFilter->SetInput(this->m_inputImage);
      rescaleFilter->SetOutputMinimum(0);
      rescaleFilter->SetOutputMaximum(maxRescaleVal);
      rescaleFilter->Update();

      std::vector< double > intensities;
      itk::ImageRegionConstIterator <TImageType> 
        imageIt(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion()),
        maskIt(this->m_Mask, this->m_Mask->GetLargestPossibleRegion());
      imageIt.GoToBegin();
      maskIt.GoToBegin();

      while (!imageIt.IsAtEnd())
      {
        if (maskIt.Get() > 0)
          intensities.push_back(std::round(imageIt.Get()));
        ++imageIt;
        ++maskIt;
      }

      //-------------------------------
      double interval = maxRescaleVal / m_Bins;
      double final_interval = (int)(interval * 100);
      final_interval = (double)final_interval / 100;


      std::vector<double> finalBins;
      std::vector<std::vector<double>> Ranges;
      double current_index = 0;
      for (int i = 0; i < m_Bins; i++)
      {
        std::vector<double> onerange;
        onerange.push_back(current_index);
        current_index = current_index + final_interval;
        onerange.push_back(current_index);
        Ranges.push_back(onerange);

        if (static_cast<int>(Ranges.size()) == m_Bins)
          Ranges[Ranges.size() - 1][1] = maxRescaleVal;
      }
      //toadd the last one
      for (unsigned int j = 0; j < Ranges.size(); j++)
      {
        std::vector<double> onerange = Ranges[j];
        int counter = 0;
        for (unsigned int i = 0; i < intensities.size(); i++)
        {
          if (onerange[0] == 0)
          {
            if (intensities[i] >= onerange[0] && intensities[i] <= onerange[1])
              counter = counter + 1;
          }
          else
          {
            if (intensities[i] > onerange[0] && intensities[i] <= onerange[1])
              counter = counter + 1;
          }
        }
        finalBins.push_back(counter);
      }
      for (unsigned int j = 0; j < finalBins.size(); j++)
      {
        finalBins[j] = (finalBins[j] * 100) / intensities.size();
        this->m_features["Bin_" + std::to_string(j)] = finalBins[j];
        this->m_features["BinEndIntensity_" + std::to_string(j)] = Ranges[j][1];
      }
      this->m_algorithmDone = true;
    }
  };

  /**
  \brief return the map of feature names and feature values
  **/
  std::map< std::string, double > GetOutput()
  {
    if (!this->m_algorithmDone)
    {
      Update();
    }
    return this->m_features;
  };

private:

  unsigned int m_Bins;
  typename TImageType::PixelType m_minimum = 0, m_maximum = 0;

  itk::Statistics::Histogram< double >::Pointer m_histogram;
};
