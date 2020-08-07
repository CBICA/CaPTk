#pragma once

#include "itkImage.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"

#include <map>
#include <string>
#include <numeric>
#include <vector>

#include "FeatureBase.h"
#include "TextureFeatureBase.h"

#include "cbicaStatistics.h"

template< typename TImageType >
class GLRLMFeatures : 
  public FeatureBase< TImageType >, public TextureFeatureBase < TImageType >
{
public:

  //! Default constructor
  GLRLMFeatures() { };

  //! Default destructor
  ~GLRLMFeatures() { };

  //! Offset selector type
  void SetOffsetSelectorType(const std::string inputOffsetSelector)
  {
    m_offsetSelector = inputOffsetSelector;
  }

  void SetDistanceMax(const double maxDistance)
  {
    m_maxDistance = maxDistance;
  }
  
  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;

      //MITK Version
      //using RunLengthFilterType = itk::Statistics::EnhancedScalarImageToRunLengthFeaturesFilter< TImageType, HistogramFrequencyContainerType >;

      //ITK Version
      using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter< TImageType, HistogramFrequencyContainerType >;

      using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
      using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;

      ////TBD
      //typename RunLengthFilterType::Pointer wrapper_generator = RunLengthFilterType::New();
      //wrapper_generator->SetInput(this->m_inputImage);
      //wrapper_generator->SetMaskImage(this->m_Mask);
      //wrapper_generator->SetInsidePixelValue(1);
      //wrapper_generator->SetPixelValueMinMax(this->m_minimum, this->m_maximum);
      ////TBD
      
      typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
      matrix_generator->SetInput(this->m_inputImage);
      matrix_generator->SetMaskImage(this->m_Mask);
      matrix_generator->SetInsidePixelValue(1);
      
      matrix_generator->SetNumberOfBinsPerAxis(this->m_Bins);
      matrix_generator->SetPixelValueMinMax(this->m_minimum, this->m_maximum);
      if (this->m_histogramBinningType != HistogramBinningType::Equal)
      {
        matrix_generator->SetPixelValueMinMax(this->m_minimum, this->m_maximum);
      }

      if (m_maxDistance != -1) // this will happen for lattice
      {
        matrix_generator->SetDistanceValueMinMax(0, m_maxDistance);
      }

      //wrapper_generator->SetNumberOfBinsPerAxis(this->m_Bins);
      //std::cout << "\n[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - this->m_Bins = " << this->m_Bins << std::endl;

      //TBD
      typename RunLengthFilterType::FeatureNameVectorPointer requestedFeatures = RunLengthFilterType::FeatureNameVector::New();
      typedef typename RunLengthFilterType::RunLengthFeaturesFilterType TextureFilterType;
      requestedFeatures->push_back(TextureFilterType::ShortRunEmphasis);
      requestedFeatures->push_back(TextureFilterType::LongRunEmphasis);
      requestedFeatures->push_back(TextureFilterType::GreyLevelNonuniformity);
      //requestedFeatures->push_back(TextureFilterType::GreyLevelNonuniformityNormalized);
      requestedFeatures->push_back(TextureFilterType::RunLengthNonuniformity);
      //requestedFeatures->push_back(TextureFilterType::RunLengthNonuniformityNormalized);
      requestedFeatures->push_back(TextureFilterType::LowGreyLevelRunEmphasis);
      requestedFeatures->push_back(TextureFilterType::HighGreyLevelRunEmphasis);
      requestedFeatures->push_back(TextureFilterType::ShortRunLowGreyLevelEmphasis);
      requestedFeatures->push_back(TextureFilterType::ShortRunHighGreyLevelEmphasis);
      requestedFeatures->push_back(TextureFilterType::LongRunLowGreyLevelEmphasis);
      requestedFeatures->push_back(TextureFilterType::LongRunHighGreyLevelEmphasis);
      //requestedFeatures->push_back(TextureFilterType::RunPercentage);
      //requestedFeatures->push_back(TextureFilterType::NumberOfRuns);
      //requestedFeatures->push_back(TextureFilterType::GreyLevelVariance);
      //requestedFeatures->push_back(TextureFilterType::RunLengthVariance);
      //requestedFeatures->push_back(TextureFilterType::RunEntropy);
      //wrapper_generator->SetRequestedFeatures(requestedFeatures);
      //TBD

      typename  RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
      typename  RunLengthFeatures::Pointer runLengthFeaturesCalculator = RunLengthFeatures::New();

      size_t offsetNum = 0;


      //std::cout << "\n[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - Set offsetNum = " << offsetNum << std::endl;

      auto size = this->m_inputImage->GetBufferedRegion().GetSize();
      double size_total = size[0];
      for (size_t d = 1; d < TImageType::ImageDimension; d++)
      {
        size_total *= size[d];
      }
      //TBD
      //std::cout << "\n[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - Set size_total = " << size_total << std::endl;
      //matrix_generator->SetDistanceValueMinMax(0, size_total); // TBD: TOCHECK - how is this affecting the computation?
      //TBD

      //TBD - Testing
      //std::cout << "\n[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - SetDistanceValueMinMax - [0, m_Range] = [" << 0 << ", " << size_total << "]" << std::endl;
      //matrix_generator->SetDistanceValueMinMax(0, m_Range); // TBD: TOCHECK - how is this affecting the computation?
      //TBD - Testing



      if ((m_offsetSelector == "Average") || (m_offsetSelector == "Individual"))
      {
        double sre = 0, lre = 0, gln = 0, glnn = 0, rln = 0, rlnn = 0, rp = 0, lglre = 0, hglre = 0, srlgle = 0, srhgle = 0, lrlgle = 0, lrhgle = 0,
          runs = 0, glv = 0, rlv = 0, re = 0;

        //TBD
        int count_offset = 0;
        //std::cout << "\n[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - Average OR Individual  - Set count_offset = " << count_offset << std::endl;
        //TBD

        for (auto offsetIt = this->m_offsets->Begin(); offsetIt != this->m_offsets->End(); offsetIt++, offsetNum++)
        {
          //TBD
          //std::cout << "\n\n";
          //std::cout << "[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - Average OR Individual  - GLRLM Matrix: Offset: " << offsetIt.Value() << "\n";
          //TBD

          if (m_maxDistance == -1) // this will happen for non-lattice
          {
            matrix_generator->SetDistanceValueMinMax(0, GetDistanceFromOffset(offsetIt.Value()));
          }

          matrix_generator->SetOffset(offsetIt.Value());
          matrix_generator->Update();



          //TBD
          //std::cout << "\n\n";
          //std::cout << "[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - Average OR Individual  - GLRLM Matrix: Offset: " << offsetIt.Value() << "(" << count_offset << ")\n";
          //auto temp = matrix_generator->GetOutput();
          //for (auto iter = temp->Begin(); iter != temp->End(); ++iter)
          //{
          //  std::cout << "\tIndex = " << iter.GetIndex() << " ||| Frequency = " << iter.GetFrequency() << std::endl;
          //}

          //std::cout << "[DEBUG] GLRLM Matrix: Offset: " << offsetIt.Value() << "\n";
          //std::cout << "\tindex\t|\t|\tfrenquency" << std::endl;
          //for (int bin_count = 0; bin_count < this->m_Bins; bin_count++)
          //{
          //  std::cout << "\t" << bin_count << "\t|\t" << temp->GetFrequency(bin_count) << std::endl;
          //}
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - GLRLM Matrix: Offset: " << offsetIt.Value() << " | Distance[" << matrix_generator->GetMinDistance() << ", " << matrix_generator->GetMaxDistance() << "] | Pixel [" << matrix_generator->GetMin() << ", " << matrix_generator->GetMax() << "] \n" << std::endl;
          //TBD


          ////Header for easier DEBUG std:cout interpretation
          //std::cout << "\tGLRLM Measurement Vectors -> intensity[]\t|\trunLength[] = midpoint[int,run]\t|\tfrenquency" << std::endl;

          //TBD
          auto temp = matrix_generator->GetOutput();
          for (auto iter = temp->Begin(); iter != temp->End(); ++iter)
          {
            auto temp_index = iter.GetIndex();
            auto temp_min_vector = temp->GetHistogramMinFromIndex(temp_index);
            auto temp_max_vector = temp->GetHistogramMaxFromIndex(temp_index);
            //std::cout << " \tGLRLM Measurement Vectors -> Intensity[" << temp_min_vector[0] << ", " << temp_max_vector[0] << "]\t|\tRun[" << temp_min_vector[1] << ", " << temp_max_vector[1] << "] = [" << temp->GetMeasurementVector(temp_index) << "]\t|\tFrequency = " << temp->GetFrequency(temp_index) << std::endl;
            //if (temp->GetFrequency(temp_index) > 0) {
            //  std::cout << " \tGLRLM Measurement Vectors -> Intensity[" << temp_min_vector[0] << ", " << temp_max_vector[0] << "]\t|\tRun[" << temp_min_vector[1] << ", " << temp_max_vector[1] << "] | Index[" << iter.GetIndex() << "] = " << temp->GetFrequency(temp_index) << std::endl;
            //}
          }
          //TBD

          runLengthFeaturesCalculator->SetInput(matrix_generator->GetOutput());
          runLengthFeaturesCalculator->Update();

          //TBD - print out matrix
          //std::cout << "\n\n\n";
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - matrix_generator->GetNumberOfBinsPerAxis() = " << matrix_generator->GetNumberOfBinsPerAxis() << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - matrix_generator->GetMinDistance() = " << matrix_generator->GetMinDistance() << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - matrix_generator->GetMaxDistance() = " << matrix_generator->GetMaxDistance() << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - matrix_generator->GetOffsets() = " << matrix_generator->GetOffsets() << std::endl;

          //TBD - Only prints out first set, may not be complete - std cout codes above for full matrix
          //std::cout << "\n[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Individual -> " << offsetIt.Value() << " - Matrix = \n" << std::endl;
          //std::cout << "\tindex\t|\t[min, max]\t|\tfrenquency" << std::endl;
          //for (int bin_count = 0; bin_count < this->m_Bins; bin_count++) 
          //{
          //  auto index_temp = matrix_generator->GetOutput()->GetIndex(bin_count);
          //  auto min_temp = matrix_generator->GetOutput()->GetHistogramMinFromIndex(index_temp);
          //  auto max_temp = matrix_generator->GetOutput()->GetHistogramMaxFromIndex(index_temp);
          //  auto min_temp2 = roundf(min_temp[0] * 100) / 100;
          //  auto max_temp2 = roundf(max_temp[0] * 100) / 100;
          //  auto freq_temp = matrix_generator->GetOutput()->GetFrequency(index_temp);
          //  std::cout << "\t" << "index[" << index_temp << "] \t | \t" << bin_count << "\t | \t[" << min_temp2 << ", " << max_temp2 << "]\t | \t" << freq_temp << std::endl;
          //}
          //TBD - Only prints out first set, may not be complete

          //std::cout << "\tLowGreyLevelRunEmphasis = " << runLengthFeaturesCalculator->GetLowGreyLevelRunEmphasis() << std::endl;
          //std::cout << "\tHighGreyLevelRunEmphasis = " << runLengthFeaturesCalculator->GetHighGreyLevelRunEmphasis() << std::endl;
          //std::cout << "\tShortRunLowGreyLevelEmphasis = " << runLengthFeaturesCalculator->GetShortRunLowGreyLevelEmphasis() << std::endl;
          //std::cout << "\tShortRunHighGreyLevelEmphasis = " << runLengthFeaturesCalculator->GetShortRunHighGreyLevelEmphasis() << std::endl;
          //std::cout << "\tLongRunLowGreyLevelEmphasis = " << runLengthFeaturesCalculator->GetLongRunLowGreyLevelEmphasis() << std::endl;
          //std::cout << "\tLongRunHighGreyLevelEmphasis = " << runLengthFeaturesCalculator->GetLongRunHighGreyLevelEmphasis() << std::endl;
          //TBD - print out matrix

          //TBD
          count_offset++;
          //TBD

          if (m_offsetSelector == "Average")
          {
            sre += runLengthFeaturesCalculator->GetShortRunEmphasis();
            lre += runLengthFeaturesCalculator->GetLongRunEmphasis();
            gln += runLengthFeaturesCalculator->GetGreyLevelNonuniformity();
            rln += runLengthFeaturesCalculator->GetRunLengthNonuniformity();
            lglre += runLengthFeaturesCalculator->GetLowGreyLevelRunEmphasis();
            hglre += runLengthFeaturesCalculator->GetHighGreyLevelRunEmphasis();
            srlgle += runLengthFeaturesCalculator->GetShortRunLowGreyLevelEmphasis();
            srhgle += runLengthFeaturesCalculator->GetShortRunHighGreyLevelEmphasis();
            lrlgle += runLengthFeaturesCalculator->GetLongRunLowGreyLevelEmphasis();
            lrhgle += runLengthFeaturesCalculator->GetLongRunHighGreyLevelEmphasis();
            runs += runLengthFeaturesCalculator->GetTotalNumberOfRuns();
            rp += static_cast<double>(runLengthFeaturesCalculator->GetTotalNumberOfRuns()) / static_cast<double>(this->m_nonZeroIndeces.size());
            //// part of the "enhanced" set coming from MITK
            //glnn += runLengthFeaturesCalculator->GetGreyLevelNonuniformityNormalized();
            //rlnn += runLengthFeaturesCalculator->GetRunLengthNonuniformityNormalized();
            //glv += runLengthFeaturesCalculator->GetGreyLevelVariance();
            //rlv += runLengthFeaturesCalculator->GetRunLengthVariance();
            //re += runLengthFeaturesCalculator->GetRunEntropy();
          }
          else // individual
          {
            this->m_features["LowGreyLevelRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetLowGreyLevelRunEmphasis();
            this->m_features["HighGreyLevelRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetHighGreyLevelRunEmphasis();
            this->m_features["ShortRunLowGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetShortRunLowGreyLevelEmphasis();
            this->m_features["ShortRunHighGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetShortRunHighGreyLevelEmphasis();

            // TBD: this features are not IBSI-compliant, hence in debugMode; re-enabling to extract NatSci features
            //if (this->m_debugMode)
            {
              this->m_features["ShortRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetShortRunEmphasis();
              this->m_features["LongRunEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetLongRunEmphasis();
              this->m_features["GreyLevelNonuniformity_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetGreyLevelNonuniformity();
              this->m_features["RunLengthNonuniformity_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetRunLengthNonuniformity();
              this->m_features["LongRunLowGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetLongRunLowGreyLevelEmphasis();
              this->m_features["LongRunHighGreyLevelEmphasis_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetLongRunHighGreyLevelEmphasis();
              this->m_features["TotalRuns_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetTotalNumberOfRuns();
              this->m_features["RunPercentage_Offset_" + std::to_string(offsetNum)] = this->m_features["TotalRuns_Offset_" + std::to_string(offsetNum)] / static_cast<double>(this->m_nonZeroIndeces.size());
              //// part of the "enhanced" set coming from MITK
              //this->m_features["GreyLevelNonuniformityNormalized_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetGreyLevelNonuniformityNormalized();
              //this->m_features["RunLengthNonuniformityNormalized_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetRunLengthNonuniformityNormalized();
              //this->m_features["GreyLevelVariance_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetGreyLevelVariance();
              //this->m_features["RunLengthVariance_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetRunLengthVariance();
              //this->m_features["RunEntropy_Offset_" + std::to_string(offsetNum)] = runLengthFeaturesCalculator->GetRunEntropy();
            }
          }
        }

        if (m_offsetSelector == "Average")
        {
          //TBD
          //std::cout << "\n[DEBUG] - FeatureExtraction.hxx - CalculateGLRLM - this->m_offsets->size() = " << this->m_offsets->size() << std::endl;
          //TBD

          sre /= this->m_offsets->size();
          lre /= this->m_offsets->size();
          gln /= this->m_offsets->size();
          rln /= this->m_offsets->size();
          lglre /= this->m_offsets->size();
          hglre /= this->m_offsets->size();
          srlgle /= this->m_offsets->size();
          srhgle /= this->m_offsets->size();
          lrlgle /= this->m_offsets->size();
          lrhgle /= this->m_offsets->size();
          rp /= this->m_offsets->size();
          runs /= this->m_offsets->size();
          rlnn /= this->m_offsets->size();
          glnn /= this->m_offsets->size();
          glv /= this->m_offsets->size();
          rlv /= this->m_offsets->size();
          re /= this->m_offsets->size();

          this->m_features["LowGreyLevelRunEmphasis"] = lglre;
          this->m_features["HighGreyLevelRunEmphasis"] = hglre;
          this->m_features["ShortRunLowGreyLevelEmphasis"] = srlgle;
          this->m_features["ShortRunHighGreyLevelEmphasis"] = srhgle;

          // TBD: this features are not IBSI-compliant, hence in debugMode; re-enabling to extract NatSci features
          //if (this->m_debugMode)
          {
             this->m_features["ShortRunEmphasis"] = sre;
             this->m_features["LongRunEmphasis"] = lre;
             this->m_features["GreyLevelNonuniformity"] = gln;
             this->m_features["RunLengthNonuniformity"] = rln;
             this->m_features["RunPercentage"] = rp;
             this->m_features["LongRunLowGreyLevelEmphasis"] = lrlgle;
             this->m_features["LongRunHighGreyLevelEmphasis"] = lrhgle;
             this->m_features["TotalRuns"] = runs;
             /// advanced set from MITK
             //this->m_features["RunLengthNonuniformityNormalized"] = rlnn;
             //this->m_features["GreyLevelNonuniformityNormalized"] = glnn;
             //this->m_features["GreyLevelVariance_Offset"] = glv;
             //this->m_features["RunLengthVariance_Offset"] = rlv;
             //this->m_features["RunEntropy"] = re;
          }

          //TBD
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Average - LowGreyLevelRunEmphasis = " << lglre << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Average - HighGreyLevelRunEmphasis = " << hglre << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Average - ShortRunLowGreyLevelEmphasis = " << srlgle << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Average - ShortRunHighGreyLevelEmphasis = " << srhgle << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Average - LongRunLowGreyLevelEmphasis = " << lrlgle << std::endl;
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - Average - LongRunHighGreyLevelEmphasis = " << lrhgle << std::endl;
          //TBD

        }

        ////TBD - test if wrapper (image -> features) give same output as image -> matrix and matrix -> features with intermediate steps
        //if (m_offsetSelect == "Average") 
        //{
        //  std::cout << "\n[DEBUG] FeatureExtraction.hxx - CalculateGLRLM - count_offset = " << count_offset << std::endl;
        //
        //  wrapper_generator->SetOffsets(offset);
        //  wrapper_generator->SetDistanceValueMinMax(0, m_Range);
        //  wrapper_generator->Update();
        //  auto featureMeans = wrapper_generator->GetFeatureMeans();
        //  auto featureStd = wrapper_generator->GetFeatureStandardDeviations();
        //
        //  //typedef typename RunLengthFilterType::RunLengthFeaturesFilterType TextureFilterType;
        //  for (std::size_t i = 0; i < featureMeans->size(); ++i)
        //  {
        //    switch (i)
        //    {
        //    case TextureFilterType::ShortRunEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - ShortRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //
        //    case TextureFilterType::GreyLevelNonuniformity:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - GreyLevelNonuniformity: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::RunLengthNonuniformity:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - RunLengthNonuniformity: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::LowGreyLevelRunEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - LowGreyLevelRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::HighGreyLevelRunEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - HighGreyLevelRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::ShortRunLowGreyLevelEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - ShortRunLowGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::ShortRunHighGreyLevelEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - ShortRunHighGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::LongRunLowGreyLevelEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - LongRunLowGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    case TextureFilterType::LongRunHighGreyLevelEmphasis:
        //      std::cout << "\n [DEBUG] featureextraction.hxx - calculateglrlm - LongRunHighGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //      break;
        //    }
        //  }
        //}
      }
      else if ((m_offsetSelector == "ITKDefault") || (m_offsetSelector == "Combined"))
      {
        if (m_maxDistance == -1) // this will happen for non-lattice
        {
          // find the maximum distance for all offsets
          for (auto offsetIt = this->m_offsets->Begin(); offsetIt != this->m_offsets->End(); offsetIt++, offsetNum++)
          {
            auto currentDistance = GetDistanceFromOffset(offsetIt.Value());
            if (m_maxDistance < currentDistance)
            {
              m_maxDistance = currentDistance;
            }
          }
          matrix_generator->SetDistanceValueMinMax(0, m_maxDistance);
        }

        matrix_generator->SetOffsets(this->m_offsets);
        matrix_generator->Update();

        //TBD - show bin to frequency
        //auto temp = matrix_generator->GetOutput();
        //std::cout << "GLRLM Matrix:\n";
        //for (int bin_count = 0; bin_count < this->m_Bins; bin_count++) 
        //{
        //  std::cout << "\t" << bin_count << "\t|\t" << temp->GetFrequency(bin_count) << std::endl;
        //}
        //TBD - show bin to frequency

        runLengthFeaturesCalculator->SetInput(matrix_generator->GetOutput());
        runLengthFeaturesCalculator->Update();
        
        this->m_features["LowGreyLevelRunEmphasis"] = runLengthFeaturesCalculator->GetLowGreyLevelRunEmphasis();
        this->m_features["HighGreyLevelRunEmphasis"] = runLengthFeaturesCalculator->GetHighGreyLevelRunEmphasis();
        this->m_features["ShortRunLowGreyLevelEmphasis"] = runLengthFeaturesCalculator->GetShortRunLowGreyLevelEmphasis();
        this->m_features["ShortRunHighGreyLevelEmphasis"] = runLengthFeaturesCalculator->GetShortRunHighGreyLevelEmphasis();


        // TBD: this features are not IBSI-compliant, hence in debugMode; re-enabling to extract NatSci features
        //if (this->m_debugMode)
        {
           this->m_features["ShortRunEmphasis"] = runLengthFeaturesCalculator->GetShortRunEmphasis();
           this->m_features["LongRunEmphasis"] = runLengthFeaturesCalculator->GetLongRunEmphasis();
           this->m_features["GreyLevelNonuniformity"] = runLengthFeaturesCalculator->GetGreyLevelNonuniformity();
           this->m_features["RunLengthNonuniformity"] = runLengthFeaturesCalculator->GetRunLengthNonuniformity();
           this->m_features["LongRunLowGreyLevelEmphasis"] = runLengthFeaturesCalculator->GetLongRunLowGreyLevelEmphasis();
           this->m_features["LongRunHighGreyLevelEmphasis"] = runLengthFeaturesCalculator->GetLongRunHighGreyLevelEmphasis();
        }

        //TBD
        //wrapper_generator->SetOffsets(offset);
        //wrapper_generator->CombinedFeatureCalculationOn(); //function not available according to VS2017. Error Code C2039

        //wrapper_generator->Update();
        //auto featureMeans = wrapper_generator->GetFeatureMeans();
        //auto featureStd = wrapper_generator->GetFeatureStandardDeviations();

        //for (std::size_t i = 0; i < featureMeans->size(); ++i)
        //{
        //  switch (i)
        //  {
        //  case TextureFilterType::ShortRunEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - ShortRunEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetShortRunEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - ShortRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::LongRunEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LongRunEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetLongRunEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LongRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::GreyLevelNonuniformity:
        //    break;
        //  //case texturefiltertype::greylevelnonuniformitynormalized:
        //  //  break;
        //  case TextureFilterType::RunLengthNonuniformity:
        //    break;
        //  //case TextureFilterType::RunLengthNonuniformityNormalized:
        //  //  break;
        //  case TextureFilterType::LowGreyLevelRunEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LowGreyLevelRunEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetLowGreyLevelRunEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LowGreyLevelRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::HighGreyLevelRunEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - HighGreyLevelRunEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetHighGreyLevelRunEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - HighGreyLevelRunEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::ShortRunLowGreyLevelEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - ShortRunLowGreyLevelEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetShortRunLowGreyLevelEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - ShortRunLowGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::ShortRunHighGreyLevelEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - ShortRunHighGreyLevelEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetShortRunHighGreyLevelEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - ShortRunHighGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::LongRunLowGreyLevelEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LongRunLowGreyLevelEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetLongRunLowGreyLevelEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LongRunLowGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::LongRunHighGreyLevelEmphasis:
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LongRunHighGreyLevelEmphasis: runLengthMatrixCalculator = " << runLengthMatrixCalculator->GetLongRunHighGreyLevelEmphasis() << std::endl;
        //    std::cout << "\n [DEBUG] FeatureExtraction.hxx - CalculateGLRLM - LongRunHighGreyLevelEmphasis: wrapper_generator = " << featureMeans->ElementAt(i) << std::endl;
        //    break;
        //  case TextureFilterType::RunPercentage:
        //    break;
        //  case TextureFilterType::NumberOfRuns:
        //    break;
        //  case TextureFilterType::GreyLevelVariance:
        //    break;
        //  case TextureFilterType::RunLengthVariance:
        //    break;
        //  case TextureFilterType::RunEntropy:
        //    break;
        //  default:
        //    break;
        //  }
        //}
        //TBD
      }
      else
      {
        // not defined, so don't do anything to this->m_features
      }

      this->m_algorithmDone = true;
    }
  }

private:

  //! Compute the distance of the current offset from the center in world coordinate system
  double GetDistanceFromOffset(typename TImageType::OffsetType currentOffset)
  {
    auto spacing = this->m_inputImage->GetSpacing();
    double temp, distance = 0;
    for (size_t d = 0; d < TImageType::ImageDimension; d++)
    {
      temp = std::abs(currentOffset[d]) * spacing[d]; // abs is needed because offset can be negative compared to the center
      distance += temp * temp;
    }
    return std::sqrt(distance);
  }

  std::string m_offsetSelector; //! type of offset selection
  double m_maxDistance = -1;
};
