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
class GLCMFeatures : 
  public FeatureBase< TImageType >, public TextureFeatureBase < TImageType >
{
public:

  //! Default constructor
  GLCMFeatures() { };

  //! Default destructor
  ~GLCMFeatures() { };

  //! Offset selector type
  void SetOffsetSelectorType(const std::string inputOffsetSelector)
  {
    m_offsetSelector = inputOffsetSelector;
  }
  
  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter < TImageType >;
      using HistogramType = typename Image2CoOccuranceType::HistogramType;
      using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter< HistogramType >;

      double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

      auto image_wrap = this->m_inputImage;
      auto mask_wrap = this->m_Mask;

      if (m_offsetSelector == "Average")
      {
        for (size_t i = 0; i < this->m_offsets->size(); i++)
        {
          typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
          glcmGenerator->SetNumberOfBinsPerAxis(this->m_Bins); //reasonable number of bins
          if (this->m_histogramBinningType != HistogramBinningType::Equal)
          {
            glcmGenerator->SetPixelValueMinMax(this->m_minimum, this->m_maximum);
          }
          glcmGenerator->SetMaskImage(mask_wrap);
          glcmGenerator->SetInput(image_wrap);
          auto featureCalc = Hist2FeaturesType::New();

          glcmGenerator->SetOffset(this->m_offsets->at(i));
          glcmGenerator->Update();
          featureCalc->SetInput(glcmGenerator->GetOutput());
          featureCalc->Update();

          //TBD - to debug and compare vs GLRLM
          //std::cout << "[DEBUG] FeatureExtraction.hxx - CalculateGLCM - GLCM Matrix: Offset: " << this->m_offsets->at(i) << "\n";

          //std::cout << "\tindex\t|\t|\tfrenquency" << std::endl;
          //auto temp = glcmGenerator->GetOutput();
          //for (auto iter = temp->Begin(); iter != temp->End(); ++iter)
          //{
          //  std::cout << "\tGLCM Measurement vectors = " << iter.GetMeasurementVector()
          //    << "; Frequency = " << iter.GetFrequency() << std::endl;
          //}

          //TBD - to debug and compare vs GLRLM

          contrast += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Inertia));
          correl += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Correlation));
          ener += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Energy));
          homo += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment));
          entro += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Entropy));
          clustershade += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
          clusterprominance += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
          autocorr += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation));
        }

        contrast = contrast / this->m_offsets->size();
        correl = correl / this->m_offsets->size();
        ener = ener / this->m_offsets->size();
        homo = homo / this->m_offsets->size();
        entro = entro / this->m_offsets->size();
        clusterprominance = clusterprominance / this->m_offsets->size();
        clustershade = clustershade / this->m_offsets->size();
        autocorr = autocorr / this->m_offsets->size();

        //2019-05-09 - For Future Reference : correlation and autocorrelation values were not matching with IBSI, and are not being outputted to the users for now
        // 2020-03-05: these features have been re-enalbed for the NatSci paper

        this->m_features["Energy"] = ener;
        this->m_features["Entropy"] = entro;
        this->m_features["Homogeneity"] = homo; // also called "inverse difference moment"
        this->m_features["Contrast"] = contrast; // also called "inertia"
        this->m_features["ClusterShade"] = clustershade;
        this->m_features["ClusterProminence"] = clusterprominance;
        this->m_features["Correlation"] = correl;
        this->m_features["AutoCorrelation"] = autocorr; // called "haralick"
      }
      else if ((m_offsetSelector == "ITKDefault") || (m_offsetSelector == "Combined"))
      {
        typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
        glcmGenerator->SetNumberOfBinsPerAxis(this->m_Bins); //reasonable number of bins
        if (this->m_histogramBinningType != HistogramBinningType::Equal)
        {
          glcmGenerator->SetPixelValueMinMax(this->m_minimum, this->m_maximum);
        }
        glcmGenerator->SetMaskImage(mask_wrap);
        glcmGenerator->SetInput(image_wrap);
        auto featureCalc = Hist2FeaturesType::New();

        glcmGenerator->SetOffsets(this->m_offsets);
        glcmGenerator->Update();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        //2019-05-09 - For Future Reference : correlation and autocorrelation values were not matching with IBSI, and are not being outputted to the users for now

        this->m_features["Energy"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Energy));
        this->m_features["Entropy"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Entropy));
        this->m_features["Homogeneity"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment)); // also called "difference moment"
        this->m_features["Contrast"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Inertia)); // also called "inertia"
        this->m_features["ClusterShade"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
        this->m_features["ClusterProminence"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
        this->m_features["Correlation"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Correlation));
        this->m_features["AutoCorrelation"] = static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation)); // called "haralick"
      }
      else
      {
        typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
        glcmGenerator->SetNumberOfBinsPerAxis(this->m_Bins); //reasonable number of bins
        if (this->m_histogramBinningType != HistogramBinningType::Equal)
        {
          glcmGenerator->SetPixelValueMinMax(this->m_minimum, this->m_maximum);
        }
        glcmGenerator->SetMaskImage(mask_wrap);
        glcmGenerator->SetInput(image_wrap);
        auto featureCalc = Hist2FeaturesType::New();

        for (size_t i = 0; i < this->m_offsets->size(); i++)
        {
          glcmGenerator->SetOffset(this->m_offsets->at(i));
          glcmGenerator->Update();
          featureCalc->SetInput(glcmGenerator->GetOutput());
          featureCalc->Update();

        //2019-05-09 - For Future Reference : correlation and autocorrelation values were not matching with IBSI, and are not being outputted to the users for now

          auto tempStr = "_Offset_" + std::to_string(i);
          this->m_features[std::string("Energy") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Energy);
          this->m_features[std::string("Entropy") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Entropy);
          this->m_features[std::string("Homogeneity") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment); // also called "difference moment"
          this->m_features[std::string("Contrast") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Inertia); // also called "inertia"
          this->m_features[std::string("ClusterShade") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::ClusterShade);
          this->m_features[std::string("ClusterProminence") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence);
          this->m_features[std::string("Correlation") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::Correlation);
          this->m_features[std::string("AutoCorrelation") + tempStr] = featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation); // called "haralick"
        }
      }

      // TODO: Sung to add his GLCM extraction code here
      //this->m_features[std::string("Correlation]) + "_Sung"] = 0;

      this->m_algorithmDone = true;
    }
  }

private:
  std::string m_offsetSelector; //! type of offset selection
};
