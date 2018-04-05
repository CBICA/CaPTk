/**
\file FeatureExtraction.h

This file holds the declaration of the class FeatureExtraction.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html
*/
#pragma once

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
#include <sstream>
#include <tuple>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkScalarImageToRunLengthFeaturesFilter.h"
#include "itkHistogramToTextureFeaturesFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkVectorImage.h"
#include "itkScalarImageToCooccurrenceMatrixFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkVectorContainer.h"
#include "itkStatisticsImageFilter.h"
#include "itkImageToHistogramFilter.h"
#include "itkMaskedImageToHistogramFilter.h"
#include "itkShapeLabelObject.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageDuplicator.h"
#include "itkExtractImageFilter.h"

#include "itkEnhancedScalarImageToNeighbourhoodGreyLevelDifferenceFeaturesFilter.h"
#include "itkEnhancedScalarImageToSizeZoneFeaturesFilter.h"
#include "LBPFeatures2D.h"

#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"
#include "cbicaStatistics.h"

//#include "CAPTk.h"
#include "FeatureMap.h"

using ImageType2D = itk::Image< float, 2 >;
using ImageTypeFloat3D = itk::Image< float, 3 >;

//! FeatureType is the main datatype enclosing all the information about features
using FeatureType = std::map < // MAP_1: overall control map
  std::string, // feature family (intensity or volumetric or morphologic...)
  std::tuple <
  bool, // true if it the feature(s) has/have been selected manually (either via checkboxes in UI or '-f' via CLI)
  std::map< // MAP_2: this map structures has all the parameters and their values parsed from param_default file for a particular feature family
  std::string, // parameter name (for instance, "Bins")
  std::map< // MAP_3: collection of properties of the specific parameter defined above
  std::string, // name of a property (for instance Type/Range/Default)
  std::string >// corresponding value of the property name (if above is Type, this is INT and so on); close of MAP_3
  >, // close of MAP_2
  std::string, // modality name, for instance T1, T2, ...
  std::string, // roi label name, for instance ED, NCR, ...
  std::map< // MAP_4 Individual features and the ouput value mapped
  std::string, // feature name, for instance Contrast, Correlation, ...
  double> // actual output value of feature specified above; close of MAP_3
  > // close of tuple
>; // close of MAP_1

using FeatureBaseType = std::vector < std::tuple< std::string, std::string, std::string, double > >;



enum Params
{
  Dimension, Axis, Radius, Neighborhood, Bins, Directions, Offset, Range, Lattice
};
static const char ParamsString[Lattice + 1][15] =
{
  "Dimension", "Axis", "Radius", "Neighborhood", "Bins", "Directions", "Offset", "Range", "Lattice"
};
//enum FeatureGroup
//{
//  Statistics,
//  Miscellaneous,
//  Histogram,
//  Textural
//};
//
//static const char FeatureGroupString[Textural + 1][15] = { "Statistics", " Miscellaneous", "Histogram", "Textural" };

enum FeatureFamily
{
  Intensity,
  Histogram,
  Volumetric,
  Morphologic,
  GLCM,
  GLRLM,
  GLSZM,
  NGTDM,
  LBP
  //Gabor, 
  //Laws,
  //FractalDimension,   // to be added when testing completes
  //PowerSpectrun

};

static const char FeatureFamilyString[LBP + 1][15] = { "Intensity", "Histogram", "Volumetric", "Morphologic", "GLCM", "GLRLM", "GLSZM", "NGTDM", "LBP" };

enum IndividualFeatures
{
  Minimum, Maximum, Variance, Skewness, Kurtosis, STD, Correlation, Energy, Entropy, Contrast, SumAverage, GLCMVariance, GLCMMean, Homogenety, Clustershade, Clusterprominence, Autocorrelation,
  SRE, LRE, Run_GLN, RLN, LGRE, HGRE, SRLGE, SRHGE, LRLGE, LRHGE,
  SZE, LZE, GLN, ZSN, LGZE, HGZE, SZLGE, SZHGE, LZLGE, LZHGE, ZP, GLV, SZV, ZE, Coarsness, Busyness, Complexity, NeighbourContrast, Strength,
  Area, Pixels, Eccentricity, Elongation, Perimeter, Flatness, Roundness, Volume
};


static const char IndividualFeaturesString[Volume + 1][25] = { "Minimum", "Maximum", "Variance", "Skewness", "Kurtosis", "STD", "Correlation", "Energy", "Entropy",
"Contrast", "Mean", "GLCMVariance", "GLCMMean", "Homogeneity", "ClusterShade", "ClusterProminence", "Autocorrelation", "SRE", "LRE", "Run_GLN", "RLN", "LGRE", "HGRE", "SRLGE", "SRHGE", "LRLGE", "LRHGE",
"SZE", "LZE", "GLN", "ZSN", "LGZE", "HGZE", "SZLGE", "SZHGE", "LZLGE", "LZHGE", "ZP", "GLV", "SZV", "ZE", "Coarseness", "Busyness", "Complexity", "NeighbourContrast", "Strength",
"Area", "Pixels", "Eccentricity", "Elongation", "Perimeter", "Flatness", "Roundness", "Volume" };

/**
\brief FeatureExtraction Class -The main class structure enclosing all the feature calculations functions.
\ templated over Imagetype and ImageType::PixelType

*/
template< class TImageType = itk::Image< float, 3 >, class T = float>
class FeatureExtraction
{

public:
  using FeatureextractionImageType = typename TImageType::Pointer;
  using Image2CoOccuranceType = itk::Statistics::ScalarImageToCooccurrenceMatrixFilter < TImageType >;
  //typedef itk::Statistics::ScalarImageToCooccurrenceMatrixFilter<TImageType>Image2CoOccuranceType;
  using HistogramType = typename Image2CoOccuranceType::HistogramType;
  using Hist2FeaturesType = itk::Statistics::HistogramToTextureFeaturesFilter<HistogramType>;
  using OffsetType = typename TImageType::OffsetType;
  using Offsets = OffsetType; // TBD: replace 'Offsets' data type with 'OffsetType'
  using OffsetVector = itk::VectorContainer< unsigned char, OffsetType >;

  //Constructor
  FeatureExtraction()
  {
    m_cancel = false;
  }

  //! Default Destructor
  ~FeatureExtraction()
  {
    // no need to deallocate image pointers since they are smart, get it?
  }

  /**
  \brief Update function  -The call to various feature functions is done from here. This class also call the write ouput function.
  */
  void Update()
  {
    //    Check if input mask is not null
    auto minMaxCal = itk::MinimumMaximumImageCalculator< TImageType >::New();
    minMaxCal->SetImage(m_Mask);
    minMaxCal->ComputeMaximum();
    if (minMaxCal->GetMaximum() == 0)
    {
      std::string errorString = "Mask hasn't been initialized";

      auto exeName = cbica::getExecutableName();
      std::transform(exeName.begin(), exeName.end(), exeName.begin(), ::tolower);
      ////ShowErrorMessage("exeName = " + exeName);

      if (exeName.find("captk") != std::string::npos) // TBD this needs a better check than simply "captk", preferably related to qt
      {
        //  ShowErrorMessage(errorString);
        //  return m_OutputfeatureVector;
      }
      else
      {
        std::cerr << errorString << "\n";
        exit(EXIT_FAILURE);
      }
    }


    for (int i = 0; i < m_inputImages.size(); i++)
    {
      /* Iteration over requested number of roi*/
      for (int j = 0; j < m_labelname.size(); j++)
      {
        //FeatureBaseType  tempFeatureVec;
        auto selected_mask = get_selected_mask(m_inputImages[i], m_Mask, m_roi[j]);

        /* Intensity features are calculated irrespective of selection in the param file*/
        //bool intensityFeaturesCalculated = false;
        //if (!intensityFeaturesCalculated/*m_Features.find(FeatureFamilyString[Intensity]) != m_Features.end()*/)
        {
          auto temp = m_Features.find(FeatureFamilyString[Intensity]);
          std::get<2>(temp->second) = m_modality[i];
          std::get<3>(temp->second) = m_labelname[j];
          Intensity_features(m_inputImages[i], selected_mask, std::get<4>(temp->second));
          //intensityFeaturesCalculated = true;
        }


        /*Iterate over m_feture data structure */
        for (auto const &mapiterator : m_Features)
        {
          SetFeatureParam(mapiterator.first);
          if (mapiterator.first == FeatureFamilyString[Morphologic])
          {
            auto temp = m_Features.find(FeatureFamilyString[Morphologic]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              /* this dimensionality reduction applies only to shape and Volumetric features */
              if (TImageType::ImageDimension == 3)
              {
                if (m_Dimension == 2) // extracts slice with maximum area along the specified axis
                {
                  auto selected_axis_image = getSelectedslice(selected_mask, m_Axis);
                  ShapeFeatures<ImageType2D>(selected_axis_image, std::get<4>(temp->second));
                }
                else
                {
                  ShapeFeatures<TImageType>(selected_mask, std::get<4>(temp->second));
                }
              }
              else
              {
                ShapeFeatures<TImageType>(selected_mask, std::get<4>(temp->second));
              }
            }
          }

          /* Intensity features are calculated irrespective of selection in the param file*/
          if (mapiterator.first == FeatureFamilyString[Histogram])
          {
            auto temp = m_Features.find(FeatureFamilyString[Histogram]);
            if (std::get<0>(temp->second))
            {
              //auto local_map = std::get<1>(temp->second);
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];
              HistogramFeatures(m_inputImages[i], selected_mask, std::get<4>(temp->second));
            }
          }

          if (mapiterator.first == FeatureFamilyString[Volumetric])
          {
            auto temp = m_Features.find(FeatureFamilyString[Volumetric]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              /* this dimensionality reduction applies only to shape and Volumetric features */
              if (TImageType::ImageDimension == 3)
              {
                if (m_Dimension == 2)
                {
                  ImageType2D::Pointer selected_axis_image = getSelectedslice(selected_mask, m_Axis);
                  VolumetricFeatures<ImageType2D>(selected_axis_image, std::get<4>(temp->second));
                }
                else
                {
                  VolumetricFeatures<TImageType>(selected_mask, std::get<4>(temp->second));
                }
              }
              else
              {
                VolumetricFeatures<TImageType>(selected_mask, std::get<4>(temp->second));
              }
            }
          }

          if (mapiterator.first == FeatureFamilyString[GLRLM])
          {
            auto temp = m_Features.find(FeatureFamilyString[GLRLM]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              auto offsets = GetOffsetVector(m_Radius, m_Direction);

              calculateRunLength(m_inputImages[i], selected_mask, offsets, std::get<4>(temp->second));
            }
          }

          if (mapiterator.first == FeatureFamilyString[GLCM])
          {
            auto temp = m_Features.find(FeatureFamilyString[GLCM]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              auto offsets = GetOffsetVector(m_Radius, m_Direction);

              calculateGLCM(m_inputImages[i], selected_mask, offsets, std::get<4>(temp->second));
            }
          }

          if (mapiterator.first == FeatureFamilyString[GLSZM])
          {
            auto temp = m_Features.find(FeatureFamilyString[GLSZM]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              auto offsets = GetOffsetVector(m_Radius, m_Direction);

              //CalculateGrayLevelSizeZoneFeatures(m_inputImages[i], selected_mask, std::get<4>(temp->second));
            }
          }

          if (mapiterator.first == FeatureFamilyString[NGTDM])
          {
            auto temp = m_Features.find(FeatureFamilyString[NGTDM]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              auto offsets = GetOffsetVector(m_Radius, m_Direction);

              //CalculateGrayLevelNeighbourhoodGreyLevelDifferenceFeatures(m_inputImages[i], selected_mask, offsets, std::get<4>(temp->second));
            }
          }

          if (TImageType::ImageDimension == 3)
          {
            if (mapiterator.first == FeatureFamilyString[LBP])
            {
              auto temp = m_Features.find(FeatureFamilyString[LBP]);

              //m_nonZeroIndeces;
              std::set<int> x_value, y_value, z_value;
              for (unsigned int k = 0; k < m_nonZeroIndeces.size(); k++)
              {
                x_value.insert(m_nonZeroIndeces[k][0]);
                y_value.insert(m_nonZeroIndeces[k][1]);
                z_value.insert(m_nonZeroIndeces[k][2]);
              }

              const ImageTypeFloat3D::SizeType  size = { { x_value.size(), y_value.size(), z_value.size() } };
              const ImageTypeFloat3D::IndexType start = { { *x_value.begin(), *y_value.begin(), *z_value.begin() } };
              //Pad image before calculating LBP

              typename TImageType::SizeType lowerExtendRegion;
              lowerExtendRegion.Fill(m_Radius);

              typename TImageType::SizeType upperExtendRegion;
              upperExtendRegion.Fill(m_Radius);

              typename TImageType::PixelType constantPixel = 0;

              typename itk::ConstantPadImageFilter < TImageType, TImageType >::Pointer padFilter = itk::ConstantPadImageFilter < TImageType, TImageType >::New();
              padFilter->SetInput(selected_mask);
              //padFilter->SetPadBound(outputRegion); // Calls SetPadLowerBound(region) and SetPadUpperBound(region)
              padFilter->SetPadLowerBound(lowerExtendRegion);
              padFilter->SetPadUpperBound(upperExtendRegion);
              padFilter->SetConstant(constantPixel);
              padFilter->Update();
              typename TImageType::Pointer  lbproi = padFilter->GetOutput();
              const typename  TImageType::SizeType image_size = selected_mask->GetLargestPossibleRegion().GetSize();
              int x1 = -1;  int y1 = -1; int z1 = -1;
              for (unsigned int z = start[2] - m_Radius; z < start[2] + size[2] + m_Radius; z1++, z++)
              {
                y1 = -1;
                for (unsigned int y = start[1] - m_Radius; y < start[1] + m_Radius + size[1]; y1++, y++)
                {
                  x1 = -1;
                  for (unsigned int x = start[0] - m_Radius; x < start[0] + m_Radius + size[0]; x1++, x++)
                  {
                    //TImageType::IndexType ind1 = { x, y, z };
                    //TImageType::IndexType ind2 = { x1, y1, z1 };
                    //float  pixelValue;
                    //if (x < image_size[0] && y < image_size[1] && z < image_size[2])
                    //{
                    //  pixelValue = selected_mask->GetPixel(ind1);
                    //  lbproi->SetPixel(ind2, pixelValue);
                    //}
                    //else
                    //{
                    //  pixelValue = 0;  // x component
                    //  lbproi->SetPixel(ind2, pixelValue);
                    //}
                  }
                }
              }

              LBPFeatures lbpfeatures;

              lbpfeatures.calculateLBP<TImageType>(m_inputImages[i], selected_mask, lbproi, m_Radius, m_neighborhood, m_modality[i], std::to_string(m_roi[j]), std::get<4>(temp->second));
            }
          }

        } // end of feature iteration      
        m_OutputfeatureVector.push_back(m_Features);
      }/* End of iteration over roi*/


    } /*End of image iteration loop*/

    writeFeatureList(m_File, m_OutputfeatureVector/*, "csv"*/);
  }


  /**
  \brief SetFeatureParam function populates the parameters for the selected feature.

  The different parameter values are set in member variables before the calculations take place

  \param selected_feature Denotes a feature name
  */
  void SetFeatureParam(std::string selected_feature)
  {
    auto featurefamily = m_Features.find(selected_feature)->second;
    auto p = std::get<1>(featurefamily);


    if (p.find(ParamsString[Dimension]) != p.end())
    {
      auto  pp = p[ParamsString[Dimension]];
      m_Dimension = stoi(pp["Value"]);
    }
    if (p.find(ParamsString[Axis]) != p.end())
    {
      auto  pp = p[ParamsString[Axis]];
      m_Axis = pp["Value"];
    }
    if (p.find(ParamsString[Bins]) != p.end())
    {
      auto  pp = p[ParamsString[Bins]];
      m_Bins = stoi(pp["Value"]);
    }
    if (p.find(ParamsString[Radius]) != p.end())
    {
      auto  pp = p[ParamsString[Radius]];
      m_Radius = stoi(pp["Value"]);
    }
    if (p.find(ParamsString[Neighborhood]) != p.end())
    {
      auto  pp = p[ParamsString[Neighborhood]];
      m_neighborhood = stoi(pp["Value"]);
    }
    if (p.find(ParamsString[Offset]) != p.end())
    {
      auto  pp = p[ParamsString[Offset]];
      m_offsetSelect = pp["Value"];
    }
    if (p.find(ParamsString[Range]) != p.end())
    {
      auto  pp = p[ParamsString[Range]];
      m_Range = std::stoi(pp["Value"]);
    }
    if (p.find(ParamsString[Directions]) != p.end())
    {
      auto  pp = p[ParamsString[Directions]];
      m_Direction = std::stoi(pp["Value"]);
    }

  }

  /**
  \brief SetSelectedLabels function populates the member varibale m_roi and m_labelname.

  \param roi is the mask value we are interested in calculating features from.
  \param roi_labels is the corresponding string name of the provided roi numbers.

  */
  void SetSelectedLabels(std::string roi, std::string roi_labels)
  {
    if (!roi.empty())
    {
      std::vector<std::string> tempstr = cbica::stringSplit(roi, "|");
      for (size_t i = 0; i< tempstr.size(); i++)
      {
        m_roi.push_back(std::stoi(tempstr[i]));
      }
    }

    if (!roi_labels.empty())
    {
      m_labelname = cbica::stringSplit(roi_labels, "|");
    }


    if (m_labelname.size() != m_roi.size())
    {
      std::string errorString = "The selected roi and the provided roi labels doesn't match";

      auto exeName = cbica::getExecutableName();
      std::transform(exeName.begin(), exeName.end(), exeName.begin(), ::tolower);

      if (exeName.find("captk") != std::string::npos) // TBD this needs a better check than simply "captk", preferably related to qt
      {
        //ShowErrorMessage(errorString);
        //        return featurevec;
      }
      else
      {
        std::cerr << errorString << "\n";
        exit(EXIT_FAILURE);
      }
    }

  }

  /**
  \brief SetRequestedFetures function populates the main data structure m_features with the requested features and updates the bool value related to the features.

  \param filename is the default param file name provided .
  \param selected_features is the command line option given to user to pass or select required features from available ones.

  */
  void SetRequestedFetures(std::string filename, std::string selected_features)
  {
    //parse through feature param file and create a map  of features and params in file
    // make all features available in param file as true 

    std::string parsed;

    //Now check for user input if any feature has been de-selected and make it false
    if (!selected_features.empty())
    {

      std::stringstream feature_stringstream(selected_features);
      while (getline(feature_stringstream, parsed, '|'))
      {
        m_featureparams = FeatureMap(filename, parsed).getFeatureMap();
        m_Features[parsed] = std::make_tuple(true, m_featureparams, parsed, parsed, m_output);
      }
    }
    else
    {
      std::string selected_features_wrap = FeatureFamilyString[0];
      for (size_t i = 0; i < LBP; i++)
      {
        m_featureparams = FeatureMap(filename, FeatureFamilyString[i]).getFeatureMap();

        m_Features[FeatureFamilyString[i]] = std::make_tuple(!m_featureparams.empty(), m_featureparams, FeatureFamilyString[i], FeatureFamilyString[i], m_output);
      }
    }
  }

  /**
  \brief SetRequestedFetures function populates the main data structure m_features with the requested features and updates the bool value related to the features.

  This overloaded function is for the use of UI

  \param filename is the default param file name provided
  \param selected_features is the map of feature name as key and bool as the value
  */
  void SetRequestedFetures(std::string filename, std::map<std::string, bool> selected_features)
  {
    //parse through feature param file and create a map  of features and params in file
    // make all features available in param file as true 

    std::string parsed;

    //Now check for user input if any feature has been de-selected and make it false
    if (!selected_features.empty())
    {
      for (const auto &f : selected_features)
      {
        m_featureparams = FeatureMap(filename, f.first).getFeatureMap();
        m_Features[f.first] = std::make_tuple(f.second, m_featureparams, f.first, f.first, m_output);
      }
    }
  }

  /**
  \brief Populates the feature structure from the GUI
  */
  void SetRequestedFetures(std::map< std::string, std::vector< std::map<std::string, std::string> > >  featuresFromUI, std::map<std::string, bool> selected_features)
  {
    //parse through feature param file and create a map  of features and params in file
    // make all features available in param file as true 

    std::string parsed;

    for (auto &currentFeature : featuresFromUI) // iterating over the feature families, i.e., GLCM, Intensity, ...
    {
      auto selectedFeatureFlagStruct = selected_features.find(currentFeature.first); // this is to check for user input in UI (if selected or not)

      std::map<std::string, std::map<std::string, std::string>> temp;
      for (size_t j = 0; j < currentFeature.second.size(); j++)
      {
        std::map< std::string, std::string > currentFeature_ParamsAndVals;
        std::string paramName;
        for (auto &currentFeature_Parameter : currentFeature.second[j]) // each parameter within the feature family
        {
          if (currentFeature_Parameter.first != "ParamName")
          {
            currentFeature_ParamsAndVals[currentFeature_Parameter.first] = currentFeature_Parameter.second;
          }
          else
          {
            paramName = currentFeature_Parameter.second;
          }
        }
        temp[paramName] = currentFeature_ParamsAndVals;
      }

      m_Features[currentFeature.first] = std::make_tuple(selectedFeatureFlagStruct->second, // whether the feature is to be extracted or not
        temp, // parameters and respective values
        currentFeature.first, currentFeature.first, // these are the modality and roi label names, which get overwritten with the correct values in the "Update" function
        m_output);
    }
  }

  /**
  \brief writeFeatureList function is to write the calculated feature values to the ouput file given.
  \param filename is the absolute path of the csv file name
  \param featurevec is a vector of base feature type , enclosing the output to each file
  \param outputfiletype is option to switch between .xml and .csv
  */
  void writeFeatureList(std::string Filename, std::vector< FeatureType> featurevec, std::string outputfiletype = "")
  {
    std::string path, base, extension;
    cbica::splitFileName(Filename, path, base, extension);

    int current_modality = 0;
    if (extension == ".csv")
    {
      std::ofstream myfile;
      myfile.open(Filename);

      for (int j = 0; j < featurevec.size(); j++)
      {
        FeatureType tempfeature = featurevec[j];
        for (auto const &ent1 : tempfeature)
        {
          auto temptuple = ent1.second;
          auto features = std::get<4>(temptuple);
          for (auto const &f : features)
          {
            auto str_f_second = std::to_string(f.second);
            myfile << ent1.first + "_" + std::get<2>(temptuple) +"_" + std::get<3>(temptuple) +"_" + f.first + "," + str_f_second + "\n";
          }
          //myfile << "\n";
        }
      }
      myfile.close();
    }
    /// TBD: do something for XML
  }

  //! this is to be used in case this class becomes long-running
  void cancel()
  {
    m_cancel = true;
  }


  /**
  \brief SetInputImages parses the input images name and populates the member variable m_inputImages
  \param images is vector of input images.
  \param modality is a string containing the respective image modalities separated by |
  */

  void SetInputImages(std::vector< typename TImageType::Pointer > images, std::string modality)
  {
    m_inputImages = images;

    if (!modality.empty())
    {
      m_modality = cbica::stringSplit(modality, "|");
    }
  }


  /**
  \brief SetInputImages parses the input images name and populates the member variable m_inputImages
  \param images is vector of input images.
  \param modality is a vector of string containing the respective image modalities.
  */
  void SetInputImages(std::vector< typename TImageType::Pointer > images, std::vector< std::string > modality)
  {
    if (images.size() != m_modality.size())
    {
      std::cout << "Number of Images and number of modalities are not same";
    }
    m_inputImages = images;
    m_modality = modality;
  }

  /**
  \brief SetMaskImage parses the input mask and populates the member variable m_Mask
  \param image is Image pointer.
  */
  void SetMaskImage(typename TImageType::Pointer image)
  {
    m_Mask = image;
  }

  /**
  \brief Setouputfilename sets m_ouputpath value
  \param filename is the absolute path of output file.
  */
  void Setouputfilename(std::string filename)
  {
    m_File = filename;
    std::string ext, basename;
    cbica::splitFileName(m_File, m_ouputpath, basename, ext);
  }

  /**
  \brief Get the ROIs that were selected for the current process as a vector of ints
  */
  std::vector< int > GetSelectedROIs()
  {
    return m_roi;
  }

private:
  /**
  \brief Calculates the OffsetVector::Pointer based on the provided radius and directions
  */
  typename OffsetVector::Pointer GetOffsetVector(int inputRadius, int inputDirections)
  {
    itk::Neighborhood< typename TImageType::PixelType, TImageType::ImageDimension > neighborhood;
    neighborhood.SetRadius(inputRadius);
    auto size = neighborhood.GetSize();
    auto directionsToCompute = 1;
    for (size_t sizeCount = 0; sizeCount < TImageType::ImageDimension; sizeCount++)
    {
      directionsToCompute *= size[sizeCount];
    }
    if (inputDirections < directionsToCompute)
    {
      directionsToCompute = inputDirections;
    }

    typename OffsetVector::Pointer offsets = OffsetVector::New();

    for (int d = 0; d < directionsToCompute; d++)
    {
      offsets->push_back(neighborhood.GetOffset(d));
    }

    return offsets;
  }

  /**
  \brief getSelectedslice enables selection of a particular slice and direction in which the features have to be generated.
  \param mask is vector of mask images.
  \param axis represents the slice selection is in x, y or z direction.
  */
  ImageType2D::Pointer getSelectedslice(typename TImageType::Pointer mask, std::string axis)
  {
    std::vector< typename ImageType2D::Pointer > maxImageSlices;

    maxImageSlices.resize(3);
    typename TImageType::SizeType originalSize = mask->GetLargestPossibleRegion().GetSize();
    auto maxVoxels = originalSize;
    maxVoxels.Fill(0);

    if (originalSize.Dimension == 3)
    {
      for (int dim = 0; dim < 3; dim++) // dimension-loop
      {
        maxVoxels[dim] = 0;
        for (size_t i = 0; i < originalSize[dim]; i++)
        {
          typename TImageType::RegionType desiredRegion;
          typename TImageType::SizeType desiredSize = originalSize;
          desiredSize[dim] = 0;
          typename TImageType::IndexType desiredIndex;
          desiredIndex.Fill(0);
          desiredIndex[dim] = i;

          desiredRegion.SetIndex(desiredIndex);
          desiredRegion.SetSize(desiredSize);
          auto extractor = itk::ExtractImageFilter< TImageType, ImageType2D >::New();
          extractor->SetInput(mask);
          extractor->SetDirectionCollapseToIdentity();
          extractor->SetExtractionRegion(desiredRegion);
          extractor->Update();

          itk::ImageRegionConstIterator< ImageType2D > iterator(extractor->GetOutput(), extractor->GetOutput()->GetLargestPossibleRegion());
          size_t currentNonZero = 0;
          for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
          {
            if (iterator.Get() != 0)
            {
              currentNonZero++;
            }
          }

          if (currentNonZero > maxVoxels[dim])
          {
            maxVoxels[dim] = currentNonZero;

            auto duplicator = itk::ImageDuplicator< ImageType2D >::New();
            duplicator->SetInputImage(extractor->GetOutput());
            duplicator->Update();
            maxImageSlices[dim] = duplicator->GetOutput();
          }
        }
      }
    }
    if (axis == "x")
    {
      return maxImageSlices[0];
    }
    else if (axis == "y")
    {
      return maxImageSlices[1];
    }
    else
    {
      return maxImageSlices[2];
    }
  }

  /**
  \brief Calculate selected labels

  \param inputTotalmask A mask image with N different labels
  \param a int specifying the label to be extracted from mask image
  \return A mask image containing only the specific label supplied
  */
  typename TImageType::Pointer get_selected_mask(const typename TImageType::Pointer inputImage, const typename TImageType::Pointer total_mask, int roi)
  {
    auto mask = cbica::CreateImage< TImageType >(total_mask);

    TConstIteratorType  totalMaskIterator(total_mask, total_mask->GetLargestPossibleRegion());
    TConstIteratorType  imageIterator(inputImage, inputImage->GetLargestPossibleRegion());
    TIteratorType  maskIterator(mask, mask->GetLargestPossibleRegion());
    imageIterator.GoToBegin();
    totalMaskIterator.GoToBegin();
    maskIterator.GoToBegin();
    while (!totalMaskIterator.IsAtEnd())
    {
      if (totalMaskIterator.Get() == roi)
      {
        maskIterator.Set(1);
        m_nonZeroPixels.push_back(imageIterator.Get());
        m_nonZeroIndeces.push_back(totalMaskIterator.GetIndex());
      }

      ++imageIterator;
      ++maskIterator;
      ++totalMaskIterator;
    }
    return mask;
  }

  /**
  \brief Calculate intensity features
  \param image - ITKImage pointer
  \param mask - ITKImage pointer
  \param featurevec - map of Individual feature name and their value.
  */
  void Intensity_features(typename TImageType::Pointer image, typename TImageType::Pointer mask, std::map<std::string, double> &featurevec)
  {
    std::vector<float>nonzero_pix;

    if (m_nonZeroPixels.empty())
    {
      m_errormsg = true;
      return;
    }

    m_statistics.SetInput(m_nonZeroPixels);

    featurevec[IndividualFeaturesString[Minimum]] = m_statistics.GetMinimum();
    featurevec[IndividualFeaturesString[Maximum]] = m_statistics.GetMaximum();
    featurevec[IndividualFeaturesString[SumAverage]] = m_statistics.GetMean();
    featurevec[IndividualFeaturesString[Variance]] = m_statistics.GetVariance();
    featurevec[IndividualFeaturesString[STD]] = m_statistics.GetStandardDeviation();
    featurevec[IndividualFeaturesString[Skewness]] = m_statistics.GetSkewness();
    featurevec[IndividualFeaturesString[Kurtosis]] = m_statistics.GetKurtosis();

    return;
  }

  /**
  \brief Calculate GLCM Features

  This function computes features showcases the texture of a given image. this function uses ITK ScalarImageToTextureFeaturesFilter.
  The texture features have proven to be useful in image classification for biological and medical imaging.
  This function computes the texture features of a specified ROI(region of interest) got through the mask provided as an input para.
  The features that are calculated are contrast, correlation, energy, entropy, haralick correlation, cluster shde etc.
  The user can set the number of bins if required. The deafult value used is 16.

  \param image The input Image on which to calculate the features
  \param mask The image specificying the roi
  \param feature A vector holding features of each offset direction
  \param featurevec - map of Individual feature name and their value.
  */
  void calculateGLCM(typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec)
  {
    double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;

    if (m_offsetSelect == "Average")
    {

      for (size_t i = 0; i < offset->size(); i++)
      {
        typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
        glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
        glcmGenerator->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
        glcmGenerator->SetMaskImage(mask);
        glcmGenerator->SetInput(image);
        auto featureCalc = Hist2FeaturesType::New();

        glcmGenerator->SetOffset(offset->at(i));
        glcmGenerator->Update();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        contrast += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Inertia));
        correl += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Correlation));
        ener += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Energy));
        homo += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment));
        entro += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::Entropy));
        clustershade += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
        clusterprominance += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
        autocorr += static_cast<double>(featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation));
      }

      contrast = contrast / offset->size();
      correl = correl / offset->size();
      ener = ener / offset->size();
      homo = homo / offset->size();
      entro = entro / offset->size();
      clusterprominance = clusterprominance / offset->size();
      clustershade = clustershade / offset->size();
      autocorr = autocorr / offset->size();

      featurevec[IndividualFeaturesString[Correlation]] = correl;
      featurevec[IndividualFeaturesString[Contrast]] = contrast;
      featurevec[IndividualFeaturesString[Entropy]] = entro;
      featurevec[IndividualFeaturesString[Homogenety]] = homo;
      featurevec[IndividualFeaturesString[Clustershade]] = clustershade;
      featurevec[IndividualFeaturesString[Clusterprominence]] = clusterprominance;
      featurevec[IndividualFeaturesString[Autocorrelation]] = autocorr;
      featurevec[IndividualFeaturesString[Energy]] = ener;
    }
    else
    {
      typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
      glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
      glcmGenerator->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
      glcmGenerator->SetMaskImage(mask);
      glcmGenerator->SetInput(image);
      auto featureCalc = Hist2FeaturesType::New();

      for (size_t i = 0; i < offset->size(); i++)
      {
        glcmGenerator->SetOffset(offset->at(i));
        glcmGenerator->Update();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        featurevec[std::string(IndividualFeaturesString[Correlation]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::Correlation);
        featurevec[std::string(IndividualFeaturesString[Energy]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::Energy);
        featurevec[std::string(IndividualFeaturesString[Contrast]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::Inertia);
        featurevec[std::string(IndividualFeaturesString[Entropy]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::Entropy);
        featurevec[std::string(IndividualFeaturesString[Homogenety]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment);
        featurevec[std::string(IndividualFeaturesString[Clustershade]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::ClusterShade);
        featurevec[std::string(IndividualFeaturesString[Clusterprominence]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence);
        featurevec[std::string(IndividualFeaturesString[Autocorrelation]) + "_Offset_" + std::to_string(i)] = featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation);
      }
    }

    // TODO: Sung to add his GLCM extraction code here
    //featurevec[std::string(IndividualFeaturesString[Correlation]) + "_Sung"] = 0;
  }

  /**
  \brief Calculate RunLength Features
  This class computes run length descriptions from an image.By default, run length features are computed for each spatial direction and then averaged afterward,
  so it is possible to access the standard deviations of the texture features

  \param image The input Image on which to calculate the features
  \param mask The image specificying the roi
  \param feature A vector holding features of each offset direction
  \param featurevec - map of Individual feature name and their value
  */
  void calculateRunLength(typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec)
  {
    using HistogramFrequencyContainerType = itk::Statistics::DenseFrequencyContainer2;
    using RunLengthFilterType = itk::Statistics::ScalarImageToRunLengthFeaturesFilter< TImageType, HistogramFrequencyContainerType >;
    using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
    using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;
    using InternalRunLengthFeatureName = typename RunLengthFeatures::RunLengthFeatureName;

    typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
    matrix_generator->SetMaskImage(mask);
    matrix_generator->SetInput(image);
    matrix_generator->SetInsidePixelValue(1);
    matrix_generator->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
    //removed SetDistanceValueMinMax
    //matrix_generator->SetDistanceValueMinMax(0, m_Range); // TOCHECK - why is this only between 0-4? P
    matrix_generator->SetNumberOfBinsPerAxis(m_Bins); // TOCHECK - needs to be statistically significant

    typename  RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();

    typename  RunLengthFilterType::FeatureNameVectorPointer requestedFeatures = RunLengthFilterType::FeatureNameVector::New();
    requestedFeatures->push_back(RunLengthFeatures::ShortRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::LongRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::GreyLevelNonuniformity);
    requestedFeatures->push_back(RunLengthFeatures::RunLengthNonuniformity);
    requestedFeatures->push_back(RunLengthFeatures::LowGreyLevelRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::HighGreyLevelRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::ShortRunLowGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::ShortRunHighGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::LongRunLowGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeatures::LongRunHighGreyLevelEmphasis);

    std::vector <std::string> featurename;
    featurename.push_back(IndividualFeaturesString[SRE]);
    featurename.push_back(IndividualFeaturesString[LRE]);
    featurename.push_back(IndividualFeaturesString[GLN]);
    featurename.push_back(IndividualFeaturesString[RLN]);
    featurename.push_back(IndividualFeaturesString[LGRE]);
    featurename.push_back(IndividualFeaturesString[HGRE]);
    featurename.push_back(IndividualFeaturesString[SRLGE]);
    featurename.push_back(IndividualFeaturesString[SRHGE]);
    featurename.push_back(IndividualFeaturesString[LRLGE]);
    featurename.push_back(IndividualFeaturesString[LRHGE]);

    typename  OffsetVector::ConstIterator offsetIt;
    size_t offsetNum = 0, featureNum = 0;

    if (m_offsetSelect == "Average")
    {
      std::vector<double> tempfeatures/*(requestedFeatures->size(), 0)*/;
      tempfeatures.resize(requestedFeatures->size());
      //matrix_generator->SetOffsets(offset);
      for (offsetIt = offset->Begin(); offsetIt != offset->End(); offsetIt++, offsetNum++)
      {
        matrix_generator->SetOffset(offsetIt.Value());
        matrix_generator->Update();

        runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
        runLengthMatrixCalculator->Update();
        typename RunLengthFilterType::FeatureNameVector::ConstIterator fnameIt;
        //std::ostringstream ss;
        featureNum = 0;
        for (fnameIt = requestedFeatures->Begin(); fnameIt != requestedFeatures->End(); ++fnameIt, featureNum++)
        {
          //ss << offsetNum;
          tempfeatures[featureNum] = tempfeatures[featureNum] + runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
          //featurevec[featurename[featureNum]] = runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
        }
      }
      for (size_t i = 0; i < tempfeatures.size(); i++)
      {
        tempfeatures[i] = tempfeatures[i] / offset->size();
        featurevec[featurename[i]] = tempfeatures[i];
      }
    }
    else
    {
      for (offsetIt = offset->Begin(); offsetIt != offset->End(); offsetIt++, offsetNum++)
      {
        matrix_generator->SetOffset(offsetIt.Value());
        matrix_generator->Update();

        runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
        runLengthMatrixCalculator->Update();
        typename RunLengthFilterType::FeatureNameVector::ConstIterator fnameIt;
        featureNum = 0;
        for (fnameIt = requestedFeatures->Begin(); fnameIt != requestedFeatures->End(); fnameIt++, featureNum++)
        {
          featurevec[featurename[featureNum] + "_Offset_" + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
        }
      }
    }

  }

  /**
  \brief Calculate ShapeFeatures
  This class computes different shape features related to the Entropy region provided.

  \param mask The image specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  template< class TImageTypeShape = TImageType >
  void ShapeFeatures(typename TImageTypeShape::Pointer mask, std::map<std::string, double>& featurevec)
  {
    if (m_Dimension == 3)
    {
      typedef  short LabelType;
      typedef itk::Image< LabelType, TImageTypeShape::ImageDimension > OutputImageType;
      typedef itk::ShapeLabelObject< LabelType, TImageTypeShape::ImageDimension > ShapeLabelObjectType;
      typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

      typedef itk::ConnectedComponentImageFilter < TImageTypeShape, OutputImageType > ConnectedComponentImageFilterType;
      typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType > I2LType;
      typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
      connected->SetInput(mask);
      connected->Update();

      typename I2LType::Pointer i2l = I2LType::New();
      i2l->SetInput(connected->GetOutput());
      i2l->SetComputePerimeter(true);
      i2l->Update();
      typename LabelMapType::Pointer labelMap = i2l->GetOutput();
      //std::cout << " has " << labelMap->GetNumberOfLabelObjects() << " labels." << std::endl;
      //typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(0);

      //std::cout << labelObject->GetPrincipalMoments() << labelObject->GetElongation() <<
      //labelObject->GetPerimeter() << labelObject->GetRoundness() << labelObject->GetFlatness();
      auto spacing = mask->GetSpacing();
      double voxvol = spacing[0] * spacing[1] * spacing[2];

      //double volume = labelObject->GetNumberOfPixels() * voxvol;

      int numbers = labelMap->GetNumberOfLabelObjects();
      std::vector <double> eccentricity1;
      std::vector <double> eccentricity2;
      std::vector <double> roundness;
      std::vector <double> flatness;
      std::vector <double> elongation;
      std::vector <double> perimeter;

      for (int i = 1; i < numbers; i++)
      {
        typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);

        if (labelObject->GetNumberOfPixels() < 20)
          continue;
        auto princComps = labelObject->GetPrincipalMoments();
        eccentricity1.push_back(sqrt(1 - std::pow((princComps[0] / princComps[2]), 2)));
        eccentricity2.push_back(sqrt(1 - std::pow((princComps[0] / princComps[1]), 2)));
        roundness.push_back(labelObject->GetRoundness());
        flatness.push_back(labelObject->GetFlatness());
        elongation.push_back(labelObject->GetElongation());
        perimeter.push_back(labelObject->GetPerimeter());

      }
      double mean_ecc1 = 0;
      double mean_ecc2 = 0;
      double mean_round = 0;
      double mean_flat = 0;
      double mean_perim = 0;
      double mean_elong = 0;
      for (unsigned int i = 0; i < eccentricity1.size(); i++)
      {
        mean_ecc1 = mean_ecc1 + eccentricity1[i];
        mean_ecc2 = mean_ecc2 + eccentricity2[i];
        mean_round = mean_round + roundness[i];
        mean_flat = mean_flat + flatness[i];
        mean_perim = mean_perim + perimeter[i];
        mean_elong = mean_elong + elongation[i];
      }

      featurevec[IndividualFeaturesString[Eccentricity]] = mean_ecc1 / eccentricity1.size();
      featurevec[IndividualFeaturesString[Elongation]] = mean_elong / elongation.size();
      featurevec[IndividualFeaturesString[Perimeter]] = mean_perim;
      featurevec[IndividualFeaturesString[Roundness]] = mean_round / roundness.size();
      featurevec[IndividualFeaturesString[Flatness]] = mean_flat / flatness.size();
    }
    if (m_Dimension == 2)
    {
      typedef  short LabelType;
      typedef itk::Image< LabelType, TImageTypeShape::ImageDimension > OutputImageType;
      typedef itk::ShapeLabelObject< LabelType, TImageTypeShape::ImageDimension > ShapeLabelObjectType;
      typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

      typedef itk::ConnectedComponentImageFilter < TImageTypeShape, OutputImageType > ConnectedComponentImageFilterType;
      typedef itk::LabelImageToShapeLabelMapFilter < OutputImageType, LabelMapType > I2LType;
      typename ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
      connected->SetInput(mask);
      connected->Update();

      typename I2LType::Pointer i2l = I2LType::New();
      i2l->SetInput(connected->GetOutput());
      i2l->SetComputePerimeter(true);
      i2l->Update();
      typename LabelMapType::Pointer labelMap = i2l->GetOutput();
      typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(0);
      double pixel_spacing = mask->GetSpacing().GetElement(0)*mask->GetSpacing().GetElement(1);
      double area = labelObject->GetNumberOfPixels() * pixel_spacing;

      int numbers = labelMap->GetNumberOfLabelObjects();
      std::vector <double> eccentricity1;
      std::vector <double> eccentricity2;
      std::vector <double> roundness;
      std::vector <double> flatness;
      std::vector <double> elongation;
      std::vector <double> perimeter;

      for (int i = 1; i < numbers; i++)
      {
        typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(i);

        if (labelObject->GetNumberOfPixels() < 20) // this is to ensure really small portions of an ROI do not get picked 
          continue;

        auto princComps = labelObject->GetPrincipalMoments();
        eccentricity1.push_back(sqrt(1 - std::pow((princComps[0] / princComps[2]), 2)));
        eccentricity2.push_back(sqrt(1 - std::pow((princComps[0] / princComps[1]), 2)));
        roundness.push_back(labelObject->GetRoundness());
        flatness.push_back(labelObject->GetFlatness());
        elongation.push_back(labelObject->GetElongation());
        perimeter.push_back(labelObject->GetPerimeter());
      }

      double mean_ecc1 = 0;
      double mean_ecc2 = 0;
      double mean_round = 0;
      double mean_flat = 0;
      double mean_perim = 0;
      double mean_elong = 0;
      for (unsigned int i = 0; i < eccentricity1.size(); i++)
      {
        mean_ecc1 = mean_ecc1 + eccentricity1[i];
        mean_ecc2 = mean_ecc2 + eccentricity2[i];
        mean_round = mean_round + roundness[i];
        mean_flat = mean_flat + flatness[i];
        mean_perim = mean_perim + perimeter[i];
        mean_elong = mean_elong + elongation[i];
      }

      featurevec[IndividualFeaturesString[Eccentricity]] = mean_ecc1 / eccentricity1.size();
      featurevec[IndividualFeaturesString[Elongation]] = mean_elong / elongation.size();
      featurevec[IndividualFeaturesString[Perimeter]] = mean_perim;
      featurevec[IndividualFeaturesString[Roundness]] = mean_round / roundness.size();
      featurevec[IndividualFeaturesString[Flatness]] = mean_flat / flatness.size();
    }
  }


  /**
  \brief Calculate VolumetricFeatures
  This class computes volumetric features related to the roi region provided.

  \param mask The image specificying the roi
  \param featurevec - map of Individual feature name and their value
  */
  template< class TImageTypeVolumetric = TImageType >
  void VolumetricFeatures(typename TImageTypeVolumetric::Pointer mask, std::map<std::string, double>& featurevec)
  {
    int count = 0;
    itk::ImageRegionIteratorWithIndex< TImageTypeVolumetric > interIt(mask, mask->GetLargestPossibleRegion());
    for (interIt.GoToBegin(); !interIt.IsAtEnd(); ++interIt)
    {
      if (interIt.Get() > 0)
        count++;
    }
    auto spacing = mask->GetSpacing();
    double voxvol = 1;
    for (size_t i = 0; i < TImageTypeVolumetric::ImageDimension; i++)
    {
      voxvol *= spacing[i];
    }
    double volume = count * voxvol;
    featurevec[IndividualFeaturesString[Pixels]] = count;

    if (TImageTypeVolumetric::ImageDimension == 3)
    {
      featurevec[IndividualFeaturesString[Volume]] = volume;
    }
    else if (TImageTypeVolumetric::ImageDimension == 2)
    {
      featurevec[IndividualFeaturesString[Area]] = volume; // this is actually the area since the computation has happened in 2D
    }
    else
    {
      // do nothing
    }
  }


  /**
  \brief Calculate HistogramFeatures
  This class computes Histogram feature  the bin frequency given a lower and upper bound and the bin number.
  \param image The input image
  \param mask The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void HistogramFeatures(typename TImageType::Pointer image, typename TImageType::Pointer mask, std::map< std::string, double > &featurevec)
  {
    const float maxRescaleVal = 1000;
    typedef itk::RescaleIntensityImageFilter< TImageType, TImageType > RescaleFilterType;
    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
    rescaleFilter->SetInput(image);
    rescaleFilter->SetOutputMinimum(0);
    rescaleFilter->SetOutputMaximum(maxRescaleVal);
    rescaleFilter->Update();

    std::vector<double> intensities;
    itk::ImageRegionConstIterator <TImageType> imageIt(rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion());
    itk::ImageRegionConstIterator <TImageType> maskIt(mask, mask->GetLargestPossibleRegion());
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
    for (double i = 0; i < m_Bins; i++)
    {
      std::vector<double> onerange;
      onerange.push_back(current_index);
      current_index = current_index + final_interval;
      onerange.push_back(current_index);
      Ranges.push_back(onerange);

      if (Ranges.size() == m_Bins)
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
      featurevec["Bin_" + std::to_string(j)] = finalBins[j];
      featurevec["BinEndIntensity_" + std::to_string(j)] = Ranges[j][1];
    }
  }


  /**
  \brief Calculate GLSZM features
  A Gray Level Size Zone(GLSZM) quantifies gray level zones in an image.A gray level zone is defined as a the number of connected voxels that share the same gray level intensity.

  \param image The input image
  \param mask The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGrayLevelSizeZoneFeatures(typename TImageType::Pointer itkImage, typename TImageType::Pointer maskImage, std::map<std::string, double>& featurevec)
  {
    //typedef TImageType ImageType;
    //typedef TImageType MaskType;
    typedef itk::Statistics::EnhancedScalarImageToSizeZoneFeaturesFilter< TImageType > FilterType;
    typedef itk::MinimumMaximumImageCalculator< TImageType > MinMaxComputerType;
    typedef typename FilterType::SizeZoneFeaturesFilterType TextureFilterType;

    typename FilterType::Pointer filter = FilterType::New();

    typename FilterType::OffsetVector::Pointer newOffset = FilterType::OffsetVector::New();
    auto oldOffsets = filter->GetOffsets();
    auto oldOffsetsIterator = oldOffsets->Begin();
    while (oldOffsetsIterator != oldOffsets->End())
    {
      bool continueOuterLoop = false;
      typename FilterType::OffsetType offset = oldOffsetsIterator->Value();
      for (size_t i = 0; i < TImageType::ImageDimension; ++i)
      {
        if (/*params.m_Direction == i + 2 &&*/ offset[i] != 0)
        {
          continueOuterLoop = true;
        }
      }
      /*if (params.m_Direction == 1)
      {
      offset[0] = 0;
      offset[1] = 0;
      offset[2] = 1;
      newOffset->push_back(offset);
      break;
      }*/

      oldOffsetsIterator++;
      if (continueOuterLoop)
        newOffset->push_back(offset);
    }
    filter->SetOffsets(newOffset);


    // All features are required
    typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();
    requestedFeatures->push_back(TextureFilterType::SmallZoneEmphasis);
    requestedFeatures->push_back(TextureFilterType::LargeZoneEmphasis);
    requestedFeatures->push_back(TextureFilterType::GreyLevelNonuniformity);
    requestedFeatures->push_back(TextureFilterType::SizeZoneNonuniformity);
    requestedFeatures->push_back(TextureFilterType::LowGreyLevelZoneEmphasis);
    requestedFeatures->push_back(TextureFilterType::HighGreyLevelZoneEmphasis);
    requestedFeatures->push_back(TextureFilterType::SmallZoneLowGreyLevelEmphasis);
    requestedFeatures->push_back(TextureFilterType::SmallZoneHighGreyLevelEmphasis);
    requestedFeatures->push_back(TextureFilterType::LargeZoneLowGreyLevelEmphasis);
    requestedFeatures->push_back(TextureFilterType::LargeZoneHighGreyLevelEmphasis);
    requestedFeatures->push_back(TextureFilterType::ZonePercentage);
    requestedFeatures->push_back(TextureFilterType::GreyLevelVariance);
    requestedFeatures->push_back(TextureFilterType::SizeZoneVariance);
    requestedFeatures->push_back(TextureFilterType::ZoneEntropy);

    filter->SetInput(itkImage);
    filter->SetMaskImage(maskImage);
    filter->SetRequestedFeatures(requestedFeatures);
    int rangeOfPixels = m_Range;
    /*if (rangeOfPixels < 2)
    rangeOfPixels = 256;

    if (params.m_UseCtRange)
    {
    filter->SetPixelValueMinMax((TPixel)(-1024.5),(TPixel)(3096.5));
    filter->SetNumberOfBinsPerAxis(3096.5+1024.5);
    } else*/
    {
      filter->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
      filter->SetNumberOfBinsPerAxis(m_Bins);
    }

    //filter->SetDistanceValueMinMax(0, rangeOfPixels);

    filter->Update();

    auto featureMeans = filter->GetFeatureMeans();
    auto featureStd = filter->GetFeatureStandardDeviations();

    //filter->Delete();

    std::ostringstream  ss;
    ss << rangeOfPixels;
    std::string strRange = ss.str();
    for (std::size_t i = 0; i < featureMeans->size(); ++i)
    {
      switch (i)
      {
      case TextureFilterType::SmallZoneEmphasis:
        featurevec[IndividualFeaturesString[ZE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::LargeZoneEmphasis:
        featurevec[IndividualFeaturesString[LZE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::GreyLevelNonuniformity:
        featurevec[IndividualFeaturesString[GLN]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::SizeZoneNonuniformity:
        featurevec[IndividualFeaturesString[ZSN]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::LowGreyLevelZoneEmphasis:
        featurevec[IndividualFeaturesString[LGZE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::HighGreyLevelZoneEmphasis:
        featurevec[IndividualFeaturesString[HGZE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::SmallZoneLowGreyLevelEmphasis:
        featurevec[IndividualFeaturesString[SZLGE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::SmallZoneHighGreyLevelEmphasis:
        featurevec[IndividualFeaturesString[SZHGE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::LargeZoneLowGreyLevelEmphasis:
        featurevec[IndividualFeaturesString[LZLGE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::LargeZoneHighGreyLevelEmphasis:
        featurevec[IndividualFeaturesString[LZHGE]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::ZonePercentage:
        featurevec[IndividualFeaturesString[ZP]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::GreyLevelVariance:
        featurevec[IndividualFeaturesString[GLV]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::SizeZoneVariance:
        featurevec[IndividualFeaturesString[SZV]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::ZoneEntropy:
        featurevec[IndividualFeaturesString[ZE]] = featureMeans->ElementAt(i);
        break;
      default:
        break;
      }
    }
  }


  /**
  \brief Calculate NGTDM features

  \param image The input image
  \param mask The mask specificying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGrayLevelNeighbourhoodGreyLevelDifferenceFeatures(typename TImageType::Pointer itkImage,
    typename TImageType::Pointer maskImage, OffsetVector *offset, std::map<std::string, double>& featurevec)
  {
    //typedef TImageType MaskType;
    typedef itk::Statistics::EnhancedScalarImageToNeighbourhoodGreyLevelDifferenceFeaturesFilter< TImageType > FilterType;
    typedef itk::MinimumMaximumImageCalculator< TImageType > MinMaxComputerType;
    typedef typename FilterType::NeighbourhoodGreyLevelDifferenceFeaturesFilterType TextureFilterType;

    typename FilterType::Pointer filter = FilterType::New();

    typename FilterType::OffsetVector::Pointer newOffset = FilterType::OffsetVector::New();
    auto oldOffsets = filter->GetOffsets();
    auto oldOffsetsIterator = oldOffsets->Begin();
    while (oldOffsetsIterator != oldOffsets->End())
    {
      bool continueOuterLoop = false;
      typename FilterType::OffsetType offset = oldOffsetsIterator->Value();
      for (size_t i = 0; i < TImageType::ImageDimension; ++i)
      {
        if (/*params.m_Direction == i + 2 &&*/ offset[i] != 0)
        {
          continueOuterLoop = true;
        }
      }
      oldOffsetsIterator++;
      if (continueOuterLoop)
        newOffset->push_back(offset);
    }
    filter->SetOffsets(newOffset);


    // All features are required
    typename FilterType::FeatureNameVectorPointer requestedFeatures = FilterType::FeatureNameVector::New();
    requestedFeatures->push_back(TextureFilterType::Coarseness);
    requestedFeatures->push_back(TextureFilterType::Contrast);
    requestedFeatures->push_back(TextureFilterType::Busyness);
    requestedFeatures->push_back(TextureFilterType::Complexity);
    requestedFeatures->push_back(TextureFilterType::Strength);

    typename MinMaxComputerType::Pointer minMaxComputer = MinMaxComputerType::New();
    minMaxComputer->SetImage(itkImage);
    minMaxComputer->Compute();

    auto duplicator_image = itk::ImageDuplicator< TImageType >::New();
    auto duplicator_mask = itk::ImageDuplicator< TImageType >::New();
    duplicator_image->SetInputImage(itkImage);
    duplicator_image->Update();
    duplicator_mask->SetInputImage(maskImage);
    duplicator_mask->Update();

    auto testimage = duplicator_image->GetOutput();
    auto testmask = duplicator_mask->GetOutput();
    filter->SetInput(testimage);
    filter->SetMaskImage(testmask);
    filter->SetRequestedFeatures(requestedFeatures);
    int rangeOfPixels = m_Range;
    {
      filter->SetPixelValueMinMax(m_statistics.GetMinimum(), m_statistics.GetMaximum());
      filter->SetNumberOfBinsPerAxis(m_Bins);
    }

    //filter->SetDistanceValueMinMax(0, rangeOfPixels);

    filter->Update();

    auto featureMeans = filter->GetFeatureMeans();
    auto featureStd = filter->GetFeatureStandardDeviations();

    //std::ostringstream  ss;
    //ss << rangeOfPixels;
    //std::string strRange = ss.str();
    for (std::size_t i = 0; i < featureMeans->size(); ++i)
    {
      switch (i)
      {
      case TextureFilterType::Coarseness:
        featurevec[IndividualFeaturesString[Coarsness]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::Contrast:
        featurevec[IndividualFeaturesString[NeighbourContrast]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::Busyness:
        featurevec[IndividualFeaturesString[Busyness]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::Complexity:
        featurevec[IndividualFeaturesString[Complexity]] = featureMeans->ElementAt(i);
        break;
      case TextureFilterType::Strength:
        featurevec[IndividualFeaturesString[Strength]] = featureMeans->ElementAt(i);
        break;
      default:
        break;
      }
    }
  }

  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;
  cbica::Statistics< typename TImageType::PixelType > m_statistics;
  std::vector< float > m_nonZeroPixels;
  std::vector< typename TImageType::IndexType > m_nonZeroIndeces;
  typename TImageType::Pointer m_Mask;
  std::vector< int > m_roi;
  //private: // TBD: this needs to be private in a proper OOP architecture
  std::vector <typename TImageType::Pointer> m_inputImages;
  std::vector<std::string> m_modality; // size corresponds to m_inputImages.size()

  std::vector<std::string> m_labels; // ROI labels
  std::vector<std::string> m_labelname; // ROI label names and size should correspond to m_labels.size()
  std::map<std::string, bool> m_selectedFeatures;
  std::map<std::string, std::map<std::string, std::string>> m_featureparams;
  FeatureType m_Features;
  std::map < std::string, double> m_output;
  std::vector< FeatureType> m_OutputfeatureVector;
  std::string m_File; // output file
  std::string  m_ouputpath; // this is the output directory
  bool m_cancel;

  bool m_errormsg = false;
  int m_Radius = 0, m_Bins = 0, m_Dimension = 0, m_Direction = 0, m_Range = 0, m_neighborhood = 0;
  std::string m_Axis, m_offsetSelect;

};
