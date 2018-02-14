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
#include "cbicaLogging.h"

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
  Minimum, Maximum, Variance, Skewness, Kurtosis, STD, Correlation, Energy, Entrophy, Contrast, SumAverage, GLCMVariance, GLCMMean, Homogenety, Clustershade, Clusterprominence, Autocorrelation,
  SRE, LRE, Run_GLN, RLN, LGRE, HGRE, SRLGE, SRHGE, LRLGE, LRHGE,
  SZE, LZE, GLN, ZSN, LGZE, HGZE, SZLGE, SZHGE, LZLGE, LZHGE, ZP, GLV, SZV, ZE, Coarsness, Busyness, Complexity, NeighbourContrast, Strength,
  Area, Pixels, Eccentricity, Elongation, Perimeter, Flatness, Roundness, Volume
};


static const char IndividualFeaturesString[Volume + 1][25] = { "Minimum", " Maximum", "Variance", "Skewness", "Kurtosis", "STD", "Correlation", "Energy", "Entrophy",
"Contrast", "SumAverage", "GLCMVariance", "GLCMMean", "Homogenity", "ClusterShade", "Clusterprominence", "Autocorrelation", "SRE", "LRE", "Run_GLN", "RLN", "LGRE", "HGRE", "SRLGE", "SRHGE", "LRLGE", "LRHGE",
"SZE", "LZE", "GLN", "ZSN", "LGZE", "HGZE", "SZLGE", "SZHGE", "LZLGE", "LZHGE", "ZP", "GLV", "SZV", "ZE", "Coarsness", "Busyness", "Complexity", "NeighbourContrast", "Strength",
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
        FeatureBaseType  tempFeatureVec;
        auto selected_mask = get_selected_mask(m_Mask, m_roi[j]);


        /* Intensity features are calculated irrespective of selection in the param file*/
        if (m_Features.find(FeatureFamilyString[Intensity]) != m_Features.end())
        {
          auto temp = m_Features.find(FeatureFamilyString[Intensity]);
          std::get<2>(temp->second) = m_modality[i];
          std::get<3>(temp->second) = m_labelname[j];
          Intensity_features(m_inputImages[i], selected_mask, std::get<4>(temp->second));
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
                if (m_Dimension == 2)
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

              itk::Neighborhood< typename TImageType::PixelType, TImageType::ImageDimension > neighborhood;
              neighborhood.SetRadius(m_Radius);

              typename OffsetVector::Pointer offsets = OffsetVector::New();

              for (int d = 0; d < m_Direction; d++)
              {
                offsets->push_back(neighborhood.GetOffset(d));
              }
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

              itk::Neighborhood< typename TImageType::PixelType, TImageType::ImageDimension > neighborhood;
              neighborhood.SetRadius(m_Radius);

              typename OffsetVector::Pointer offsets = OffsetVector::New();

              //// TBD: check this section of code; this is for the direction
              //if (m_Direction == 3)
              //{
              //  m_Direction = 6;
              //}
              //for (int d = 0; d < m_Direction; d++)
              //{
              //  offsets->push_back(neighborhood.GetOffset(d));
              //}
              ////

              ///// TBD: calculate the full set of offsets using this method
              //std::vector< itk::Offset< 3 > > vectorOfOffsets;
              //auto size = neighborhood.GetSize();
              //auto centerIdx = neighborhood.GetCenterNeighborhoodIndex();
              //for (size_t i = 0; i < TImageType::ImageDimension; i++)
              //{
              //  totalSize *= size[i];
              //}
              //for (size_t i = 0; i < size[0] * size[1] * size[2]; i++)
              //{
              //  // this is to ensure the middle offset (which is the center) isn't taken into account
              //  if (i != centerIdx)
              //  {
              //    offsets->emplace_back(neighborhood.GetOffset(i));
              //    //std::cout << test.GetOffset(i) << "\n";
              //  }
              //}
              /////

              for (int d = 0; d < neighborhood.GetCenterNeighborhoodIndex(); d++)
              {
                offsets->push_back(neighborhood.GetOffset(d));
                //offsets->push_back(neighborhood.GetOffset(6));
                //offsets->push_back(neighborhood.GetOffset(3));
                //offsets->push_back(neighborhood.GetOffset(0));
              }

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

              itk::Neighborhood< typename TImageType::PixelType, TImageType::ImageDimension > neighborhood;
              neighborhood.SetRadius(m_Radius);

              typename OffsetVector::Pointer offsets = OffsetVector::New();

              for (int d = 0; d < neighborhood.GetCenterNeighborhoodIndex(); d++)
              {
                offsets->push_back(neighborhood.GetOffset(d));
                //offsets->push_back(neighborhood.GetOffset(6));
                //offsets->push_back(neighborhood.GetOffset(3));
                //offsets->push_back(neighborhood.GetOffset(0));
              }
              CalculateGrayLevelSizeZoneFeatures(m_inputImages[i], selected_mask, std::get<4>(temp->second));
            }
          }

          if (mapiterator.first == FeatureFamilyString[NGTDM])
          {
            auto temp = m_Features.find(FeatureFamilyString[NGTDM]);
            if (std::get<0>(temp->second))
            {
              std::get<2>(temp->second) = m_modality[i];
              std::get<3>(temp->second) = m_labelname[j];

              itk::Neighborhood< typename TImageType::PixelType, TImageType::ImageDimension > neighborhood;
              neighborhood.SetRadius(m_Radius);

              typename OffsetVector::Pointer offsets = OffsetVector::New();

              for (int d = 0; d < neighborhood.GetCenterNeighborhoodIndex(); d++)
              {
                offsets->push_back(neighborhood.GetOffset(d));
                //offsets->push_back(neighborhood.GetOffset(6));
                //offsets->push_back(neighborhood.GetOffset(3));
                //offsets->push_back(neighborhood.GetOffset(0));
              }

              CalculateGrayLevelNeighbourhoodGreyLevelDifferenceFeatures(m_inputImages[i], selected_mask, offsets, std::get<4>(temp->second));
            }
          }
        } // end of feature iteration
      }/* End of iteration over roi*/

      m_OutputfeatureVector.push_back(m_Features);
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
      std::string errorString = "The selected roi and the provided roi lables doesn't match";

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
  }

  /**
  \brief SetRequestedFetures function populates the main data structure m_features with the requested features and updates the bool value related to the features.

  \param filename is the default param file name provided .
  \param selected_features is the map of feature name as key and bool as the value.
  This overloaded function is for the use of UI .

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
  \brief writeFeatureList function is to write the calculated feature values to the ouput file given.
  \param filename is the absolute path of the csv file name.
  \param featurevec is a vector of base feature type , enclosing the output to each file.
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
            myfile << ent1.first + "_" + std::get<2>(temptuple) + "_" + std::get<3>(temptuple) + "_" + f.first + "," + str_f_second + "\n";
          }
          myfile << "\n";
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
  \param modality is a string contating the respective image modalities seperated by |
  */

  void SetInputImages(std::vector< typename TImageType::Pointer > images, std::string modality)
  {
    m_inputImages = images;


    if (!modality.empty())
    {
      // TBD: replace 
      m_modality = cbica::stringSplit(modality, "|");
    }
  }


  /**
  \brief SetInputImages parses the input images name and populates the member variable m_inputImages
  \param images is vector of input images.
  \param modality is a vector of string contating the respective image modalities.
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
  \parma fileanme is the absolute path of ouput file.
  */
  void Setouputfilename(std::string filename)
  {
    m_File = filename;
    std::string ext, basename;
    cbica::splitFileName(m_File, m_ouputpath, basename, ext);

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

  \param inputTotalmask A mask image with N ifferent labels
  \param a int specifying the label to be extracted from mask image
  \return A mask image containing only the specific label supplied
  */
  typename TImageType::Pointer get_selected_mask(typename TImageType::Pointer total_mask, int roi)
  {
    typename  TImageType::Pointer mask = TImageType::New();
    mask->SetLargestPossibleRegion(total_mask->GetLargestPossibleRegion());
    mask->SetRequestedRegion(total_mask->GetRequestedRegion());
    mask->SetBufferedRegion(total_mask->GetBufferedRegion());
    mask->SetDirection(total_mask->GetDirection());
    mask->SetOrigin(total_mask->GetOrigin());
    mask->SetSpacing(total_mask->GetSpacing());
    mask->Allocate();
    mask->FillBuffer(0);

    typedef itk::ImageRegionIterator< TImageType > IteratorType;
    IteratorType  imageIterator(total_mask, total_mask->GetRequestedRegion());
    IteratorType  maskIterator(mask, mask->GetRequestedRegion());
    imageIterator.GoToBegin();
    while (!imageIterator.IsAtEnd())
    {
      if (imageIterator.Get() == roi)
      {
        maskIterator.SetIndex(imageIterator.GetIndex());
        maskIterator.Set(1);
      }

      ++imageIterator;
    }
    return mask;
  }

  /**
  \brief Calculate intensity features
  \param image - ITKImage pointer
  \param mask - ITKImage pointer
  \param featurevec - map of Individual feature name and thier value.
  */
  void Intensity_features(typename TImageType::Pointer image, typename TImageType::Pointer mask, std::map<std::string, double> &featurevec)
  {
    std::vector<float>nonzero_pix;
    typedef itk::ImageRegionConstIterator< TImageType > IteratorType;
    std::vector< typename TImageType::IndexType> index_vec;

    // TBD: replace this with a MaskFilter + MinMaxCalculator + cbica::Statistics by putting all non-zero values in a vector and using that to initialize cbica::Statistics class
    IteratorType inputIt(mask, mask->GetLargestPossibleRegion());
    for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    {
      if (inputIt.Get() != 0)
      {
        index_vec.push_back(inputIt.GetIndex());
        nonzero_pix.push_back(image->GetPixel(inputIt.GetIndex()));
      }
    }

    if (nonzero_pix.empty())
    {
      m_errormsg = true;
      return;
    }

    double min = *min_element(nonzero_pix.begin(), nonzero_pix.end());
    double max = *max_element(nonzero_pix.begin(), nonzero_pix.end());
    double sum = std::accumulate(nonzero_pix.begin(), nonzero_pix.end(), 0.0);
    double mean = sum / nonzero_pix.size();

    // TBD: use cbica::Statistics for all of these and ensure all these values and available throughout the class; take care of the different modalities 
    double stdev = getstd(nonzero_pix);
    double variance = FeatureExtraction::getvariance(nonzero_pix);
    double skew = FeatureExtraction::getskewness(nonzero_pix);
    double kurtosis = FeatureExtraction::getkurtosis(nonzero_pix);

    featurevec[IndividualFeaturesString[Maximum]] = max;
    featurevec[IndividualFeaturesString[SumAverage]] = mean;
    featurevec[IndividualFeaturesString[Variance]] = variance;
    featurevec[IndividualFeaturesString[STD]] = stdev;
    featurevec[IndividualFeaturesString[Skewness]] = skew;
    featurevec[IndividualFeaturesString[Kurtosis]] = kurtosis;

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
  \param featurevec - map of Individual feature name and thier value.


  */
  void calculateGLCM(typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec)
  {
    typedef itk::MinimumMaximumImageCalculator<TImageType> MinMaxComputerType;
    typename MinMaxComputerType::Pointer minMaxComputer = MinMaxComputerType::New();
    minMaxComputer->SetImage(image);
    minMaxComputer->Compute();

    typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
    std::vector<double> nonzero_pix;
    IteratorType inputIt(mask, mask->GetLargestPossibleRegion()), imageIt(image, image->GetLargestPossibleRegion());
    inputIt.GoToBegin();
    for (; !inputIt.IsAtEnd(); ++inputIt)
    {
      if (inputIt.Get() != 0)
      {
        imageIt.SetIndex(inputIt.GetIndex());
        nonzero_pix.push_back(imageIt.Get());
      }
      //if (mask->GetPixel(inputIt.GetIndex()) != 0)
      //{
      //  nonzero_pix.push_back(image->GetPixel(inputIt.GetIndex()));
      //}
    }
    // calculate mean and variance
    double sum = std::accumulate(nonzero_pix.begin(), nonzero_pix.end(), 0.0);
    double pixelMean = sum / nonzero_pix.size();
    double var = 0;
    for (size_t x = 0; x < nonzero_pix.size(); x++)
    {
      var += pow((nonzero_pix[x] - pixelMean), 2);
    }
    double pixelVariance = (var / nonzero_pix.size());

    // input image is getting rescaled/quantized to (0, m_Bins - 1) here
    typedef itk::RescaleIntensityImageFilter< TImageType, TImageType > RescaleFilterType;
    //typename   RescaleFilterType::Pointer filter = RescaleFilterType::New();
    //filter->SetInput(image);
    //filter->SetOutputMinimum(0);
    //filter->SetOutputMaximum(m_Bins - 1);
    //filter->Update();
    double contrast = 0, correl = 0, ener = 0, entro = 0, homo = 0, clustershade = 0, clusterprominance = 0, autocorr = 0;
    auto result = std::minmax_element(nonzero_pix.begin(), nonzero_pix.end());
    double m_min = *result.first;
    double m_max = *result.second;
    if (m_offsetSelect == "Average")
    {
      for (size_t i = 0; i < offset->size(); i++)
      {
        typename  Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
        glcmGenerator->SetOffset(offset->at(i));
        glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
        glcmGenerator->SetPixelValueMinMax(m_min, m_max);
        glcmGenerator->SetMaskImage(mask);
        glcmGenerator->SetInput(image);
        glcmGenerator->Update();
        auto featureCalc = Hist2FeaturesType::New();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        // cbica::WriteImage<TImageType>(image, "oriimge.nii.gz");
        //cbica::WriteImage<TImageType>(mask, "mask.nii.gz");

        //  std::cout << offset[i] << "\n";
        
        contrast = static_cast<double>(contrast + featureCalc->GetFeature(Hist2FeaturesType::Inertia));
        correl = static_cast<double>(correl + featureCalc->GetFeature(Hist2FeaturesType::Correlation));
        ener = static_cast<double>(ener + featureCalc->GetFeature(Hist2FeaturesType::Energy));
        homo = static_cast<double>(homo + featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment));
        entro = static_cast<double>(entro + featureCalc->GetFeature(Hist2FeaturesType::Entropy));
        clustershade = static_cast<double>(clustershade + featureCalc->GetFeature(Hist2FeaturesType::ClusterShade));
        clusterprominance = static_cast<double>(clusterprominance + featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence));
        autocorr = static_cast<double>(autocorr + featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation));
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
      featurevec[IndividualFeaturesString[Entrophy]] = entro;
      featurevec[IndividualFeaturesString[GLCMMean]] = pixelMean;
      featurevec[IndividualFeaturesString[GLCMVariance]] = pixelVariance;
      featurevec[IndividualFeaturesString[Homogenety]] = homo;
      featurevec[IndividualFeaturesString[Clustershade]] = clustershade;
      featurevec[IndividualFeaturesString[Clusterprominence]] = clusterprominance;
      featurevec[IndividualFeaturesString[Autocorrelation]] = autocorr;
      featurevec[IndividualFeaturesString[Energy]] = ener;
    }
    else
    {
      for (size_t i = 0; i < offset->size(); i++)
      {
        typename Image2CoOccuranceType::Pointer glcmGenerator = Image2CoOccuranceType::New();
        glcmGenerator->SetOffset(offset->at(i));
        glcmGenerator->SetNumberOfBinsPerAxis(m_Bins); //reasonable number of bins
        glcmGenerator->SetPixelValueMinMax(m_min, m_max); //for input UCHAR pixel type
        glcmGenerator->SetMaskImage(mask);
        glcmGenerator->SetInput(image);
        glcmGenerator->Update();
        auto featureCalc = Hist2FeaturesType::New();
        featureCalc->SetInput(glcmGenerator->GetOutput());
        featureCalc->Update();

        featurevec[IndividualFeaturesString[Correlation]] = featureCalc->GetFeature(Hist2FeaturesType::Correlation);
        featurevec[IndividualFeaturesString[Energy]] = featureCalc->GetFeature(Hist2FeaturesType::Energy);
        featurevec[IndividualFeaturesString[Contrast]] = featureCalc->GetFeature(Hist2FeaturesType::Inertia);
        featurevec[IndividualFeaturesString[Entrophy]] = featureCalc->GetFeature(Hist2FeaturesType::Entropy);
        featurevec[IndividualFeaturesString[GLCMMean]] = pixelMean;
        featurevec[IndividualFeaturesString[GLCMVariance]] = pixelVariance;
        featurevec[IndividualFeaturesString[Homogenety]] = featureCalc->GetFeature(Hist2FeaturesType::InverseDifferenceMoment);
        featurevec[IndividualFeaturesString[Clustershade]] = featureCalc->GetFeature(Hist2FeaturesType::ClusterShade);
        featurevec[IndividualFeaturesString[Clusterprominence]] = featureCalc->GetFeature(Hist2FeaturesType::ClusterProminence);
        featurevec[IndividualFeaturesString[Autocorrelation]] = featureCalc->GetFeature(Hist2FeaturesType::HaralickCorrelation);
      }
    }
  }

  /**
  \brief Calculate RunLength Features
  This class computes run length descriptions from an image.By default, run length features are computed for each spatial direction and then averaged afterward,
  so it is possible to access the standard deviations of the texture features

  \param image The input Image on which to calculate the features
  \param mask The image specificying the roi
  \param feature A vector holding features of each offset direction
  \param featurevec - map of Individual feature name and thier value.

  */

  void calculateRunLength(typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec)
  {

    //typename Offsets::Pointer offsets = Offsets::New();
    // offsets = offset;

    typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
    std::vector<double> nonzero_pix;
    IteratorType inputIt(mask, mask->GetLargestPossibleRegion());

    for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    {
      if (mask->GetPixel(inputIt.GetIndex()) != 0)
      {
        nonzero_pix.push_back(image->GetPixel(inputIt.GetIndex()));
      }
    }
    int min = *min_element(nonzero_pix.begin(), nonzero_pix.end());
    int max = *max_element(nonzero_pix.begin(), nonzero_pix.end());
    typedef itk::Statistics::DenseFrequencyContainer2 HistogramFrequencyContainerType;
    typedef itk::Statistics::ScalarImageToRunLengthFeaturesFilter
      <TImageType, HistogramFrequencyContainerType> RunLengthFilterType;

    using RunLengthMatrixGenerator = typename RunLengthFilterType::RunLengthMatrixFilterType;
    using RunLengthFeatures = typename RunLengthFilterType::RunLengthFeaturesFilterType;


    //typename  RunLengthFilterType::RunLengthMatrixFilterType RunLengthMatrixGenerator;
    //typename  RunLengthFilterType::RunLengthFeaturesFilterType RunLengthFeatures;
    typename  RunLengthMatrixGenerator::Pointer matrix_generator = RunLengthMatrixGenerator::New();
    using RunLengthFeatureName = typename RunLengthFilterType::RunLengthFeaturesFilterType;
    // typename RunLengthFeatureName::RunLengthFeatureName InternalRunLengthFeatureName;
    using InternalRunLengthFeatureName = typename RunLengthFeatureName::RunLengthFeatureName;

    //  
    typename  RunLengthFilterType::FeatureNameVectorPointer requestedFeatures = RunLengthFilterType::FeatureNameVector::New();

    requestedFeatures->push_back(RunLengthFeatureName::ShortRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::LongRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::GreyLevelNonuniformity);
    requestedFeatures->push_back(RunLengthFeatureName::RunLengthNonuniformity);
    requestedFeatures->push_back(RunLengthFeatureName::LowGreyLevelRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::HighGreyLevelRunEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::ShortRunLowGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::ShortRunHighGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::LongRunLowGreyLevelEmphasis);
    requestedFeatures->push_back(RunLengthFeatureName::LongRunHighGreyLevelEmphasis);


    std::vector <std::string> featurename;
    featurename.push_back(IndividualFeaturesString[SRE]);    featurename.push_back(IndividualFeaturesString[LRE]);
    featurename.push_back(IndividualFeaturesString[GLN]);    featurename.push_back(IndividualFeaturesString[RLN]);
    featurename.push_back(IndividualFeaturesString[LGRE]);   featurename.push_back(IndividualFeaturesString[HGRE]);
    featurename.push_back(IndividualFeaturesString[SRLGE]);  featurename.push_back(IndividualFeaturesString[SRHGE]);
    featurename.push_back(IndividualFeaturesString[LRLGE]);  featurename.push_back(IndividualFeaturesString[LRHGE]);

    typename  OffsetVector::ConstIterator offsetIt;
    size_t offsetNum = 0, featureNum = 0;

    if (m_offsetSelect == "Average")
    {
      std::vector<double> tempfeatures(requestedFeatures->size(), 0);
      for (offsetIt = offset->Begin();
        offsetIt != offset->End(); offsetIt++, offsetNum++)
        matrix_generator->SetOffset(offsetIt.Value());
      matrix_generator->SetMaskImage(mask);
      matrix_generator->SetInput(image);
      matrix_generator->SetInsidePixelValue(1);
      matrix_generator->SetPixelValueMinMax(min, max);
      matrix_generator->SetDistanceValueMinMax(0, m_Range); // TOCHECK - why is this only between 0-4? P
      matrix_generator->SetNumberOfBinsPerAxis(m_Bins); // TOCHECK - needs to be statistically significant
      matrix_generator->Update();

      typename  RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
      runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
      runLengthMatrixCalculator->Update();
      typename   RunLengthFilterType::FeatureNameVector::ConstIterator fnameIt;
      //std::ostringstream ss;
      featureNum = 0;
      for (fnameIt = requestedFeatures->Begin(); fnameIt != requestedFeatures->End(); fnameIt++, featureNum++)
      {
        //ss << offsetNum;
        tempfeatures[featureNum] = tempfeatures[featureNum] + runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());
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
        matrix_generator->SetMaskImage(mask);
        matrix_generator->SetInput(image);
        matrix_generator->SetInsidePixelValue(1);
        matrix_generator->SetPixelValueMinMax(min, max);
        matrix_generator->SetDistanceValueMinMax(0, m_Range); // TOCHECK - why is this only between 0-4? SP
        matrix_generator->SetNumberOfBinsPerAxis(m_Bins); // TOCHECK - needs to be statistically significant
        matrix_generator->Update();

        typename   RunLengthFeatures::Pointer runLengthMatrixCalculator = RunLengthFeatures::New();
        runLengthMatrixCalculator->SetInput(matrix_generator->GetOutput());
        runLengthMatrixCalculator->Update();
        typename RunLengthFilterType::FeatureNameVector::ConstIterator fnameIt;
        std::ostringstream ss;
        featureNum = 0;
        for (fnameIt = requestedFeatures->Begin(); fnameIt != requestedFeatures->End(); fnameIt++, featureNum++)
        {
          ss << offsetNum;
          featurevec[featurename[featureNum] + std::to_string(offsetNum)] = runLengthMatrixCalculator->GetFeature((InternalRunLengthFeatureName)fnameIt.Value());

        }

      }
    }

  }



  /**
  \brief Calculate ShapeFeatures
  This class computes different shape features related to the roi region provided.

  \param mask The image specificying the roi
  \param featurevec - map of Individual feature name and thier value.

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
      typename ShapeLabelObjectType::Pointer labelObject = labelMap->GetNthLabelObject(0);

      //std::cout << labelObject->GetPrincipalMoments() << labelObject->GetElongation() <<
      //labelObject->GetPerimeter() << labelObject->GetRoundness() << labelObject->GetFlatness();
      double voxvol = mask->GetSpacing().GetElement(0)*mask->GetSpacing().GetElement(1)*mask->GetSpacing().GetElement(2);
      double volume = labelObject->GetNumberOfPixels() * voxvol;
      //double eccentricity = std::sqrt(1 - ((labelObject->GetPrincipalMoments()[1] * labelObject->GetPrincipalMoments()[2]) / std::pow(labelObject->GetPrincipalMoments()[0], 2)));

      auto princComps = labelObject->GetPrincipalMoments();
      double checker = princComps[1] * princComps[2] / std::pow(princComps[0], 2), eccentricity;
      if (checker > 1)
      {
        eccentricity = std::sqrt(1 + checker);
      }
      else
      {
        eccentricity = std::sqrt(1 - checker);
      }
      featurevec[IndividualFeaturesString[Eccentricity]] = eccentricity;
      featurevec[IndividualFeaturesString[Elongation]] = labelObject->GetElongation();
      featurevec[IndividualFeaturesString[Perimeter]] = labelObject->GetPerimeterOnBorder();
      featurevec[IndividualFeaturesString[Roundness]] = labelObject->GetRoundness();
      featurevec[IndividualFeaturesString[Flatness]] = labelObject->GetFlatness();
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
      auto princComps = labelObject->GetPrincipalMoments();
      double checker = std::pow(princComps[1] / princComps[0], 2), eccentricity;
      if (checker > 1)
      {
        eccentricity = std::sqrt(1 + checker);
      }
      else
      {
        eccentricity = std::sqrt(1 - checker);
      }
      featurevec[IndividualFeaturesString[Eccentricity]] = eccentricity;
      featurevec[IndividualFeaturesString[Elongation]] = labelObject->GetElongation();
      featurevec[IndividualFeaturesString[Perimeter]] = labelObject->GetPerimeterOnBorder();
      featurevec[IndividualFeaturesString[Roundness]] = labelObject->GetRoundness();
      featurevec[IndividualFeaturesString[Flatness]] = labelObject->GetFlatness();
    }
  }


  /**
  \brief Calculate VolumetricFeatures
  This class computes volumetric features related to the roi region provided.

  \param mask The image specificying the roi
  \param featurevec - map of Individual feature name and thier value.

  */
  template< class TImageTypeVolumetric = TImageType >
  void VolumetricFeatures(typename TImageTypeVolumetric::Pointer mask, std::map<std::string, double>& featurevec)
  {
    if (m_Dimension == 3)
    {
      typedef  short LabelType;
      typedef itk::Image< LabelType, TImageTypeVolumetric::ImageDimension > OutputImageType;
      typedef itk::ShapeLabelObject< LabelType, TImageTypeVolumetric::ImageDimension > ShapeLabelObjectType;
      typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

      typedef itk::ConnectedComponentImageFilter < TImageTypeVolumetric, OutputImageType > ConnectedComponentImageFilterType;
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

      auto spacing = mask->GetSpacing();
      double voxvol = spacing[0] * spacing[1] * spacing[2];
      double volume = labelObject->GetNumberOfPixels() * voxvol;

      featurevec[IndividualFeaturesString[Volume]] = volume;
      featurevec[IndividualFeaturesString[Pixels]] = labelObject->GetNumberOfPixels();
    }
    else if (m_Dimension == 2)
    {
      typedef  short LabelType;
      typedef itk::Image< LabelType, TImageTypeVolumetric::ImageDimension > OutputImageType;
      typedef itk::ShapeLabelObject< LabelType, TImageTypeVolumetric::ImageDimension> ShapeLabelObjectType;
      typedef itk::LabelMap< ShapeLabelObjectType > LabelMapType;

      typedef itk::ConnectedComponentImageFilter < TImageTypeVolumetric, OutputImageType > ConnectedComponentImageFilterType;
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
      double space = mask->GetSpacing().GetElement(0)*mask->GetSpacing().GetElement(1);
      double area = labelObject->GetNumberOfPixels() * space;
      featurevec[IndividualFeaturesString[Area]] = area;
      featurevec[IndividualFeaturesString[Pixels]] = labelObject->GetNumberOfPixels();
    }
    else
    {
      // do nothing in this case
    }
  }


  /**
  \brief Calculate HistogramFeatures
  This class computes Histogram feature  the bin frequency given a lower and upper bound and the bin number.
  \param image The input image
  \param mask The mask specificying the roi
  \param featurevec - map of Individual feature name and thier value.

  */
  void HistogramFeatures(typename TImageType::Pointer image, typename TImageType::Pointer mask, std::map< std::string, double > &featurevec)
  {

    typedef itk::MaskImageFilter< TImageType, TImageType > MaskFilterType;
    typename   MaskFilterType::Pointer maskFilter = MaskFilterType::New();
    maskFilter->SetInput(image);
    maskFilter->SetMaskImage(mask);
    maskFilter->Update();

    const unsigned int MeasurementVectorSize = 1; // Grayscale

    typedef itk::Statistics::ImageToHistogramFilter< TImageType > ImageToHistogramFilterType;

    typename  ImageToHistogramFilterType::HistogramType::MeasurementVectorType lowerBound(0);
    lowerBound.Fill(0);

    typename  ImageToHistogramFilterType::HistogramType::MeasurementVectorType upperBound(m_Bins);
    upperBound.Fill(255);

    typename ImageToHistogramFilterType::HistogramType::SizeType size(MeasurementVectorSize);
    size.Fill(m_Bins);

    typename  ImageToHistogramFilterType::Pointer imageToHistogramFilter = ImageToHistogramFilterType::New();
    imageToHistogramFilter->SetInput(maskFilter->GetOutput());
    imageToHistogramFilter->SetHistogramBinMinimum(lowerBound);
    imageToHistogramFilter->SetHistogramBinMaximum(upperBound);
    imageToHistogramFilter->SetHistogramSize(size);
    imageToHistogramFilter->Update();
    typename  ImageToHistogramFilterType::HistogramType* histogram = imageToHistogramFilter->GetOutput();

    for (unsigned int i = 0; i < histogram->GetSize()[0]; ++i)
    {
      //featurevec->emplace_back("Intensity", "Histogram_bin_frequency_" + std::to_string(i), histogram->GetFrequency(i));
      featurevec["Bin_Frequency_" + std::to_string(i)] = histogram->GetFrequency(i);
      //std::cout << histogram->GetFrequency(i) << " ";
    }

  }


  /**
  \brief Calculate GLSZM features
  A Gray Level Size Zone(GLSZM) quantifies gray level zones in an image.A gray level zone is defined as a the number of connected voxels that share the same gray level intensity.

  \param image The input image
  \param mask The mask specificying the roi
  \param featurevec - map of Individual feature name and thier value.

  */
  void CalculateGrayLevelSizeZoneFeatures(typename TImageType::Pointer itkImage, typename TImageType::Pointer maskImage, std::map<std::string, double>& featurevec)
  {
    typedef TImageType ImageType;
    //typedef TImageType MaskType;
    typedef itk::Statistics::EnhancedScalarImageToSizeZoneFeaturesFilter<ImageType> FilterType;
    typedef itk::MinimumMaximumImageCalculator<ImageType> MinMaxComputerType;
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

    typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
    std::vector<double> nonzero_pix;
    IteratorType inputIt(maskImage, maskImage->GetLargestPossibleRegion());

    for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    {
      if (inputIt.Get() != 0)
      {
        nonzero_pix.push_back(itkImage->GetPixel(inputIt.GetIndex()));
      }
    }
    int min = *min_element(nonzero_pix.begin(), nonzero_pix.end());
    int max = *max_element(nonzero_pix.begin(), nonzero_pix.end());

    typename MinMaxComputerType::Pointer minMaxComputer = MinMaxComputerType::New();
    minMaxComputer->SetImage(itkImage);
    minMaxComputer->Compute();

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
      filter->SetPixelValueMinMax(min, max);
      filter->SetNumberOfBinsPerAxis(m_Bins);
    }

    filter->SetDistanceValueMinMax(0, rangeOfPixels);

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
  \param featurevec - map of Individual feature name and thier value.

  */
  void CalculateGrayLevelNeighbourhoodGreyLevelDifferenceFeatures(typename TImageType::Pointer itkImage,
    typename TImageType::Pointer maskImage, OffsetVector *offset, std::map<std::string, double>& featurevec)
  {
    //typedef TImageType MaskType;
    typedef itk::Statistics::EnhancedScalarImageToNeighbourhoodGreyLevelDifferenceFeaturesFilter<TImageType> FilterType;
    typedef itk::MinimumMaximumImageCalculator<TImageType> MinMaxComputerType;
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


    typedef itk::ImageRegionIteratorWithIndex< TImageType > IteratorType;
    std::vector<double> nonzero_pix;
    IteratorType inputIt(maskImage, maskImage->GetLargestPossibleRegion());

    for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); ++inputIt)
    {
      if (inputIt.Get() != 0)
      {
        nonzero_pix.push_back(itkImage->GetPixel(inputIt.GetIndex()));
      }
    }
    int min = *min_element(nonzero_pix.begin(), nonzero_pix.end());
    int max = *max_element(nonzero_pix.begin(), nonzero_pix.end());

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
      filter->SetPixelValueMinMax(min, max);
      filter->SetNumberOfBinsPerAxis(m_Bins);
    }

    filter->SetDistanceValueMinMax(0, rangeOfPixels);

    filter->Update();

    auto featureMeans = filter->GetFeatureMeans();
    auto featureStd = filter->GetFeatureStandardDeviations();

    std::ostringstream  ss;
    ss << rangeOfPixels;
    std::string strRange = ss.str();
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

  // These statistics has to pulled in from cbica statistics and removed
  // TBD:: replace 'T' with 'double' for statistics and replace entire calculation with cbica::Statistics class
  T getstd(std::vector<T> input_array)
  {
    T sum = std::accumulate(input_array.begin(), input_array.end(), 0.0);
    T mean = sum / input_array.size();
    T var = 0;
    for (size_t x = 0; x < input_array.size(); x++)
    {
      var += std::pow((input_array[x] - mean), 2);
    }
    var = var / input_array.size();
    return  static_cast<double>(std::sqrt(var));
  }

  // TBD:: replace 'T' with 'double' for statistics and replace entire calculation with cbica::Statistics class
  T getvariance(std::vector<T> input_array)
  {
    T sum = std::accumulate(input_array.begin(), input_array.end(), 0.0);
    T mean = sum / input_array.size();
    T var = 0;
    for (size_t x = 0; x<input_array.size(); x++)
    {
      var += std::pow((input_array[x] - mean), 2);
    }

    return static_cast<double> (var / input_array.size());
  }

  // TBD:: replace 'T' with 'double' for statistics and replace entire calculation with cbica::Statistics class
  T getkurtosis(std::vector<T> input_array)
  {
    T sum = std::accumulate(input_array.begin(), input_array.end(), 0.0);
    T mean = sum / input_array.size();
    T var = 0;
    for (size_t x = 0; x<input_array.size(); x++)
    {
      var += std::pow((input_array[x] - mean), 2);
    }
    T stdev = std::sqrt(var);

    T kurtosis = 0.0;
    for (size_t x = 0; x<input_array.size(); x++)
    {
      kurtosis += std::pow(((input_array[x] - mean) / stdev), 4);
    }
    return static_cast<double>(kurtosis / input_array.size());

  }

  // TBD:: replace 'T' with 'double' for statistics and replace entire calculation with cbica::Statistics class
  T getskewness(std::vector<T> input_array)
  {
    T sum = std::accumulate(input_array.begin(), input_array.end(), 0.0);
    T mean = sum / input_array.size();
    T var = 0;
    for (size_t x = 0; x<input_array.size(); x++)
    {
      var += std::pow((input_array[x] - mean), 2);
    }
    T stdev = std::sqrt(var);

    T skewness = 0;
    for (size_t x = 0; x<input_array.size(); x++)
    {
      skewness += std::pow(((input_array[x] - mean) / stdev), 3);
    }
    return static_cast<double>(skewness / input_array.size());

  }



  typename TImageType::Pointer m_Mask;
  std::vector<int>m_roi;
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
  int m_Radius = 0, m_Bins = 0, m_Dimension = 0, m_Direction = 0, m_Range = 0; 
  std::string m_Axis, m_offsetSelect;

};
