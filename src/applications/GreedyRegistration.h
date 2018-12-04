/**
\file FeatureExtraction.h

This file holds the declaration of the class FeatureExtraction.

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
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

#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"
#include "cbicaStatistics.h"

//! Common Definitions used in the class
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
  Dimension, Axis, Radius, Neighborhood, Bins, Directions, Offset, Range, ParamMax
};
static const char ParamsString[ParamMax + 1][15] =
{
  "Dimension", "Axis", "Radius", "Neighborhood", "Bins", "Directions", "Offset", "Range", "ParamMax"
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
  LBP,
  Lattice,
  //Gabor, 
  //Laws,
  //FractalDimension,   // to be added when testing completes
  //PowerSpectrun
  FeatureMax
};

static const char FeatureFamilyString[FeatureMax + 1][15] = { "Intensity", "Histogram", "Volumetric", "Morphologic", "GLCM", "GLRLM", "GLSZM", "NGTDM", "LBP", "Lattice", "FeatureMax" };

/**
\brief FeatureExtraction Class -The main class structure enclosing all the feature calculations functions.
\ templated over Imagetype and ImageType::PixelType

*/
template< class TImageType = ImageTypeFloat3D >
class FeatureExtraction
{
  /**
  \brief Helper structure to store the ROI properties for Feature Extraction
  */
  struct ROIProperties
  {
    int value; // the value of the ROI coming from the complete mask image
    typename TImageType::Pointer maskImage; // the mask image which is always initiated as a binary image
    std::vector< typename TImageType::PixelType > nonZeroPixelValues; // for statistical computations
    std::vector< typename TImageType::IndexType > m_nonZeroIndeces;
  };

public:
  // useful aliases
  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;
  using OffsetType = typename TImageType::OffsetType;
  //using Offsets = OffsetType; // TBD: replace 'Offsets' data type with 'OffsetType'
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
  void Update();

  //! Adds extra debugging information to output message
  void SetDebugMode()
  {
    m_debug = true;
  }

  /**
  \brief SetSelectedLabels function populates the member varibale m_roi and m_labelname.

  \param roi is the mask value we are interested in calculating features from.
  \param roi_labels is the corresponding string name of the provided roi numbers.

  */
  void SetSelectedLabels(std::string roi, std::string roi_labels);

  /**
  \brief SetSelectedLabels function populates the member varibale m_roi and m_labelname.

  \param roi is the mask value we are interested in calculating features from.
  \param roi_labels is the corresponding string name of the provided roi numbers.

  */
  void SetSelectedLabels(std::vector< std::string > roi, std::vector< std::string > roi_labels);

  /**
  \brief SetRequestedFeatures function populates the main data structure m_features with the requested features and updates the bool value related to the features.

  \param filename is the default param file name provided .
  \param selected_features is the command line option given to user to pass or select required features from available ones.

  */
  void SetRequestedFeatures(std::string filename, std::string selected_features = "");

  /**
  \brief SetRequestedFeatures function populates the main data structure m_features with the requested features and updates the bool value related to the features.

  This overloaded function is for the use of UI

  \param filename is the default param file name provided
  \param selected_features is the map of feature name as key and bool as the value
  */
  void SetRequestedFeatures(std::string filename, std::map<std::string, bool> selected_features);

  /**
  \brief Populates the feature structure from the GUI
  */
  void SetRequestedFeatures(std::map< std::string, std::vector< std::map<std::string, std::string> > >  featuresFromUI, std::map<std::string, bool> selected_features);

  /**
  \brief WriteFeatureList function is to write the calculated feature values to the ouput file given.
  \param filename is the absolute path of the csv file name
  \param featurevec is a vector of base feature type , enclosing the output to each file
  \param outputfiletype is option to switch between .xml and .csv
  */
  void WriteFeatureList(std::string Filename, std::vector< FeatureType > featurevec, std::string outputfiletype = "");

  //! this is to be used in case this class becomes long-running
  void cancel()
  {
    m_cancel = true;
  }

  //! Does as expected
  void SetPatientID(const std::string &patientID)
  {
    m_patientID = patientID;
  }

  /**
  \brief SetInputImages parses the input images name and populates the member variable m_inputImages
  \param images is vector of input images.
  \param modality is a string containing the respective image modalities separated by |
  */
  void SetInputImages(std::vector< typename TImageType::Pointer > images, std::string modality);


  /**
  \brief SetInputImages parses the input images name and populates the member variable m_inputImages
  \param images is vector of input images.
  \param modality is a vector of string containing the respective image modalities.
  */
  void SetInputImages(std::vector< typename TImageType::Pointer > images, std::vector< std::string >& modality);

  /**
  \brief SetMaskImage parses the input mask and populates the member variable m_Mask
  \param image is Image pointer.
  */
  void SetMaskImage(typename TImageType::Pointer image)
  {
    m_Mask = image;
  }

  /**
  \brief SetOutputFilename sets m_ouputpath value
  \param filename is the absolute path of output file.
  */
  void SetOutputFilename(std::string filename)
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
  \brief SetFeatureParam function populates the parameters for the selected feature.

  The different parameter values are set in member variables before the calculations take place

  \param selected_feature Denotes a feature name
  */
  void SetFeatureParam(std::string selected_feature);

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
    auto centerIndex = neighborhood.GetCenterNeighborhoodIndex();

    for (int d = 0; d < directionsToCompute; d++)
    {
      if (d != static_cast<int>(centerIndex))
      {
        offsets->push_back(neighborhood.GetOffset(d));
      }
    }

    return offsets;
  }

  /**
  \brief GetSelectedSlice enables selection of a particular slice and direction in which the features have to be generated.
  \param mask is vector of mask images.
  \param axis represents the slice selection is in x, y or z direction.
  */
  ImageType2D::Pointer GetSelectedSlice(typename TImageType::Pointer mask, std::string axis);

  /**
  \brief Calculates individual mask images for each ROI provided, if m_roi is empty, mask images are estimated for every ROI present
  */
  std::vector< ROIProperties > GetAllROIs(const typename TImageType::Pointer inputImage, const typename TImageType::Pointer total_mask);

  /**
  \brief Calculate intensity features
  \param image - ITKImage pointer
  \param mask - ITKImage pointer
  \param featurevec - map of Individual feature name and their value.
  */
  void CalculateIntensity(std::vector< typename TImageType::PixelType >& nonZeroVoxels, std::map< std::string, double > &featurevec);

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
  void CalculateGLCM(typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec);

  /**
  \brief Calculate RunLength Features
  This class computes run length descriptions from an image.By default, run length features are computed for each spatial direction and then averaged afterward,
  so it is possible to access the standard deviations of the texture features

  \param image The input Image on which to calculate the features
  \param mask The image specificying the roi
  \param feature A vector holding features of each offset direction
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGLRLM(typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVector *offset, std::map<std::string, double> &featurevec);

  /**
  \brief Calculate CalculateShape
  This class computes different shape features related to the Entropy region provided.

  \param mask The image specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  template< class TImageTypeShape = TImageType >
  void CalculateShape(typename TImageTypeShape::Pointer mask, std::map<std::string, double>& featurevec);


  /**
  \brief Calculate CalculateVolumetric
  This class computes volumetric features related to the roi region provided.

  \param mask The image specificying the roi
  \param featurevec - map of Individual feature name and their value
  */
  template< class TImageTypeVolumetric = TImageType >
  void CalculateVolumetric(typename TImageTypeVolumetric::Pointer mask, std::map<std::string, double>& featurevec);


  /**
  \brief Calculate CalculateHistogram
  This class computes Histogram feature  the bin frequency given a lower and upper bound and the bin number.
  \param image The input image
  \param mask The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateHistogram(typename TImageType::Pointer image, typename TImageType::Pointer mask, std::map< std::string, double > &featurevec);


  /**
  \brief Calculate GLSZM features
  A Gray Level Size Zone(GLSZM) quantifies gray level zones in an image.A gray level zone is defined as a the number of connected voxels that share the same gray level intensity.

  \param image The input image
  \param mask The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGLSZM(typename TImageType::Pointer itkImage, typename TImageType::Pointer maskImage, std::map<std::string, double>& featurevec);


  /**
  \brief Calculate NGTDM features

  \param image The input image
  \param mask The mask specificying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateNGTDM(typename TImageType::Pointer itkImage, typename TImageType::Pointer maskImage, OffsetVector *offset, std::map< std::string, double >& featurevec);

  cbica::Statistics< typename TImageType::PixelType > m_statistics;
  typename TImageType::Pointer m_Mask; //! the overall mask file; this is split into different ones in the GetAllROIs() function
  std::vector< int > m_roi; //! the ROI values on which calculation needs to happen
  std::vector< std::string > m_labelname; //! ROI label names 
  std::vector< typename TImageType::Pointer > m_inputImages; //! images on which to do computation
  std::vector< std::string > m_modality; //! size corresponds to m_inputImages.size()
  std::map< std::string, bool > m_selectedFeatures; //! selected features from GUI
  //FeatureMap::Container m_featureParams; //! all selected features and their parameters
  FeatureType m_Features; //! all selected features, their parameters and outputs including whether they are enabled or not from GUI
  std::vector< FeatureType> m_outputFeatureVector; //! the final output for all requested features for all ROIs -- may become redundant after writing is moved to a per-feature basis TBD
  std::string m_File; //! output file
  std::string  m_ouputpath; //! this is the output directory and can be used to save intermediate files, if required
  bool m_cancel; //! unused right now but scope for extension in the future
  bool m_debug = false; //! extra debugging information
  std::string m_patientID; //! used to write first field of CSV

  // the parameters that keep changing on a per-feature basis
  int m_Radius = 0, m_Bins = 0, m_Dimension = 0, m_Direction = 0, m_Range = 0, m_neighborhood = 0;
  std::string m_Axis, m_offsetSelect;
};

