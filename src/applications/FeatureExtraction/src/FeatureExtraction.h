/**
\file FeatureExtraction.h

This file holds the declaration of the class FeatureExtraction.

https://www.med.upenn.edu/sbia/software/ <br>
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
#include <cmath> 

#include "itkImage.h"
#include "itkPoint.h"
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
#include "itkImageDuplicator.h"
#include "itkExtractImageFilter.h"

#include "itkEnhancedScalarImageToNeighbourhoodGreyLevelDifferenceFeaturesFilter.h"
#include "itkEnhancedScalarImageToSizeZoneFeaturesFilter.h"
#include "ROIConstruction.h"
//#include "LBP/LBPFeatures2D.h"

#include "cbicaUtilities.h"
#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaITKUtilities.h"
#include "cbicaLogging.h"
#include "cbicaStatistics.h"

//#include "CAPTk.h"
#include "FeatureMap.h"
#include "TextureFeatureBase.h"

//! Common Definitions used in the class
using ImageType2D = itk::Image< float, 2 >;
using ImageTypeFloat3D = itk::Image< float, 3 >;

//! FeatureType is the main datatype enclosing all the information about features
using FeatureType = std::map < // MAP_1: overall control map
  std::string, // feature family (intensity or volumetric or morphologic...)
  std::tuple <
  bool, // true if it the feature(s) has/have been selected manually (via checkboxes in UI)
  std::map< // MAP_2: this map structures has all the parameters and their values parsed from param_default file for a particular feature family
  std::string, // parameter name (for instance, "Bins")
  std::map< // MAP_3: collection of properties of the specific parameter defined above
  std::string, // name of a property (for instance Type/Range/Default)
  std::string >// corresponding value of the property name (if above is Type, this is INT and so on); close of MAP_3
  >, // close of MAP_2
  std::string, // modality names, for instance T1, T2, ...
  std::string, // roi label names, for instance ED, NCR, ...
  std::map< // MAP_4 Individual features and the ouput value mapped
  std::string, // feature names, for instance Contrast, Correlation, ...
  double> // actual output value of feature specified above; close of MAP_3
  > // close of tuple
>; // close of MAP_1

using FeatureBaseType = std::vector < std::tuple< std::string, std::string, std::string, double > >;

/*
GaborWavelets,Radius,int,0:50,16,Radius of the Gabor feature calculation
GaborWavelets,FMax,float,0:10,0.25,
GaborWavelets,Gamma,float,0:10,1.414,
GaborWavelets,Directions,Int,03:13,8,The number of directions around the center voxel to calculate features on
*/
enum Params
{
  Dimension, Axis, Radius, Neighborhood, Bins, Bins_Min, Directions, Offset, Range,
  LatticeWindow, LatticeStep, LatticeBoundary, LatticePatchBoundary, LatticeWeight, 
  LatticeFullImage, GaborFMax, GaborGamma, GaborLevel, EdgesETA, EdgesEpsilon, 
  QuantizationExtent, QuantizationType, Resampling, ResamplingInterpolator_Image, 
  ResamplingInterpolator_Mask, SliceComputation, NaNHandling, 
  LBPStyle, MorphologicFeret, ParamMax
};
static const char ParamsString[ParamMax + 1][30] =
{
  "Dimension", "Axis", "Radius", "Neighborhood", "Bins", "Bins_Min", "Directions", "Offset", "Range",
  "Window", "Step", "Boundary", "PatchBoundary", "Weight", 
  "FullImage", "FMax", "Gamma", "Level", "ETA", "Epsilon", 
  "Quantization_Extent", "Quantization_Type", "Resampling", "ResamplingInterpolator_Image", 
  "ResamplingInterpolator_Mask", "SliceComputation", "NaNHandling",
  "LBPStyle", "Feret", "ParamMax"
};

enum FeatureFamily
{
  Generic,
  Intensity,
  Histogram,
  Volumetric,
  Morphologic,
  GLCM,
  GLRLM,
  GLSZM,
  NGTDM,
  NGLDM,
  LBP,
  Lattice,
  FractalDimension,
  Gabor,
  Laws,
  Edges,
  Power,
  FeatureMax
};

static const char FeatureFamilyString[FeatureMax + 1][20] =
{ "Generic", "Intensity", "Histogram", "Volumetric", "Morphologic", "GLCM", "GLRLM", "GLSZM", "NGTDM", "NGLDM", "LBP",
"Lattice", "FractalDimension", "GaborWavelets", "Laws", "EdgeEnhancement", "PowerSpectrum", "FeatureMax" };

/**
\brief FeatureExtraction Class -The main class structure enclosing all the feature calculations functions.
\ templated over Imagetype and ImageType::PixelType

*/
template< class TImageType = ImageTypeFloat3D >
class FeatureExtraction
{
public:
  // useful aliases
  using TIteratorType = itk::ImageRegionIteratorWithIndex< TImageType >;
  using TConstIteratorType = itk::ImageRegionConstIterator< TImageType >;
  using OffsetType = typename TImageType::OffsetType;
  //using Offsets = OffsetType; // TBD: replace 'Offsets' data type with 'OffsetType'
  using OffsetVector = itk::VectorContainer< unsigned char, typename TImageType::OffsetType >;
  using OffsetVectorPointer = typename OffsetVector::Pointer;

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
  void EnableDebugMode()
  {
    m_debug = true;
  }

  //! Adds extra debugging information to output message
  void EnableWritingOfIntermediateFiles()
  {
    m_writeIntermediateFiles = true;
    if (!m_outputPath.empty())
    {
      m_outputIntermediatePath = m_outputPath + "/intermediate/";
      cbica::createDir(m_outputIntermediatePath);
    }
  }

  /**
  \brief SetSelectedROIsAndLabels function populates the member varibale m_roi and m_roiLabels.

  \param roi is the mask value we are interested in calculating features from.
  \param roi_labels is the corresponding string name of the provided roi numbers.

  */
  void SetSelectedROIsAndLabels(std::string roi, std::string roi_labels);

  /**
  \brief SetSelectedROIsAndLabels function populates the member varibale m_roi and m_roiLabels.

  \param roi is the mask value we are interested in calculating features from.
  \param roi_labels is the corresponding string name of the provided roi numbers.

  */
  void SetSelectedROIsAndLabels(std::vector< std::string > roi, std::vector< std::string > roi_labels);

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
  \brief This function is used to populate the variables throughout the Update() step which are then used to write to a file

  Writes in the format: SubID,Modality,ROI,FeatureFamily,Feature,Value,Parameters

  \param modality The modality of the subject on which computation is happening
  \param label The label of the ROI on which computation is happening
  \param featureFamily Categorization of the features
  \param featureList Individual features and their values
  \param parameters The Parameters to put in the comments section
  */
  void WriteFeatures(const std::string &modality, const std::string &label, const std::string &featureFamily, const std::map< std::string, double > &featureList, const std::string &parameters,
    typename TImageType::IndexType centerWorld, bool featureMapWriteForLattice = false, float weight = 1.0);

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
  \brief SetOutputFilename sets m_outputPath value
  \param filename is the absolute path of output file.
  */
  void SetOutputFilename(std::string filename)
  {
    auto ext = cbica::getFilenameExtension(filename, false);
    if (cbica::isDir(filename))
    {
      m_outputPath = filename;
      m_outputFile = m_outputPath + "/results_" + cbica::getCurrentProcessID() + "-" + cbica::getCurrentLocalTimestamp() + ".csv";
    }
    else if (ext.empty()) // this is a directory, check if it is present and if not, create it
    {
      m_outputPath = filename;
      cbica::createDir(m_outputPath);
      m_outputFile = m_outputPath + "/results" + cbica::getCurrentProcessID() + "-" + cbica::getCurrentLocalTimestamp() + ".csv";
    }
    else
    {
      m_outputFile = filename;
      m_outputPath = cbica::getFilenamePath(m_outputFile, false);
    }

    if (!cbica::isDir(m_outputPath))
    {
      cbica::createDir(m_outputPath);
    }
    if (!cbica::isFile(m_outputFile) && !m_Features.empty())
    {
//      if (m_outputVerticallyConcatenated)
//      {
//        std::ofstream myfile;
//        myfile.open(m_outputFile, std::ofstream::out | std::ofstream::app);
//        // check for locks in a cluster environment
//        while (!myfile.is_open())
//        {
//          cbica::sleep(100);
//          myfile.open(m_outputFile, std::ios_base::out | std::ios_base::app);
//        }
//        myfile << "SubID,Modality,ROI,FeatureFamily,Feature,Value,Parameters\n";
//#ifndef WIN32
//        myfile.flush();
//#endif
//        myfile.close();
//      }
    }
    if (m_writeIntermediateFiles)
    {
      m_outputIntermediatePath = m_outputPath + "/intermediate/";
      cbica::createDir(m_outputIntermediatePath);
    }
  }

  /**
  \brief Get the complete output filename
  */
  std::string GetOutputFile() { return m_outputFile; };

  /**
  \brief Enable vertically-concatenated output in CSV
  */
  void SetVerticallyConcatenatedOutput(bool flag) 
  { 
    m_outputVerticallyConcatenated = flag; 
  }

  /**
  \brief Enable/disable feature map writing
  */
  void SetWriteFeatureMaps(bool flag) 
  { 
    m_writeFeatureMaps = flag; 
  }

  /**
  \brief Get the ROIs that were selected for the current process as a vector of ints
  */
  std::vector< int > GetSelectedROIs() { return m_roi; }

  //! Set new logger file
  void SetNewLogFile(const std::string &logFile);

  //! Set mask validation
  void SetValidMask() { m_maskValidated = true; };

  //! Set number of (OpenMP) threads for FE
  void SetNumberOfThreads(int threads) { m_threads = threads; };

  //! This is used to customize the offsets being used by the class for computation
  void SetOffsetString(std::string &offsetString) 
  {
    if (offsetString.find("|") != std::string::npos)
    {
      m_offsetString = cbica::stringSplit(offsetString, "|");
    }
    else if (offsetString.find(",") != std::string::npos)
    {
      m_offsetString = cbica::stringSplit(offsetString, ",");
    }
  }

private:

  /**
  \brief SetFeatureParam function populates the parameters for the selected feature.

  The different parameter values are set in member variables before the calculations take place

  \param selected_feature Denotes a feature name
  */
  void SetFeatureParam(std::string featureFamily);

  /**
  \brief GetRadiusInImageCoordinates function gets the radius in image coordinates

  \param radiusInWorldCoordinates Input radius in world coordinates
  \return Corresponding radius in image coordinates
  */
  int GetRadiusInImageCoordinates(float radiusInWorldCoordinates);

  /**
  \brief Calculates the OffsetVectorPointer based on the provided radius in mm and directions
  */
  std::vector< OffsetVectorPointer > GetOffsetVector(float inputRadius, int inputDirections)
  {
    auto spacing = m_Mask->GetSpacing();
    itk::Size< TImageType::ImageDimension > radius; // radius value along individual axes in image coordinates
    int actualRadius = std::numeric_limits< int >::max();
    for (size_t i = 0; i < TImageType::ImageDimension; i++)
    {
      auto temp = inputRadius / spacing[i];
      if ((temp < 1) && (temp > 0)) // this is a contingency in cases where the radius has been initialized to be less than the pixel spacing
      {
        radius[i] = 1;
      }
      else
      {
        radius[i] = std::round(temp);
      }

      auto temp_rad = static_cast<int>(radius[i]);

      if (temp_rad < actualRadius)
      {
        actualRadius = temp_rad;
      }
    }

    return GetOffsetVector(actualRadius, inputDirections);
  }

  /**
  \brief Calculates the OffsetVectorPointer based on the provided radius and directions
  */
  std::vector< OffsetVectorPointer > GetOffsetVector(int inputRadius, int inputDirections)
  {
    std::vector< OffsetVectorPointer > allOffsets;
    if (m_offsetString.empty())
    {
      if (inputRadius == -1) // this is in the contingency when someone has used mm in the param file
      {
        return GetOffsetVector(m_Radius_float, inputDirections);
      }

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

      OffsetVectorPointer offsets_total = OffsetVector::New(),
        offsets_x = OffsetVector::New(),
        offsets_y = OffsetVector::New(),
        offsets_z = OffsetVector::New();
      auto centerIndex = neighborhood.GetCenterNeighborhoodIndex();

      for (int d = directionsToCompute - 1; d >= 0; d--)
      {
        if (d != static_cast<int>(centerIndex))
        {
          offsets_total->push_back(neighborhood.GetOffset(d));
          offsets_x->push_back(neighborhood.GetOffset(d));
          offsets_y->push_back(neighborhood.GetOffset(d));
          offsets_z->push_back(neighborhood.GetOffset(d));
        }
      }
      allOffsets.push_back(offsets_total);

      if (TImageType::ImageDimension == 3)
      {
        for (size_t i = 0; i < offsets_total->size(); i++)
        {
          offsets_x->at(i)[0] = 0;
          offsets_y->at(i)[1] = 0;
          offsets_z->at(i)[2] = 0;
        }
        allOffsets.push_back(offsets_x);
        allOffsets.push_back(offsets_y);
        allOffsets.push_back(offsets_z);
      }

      return allOffsets;
    }
    else // customized offset values
    {
      itk::Size< TImageType::ImageDimension > radius; // radius value along individual axes in image coordinates
      if (inputRadius == -1) // this is in the contingency when someone has used mm in the param file
      {
        auto spacing = m_Mask->GetSpacing();
        for (size_t d = 0; d < TImageType::ImageDimension; d++)
        {
          auto temp = m_Radius_float / spacing[d];
          if ((temp < 1) && (temp > 0)) // this is a contingency in cases where the radius has been initialized to be less than the pixel spacing
          {
            radius[d] = 1;
          }
          else
          {
            radius[d] = std::round(temp);
          }
        }
      }
      else
      {
        radius.Fill(inputRadius);
      }
      OffsetVectorPointer offsets = OffsetVector::New();
      for (size_t i = 0; i < m_offsetString.size(); i++)
      {
        auto tempCurrentOffset = cbica::stringSplit(m_offsetString[i], "x");
        if (tempCurrentOffset.size() != TImageType::ImageDimension)
        {
          std::cerr << "Offset provided does not match input image/mask dimension; SubjectID: '" << m_patientID << "'\n";
          //exit(EXIT_FAILURE);
          WriteErrorFile("Offset provided does not match input image/mask dimension");
          allOffsets.push_back(offsets);
          return allOffsets;
        }
        else
        {
          OffsetType tempOffset;
          for (size_t j = 0; j < tempCurrentOffset.size(); j++)
          {
            tempOffset[j] = std::atoi(tempCurrentOffset[j].c_str()) * radius[j]; // offset has been scaled
          }
          offsets->push_back(tempOffset);
        }
      }
      allOffsets.push_back(offsets);
      return allOffsets;
    }  
  }

  /**
  \brief GetSelectedSlice enables selection of a particular slice and direction in which the features have to be generated.
  \param mask is vector of mask images.
  */
  std::vector< typename TImageType::Pointer > GetSelectedSlice(typename TImageType::Pointer mask);

  /**
  \brief GetSelectedSlice enables selection of a particular slice and direction in which the features have to be generated.
  \param mask is vector of mask images.
  \param axis represents the slice selection is in x, y or z direction.
  */
  typename TImageType::Pointer GetSelectedSlice(typename TImageType::Pointer mask, std::string axis);

  /**
  \brief Calculate intensity features
  \param image - ITKImage pointer
  \param mask - ITKImage pointer
  \param featurevec - map of Individual feature name and their value.
  */
  void CalculateIntensity(std::vector< typename TImageType::PixelType >& nonZeroVoxels, std::map< std::string, double > &featurevec, bool latticePatch = false);

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
  \param featurevec Map of Individual feature name and their value
  \param latticePatch Whether the computation is happening on a lattice patch or not
  */
  void CalculateGLCM(const typename TImageType::Pointer image, const typename TImageType::Pointer mask, OffsetVectorPointer offset, std::map<std::string, double> &featurevec, bool latticePatch = false);

  /**
  \brief Calculate RunLength Features
  This class computes run length descriptions from an image.By default, run length features are computed for each spatial direction and then averaged afterward,
  so it is possible to access the standard deviations of the texture features

  \param image The input Image on which to calculate the features
  \param mask The image specificying the roi
  \param feature A vector holding features of each offset direction
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGLRLM(const typename TImageType::Pointer image, typename TImageType::Pointer mask, OffsetVectorPointer offset, std::map<std::string, double> &featurevec, bool latticePatch = false);

  /**
  \brief Calculate CalculateMorphologic
  This class computes different shape features related to the Entropy region provided.

  \param mask The image specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  template< class TImageTypeShape = TImageType >
  void CalculateMorphologic(const typename TImageType::Pointer image, const typename TImageTypeShape::Pointer mask1, const typename TImageType::Pointer mask2, std::map<std::string, double>& featurevec);


  /**
  \brief Calculate CalculateVolumetric
  This class computes volumetric features related to the roi region provided.

  \param mask The image specificying the roi
  \param featurevec - map of Individual feature name and their value
  */
  template< class TImageTypeVolumetric = TImageType >
  void CalculateVolumetric(const typename TImageTypeVolumetric::Pointer mask, std::map<std::string, double>& featurevec);


  /**
  \brief Calculate CalculateHistogram
  This class computes Histogram feature  the bin frequency given a lower and upper bound and the bin number.

  \param image The input image
  \param mask The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateHistogram(const typename TImageType::Pointer image, const typename TImageType::Pointer mask, std::map< std::string, double > &featurevec, bool latticePatch = false);


  /**
  \brief Calculate GLSZM features
  A Gray Level Size Zone(GLSZM) quantifies gray level zones in an image.A gray level zone is defined as a the number of connected voxels that share the same gray level intensity.

  \param itkImage The input image
  \param maskImage The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGLSZM(const typename TImageType::Pointer itkImage, const typename TImageType::Pointer maskImage, OffsetVectorPointer offset, std::map<std::string, double>& featurevec);


  /**
  \brief Calculate LBP features
  
  \param itkImage The input image
  \param maskImage The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateLBP(const typename TImageType::Pointer itkImage, const typename TImageType::Pointer maskImage, std::map<std::string, double>& featurevec);

  /**
  \brief Calculate NGLDM features

  \param itkImage The input image
  \param maskImage The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateNGLDM(const typename TImageType::Pointer itkImage, const typename TImageType::Pointer maskImage, OffsetVectorPointer offset, std::map< std::string, double >& featurevec);

  /**
  \brief Calculate NGTDM features

  \param itkImage The input image
  \param maskImage The mask specifying the roi
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateNGTDM(const typename TImageType::Pointer itkImage, const typename TImageType::Pointer maskImage, OffsetVectorPointer offset, std::map< std::string, double >& featurevec);

  /**
  \brief Calculate Fractal Dimension features (box count and minkovski)

  \param itkImage The input image
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateFractalDimensions(const typename TImageType::Pointer itkImage, std::map< std::string, double > &featurevec, bool latticePatch = false);

  /**
  \brief Calculate Laws Measures Features

  \param itkImage The input image
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateLawsMeasures(const typename TImageType::Pointer itkImage, std::map< std::string, double > &featurevec);

  /**
  \brief Calculate Edge Enhancement Features

  \param itkImage The input image
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateEdgeEnhancement(const typename TImageType::Pointer itkImage, std::map< std::string, double > &featurevec);

  /**
  \brief Calculate Power Spectrum Features

  \param itkImage The input image
  \param featurevec - map of Individual feature name and their value
  */
  void CalculatePowerSpectrum(const typename TImageType::Pointer itkImage, std::map< std::string, double > &featurevec);

  /**
  \brief Calculate Gabor Wavelets features (box count and minkovski)

  \param itkImage The input image
  \param featurevec - map of Individual feature name and their value
  */
  void CalculateGaborWavelets(const typename TImageType::Pointer itkImage, std::map< std::string, double > &featurevec, bool latticePatch = false);

  /**
  \brief Get the patched image; useful only for features that do not need the overall image coordinate system defined

  Returns image and corresponding mask
  */
  typename TImageType::Pointer GetPatchedImage(const typename TImageType::Pointer inputImage);

  /**
  \brief Write the error file with some comments
  */
  void WriteErrorFile(const std::string &comments)
  {
    std::ofstream myfile(m_outputPath + "/" + m_patientID + ".error");
    myfile << comments << "\n";
    myfile.close();
  }

  //! This is used to calculate the maximum possible distance in the defined ROI - useful for GLRLM and NGLDM calculations
  float GetMaximumDistanceWithinTheDefinedROI(const typename TImageType::Pointer itkImage, const typename TImageType::Pointer maskImage);
    
  // member variables
  cbica::Statistics< typename TImageType::PixelType > m_statistics_local; //! this is for the intensity features
  typename TImageType::PixelType m_minimumToConsider, m_maximumToConsider;
  std::map< int, cbica::Statistics< typename TImageType::PixelType > > m_statistics_global; //! this is for global statistic features
  typename TImageType::Pointer m_Mask; //! the overall mask file; this is split into different ones in the GetAllROIs() function
  std::vector< int > m_roi; //! the ROI values on which calculation needs to happen
  std::vector< std::string > m_roiLabels; //! ROI label names 
  std::vector< typename TImageType::Pointer > m_inputImages; //! images on which to do computation
  std::vector< std::string > m_modality; //! size corresponds to m_inputImages.size()
  std::map< std::string, bool > m_selectedFeatures; //! selected features from GUI
  FeatureMap::Container m_featureParams; //! all selected features and their parameters
  FeatureType m_Features; //! all selected features, their parameters and outputs including whether they are enabled or not from GUI
  std::string m_outputFile; //! output file
  std::string  m_outputPath, //! this is the output directory and can be used to save intermediate files, if required
    m_outputIntermediatePath; //! store intermediate files (if any)
  bool m_outputVerticallyConcatenated = false; //! flag to check how to write the output file (whether in individual fields or vertically concatenated), defaults to horizontal-concatenation
  std::string m_finalOutputToWrite; //! this gets populated with the feature values to write at the end, TBD: needs to change when YML format is incorporated
  bool m_cancel; //! unused right now but scope for extension in the future
  bool m_debug = false; //! extra debugging information
  bool m_writeIntermediateFiles = false; //! write intermediate files in feature extraction, done in m_outputPath/intermediate
  bool m_algorithmDone = false; // check if processing is finished or not
  bool m_maskValidated = false; // whether the mask image has been validated or not
  std::string m_patientID; //! used to write first field of CSV
  std::string m_trainingFile_features; //! this is the string to populate with features for the current subject
  std::string m_trainingFile_featureNames; //! this contains all the feature names for the training file
  ROIConstruction< TImageType > m_roiConstructor; //! sets up the ROIs
  cbica::Logging m_logger;
  std::string m_separator = ","; //! the separator used during writing the output
  std::vector< std::string > m_offsetString; //! in case the user wants to customize offsets
  int m_threads = -1; //! number of OpenMP threads for FE
  int m_currentROIValue; //! the original value of the overall mask 

  std::vector< typename TImageType::PixelType > m_currentNonZeroImageValues; //! this contains the non-zero voxel values for the current ROI and always keeps changing; intended for use for different features in an easier manner
  typename TImageType::IndexType m_currentLatticeCenter, //! this is the current lattice patch center
    m_currentLatticeStart; //! this is the starting index of the current lattice patch
  bool m_LatticeComputation = false; //! flag to check if lattice-based computation has been enabled or not
  bool m_writeFeatureMaps = false; //! flag to check to write feature maps or not
  bool m_keepNaNs = true; //! whether to keep the nan values or not
  float m_latticeWindow = 0, m_latticeStep = 0; //! these are defined in mm
  typename TImageType::IndexType m_latticeStepImage; //! lattice step as defined in the image coordinates
  typename TImageType::SizeType m_latticeSizeImage; //! lattice size in image space
  bool m_fluxNeumannEnabled = false, m_zeroPaddingEnabled = false; //! the boundary conditions
  bool m_patchOnRoiEnabled = false; //! whether to pull the entire patch or only along the ROI
  bool m_patchBoundaryDisregarded = false; //! only considers patches with all pixels != 0
  bool m_patchFullImageComputation = false; //! whether computations across the entire image need to happen in addition to lattice or not
  bool m_morphologicCalculateFeret = false; //! controls calculation of feret diameter
  typename TImageType::Pointer m_featureMapBaseImage; //! the feature map base: this is only used as the base image for the lattice feature maps
  std::map< std::string, // FeatureFamily_FeatureName
    typename TImageType::Pointer > m_downscaledFeatureMaps; // each feature map (represented by the string in key)
  std::map< std::string, // FeatureFamily_FeatureName
    std::vector< double > > m_LatticeFeatures; // for lattice only to ensure consistent dimensions
  std::string m_centerIndexString; //! the center index of the current lattice in string

  // the parameters that keep changing on a per-feature basis
  int m_Radius = 0, m_Bins = 0, m_Dimension = 0, m_Direction = 0, m_neighborhood = 0, m_LBPStyle = 0;
  std::vector< int > m_Bins_range, //! range of bins to calculate features on
    m_Radius_range; //! range of radii to calculate features on
  float m_Radius_float = 0.0, m_Range = 0, 
    m_Bins_min = std::numeric_limits<float>::max(); //! the starting index of the histogram binning
  std::string m_Axis, m_offsetSelect; //! these are string based parameters
  int m_histogramBinningType = HistogramBinningType::FixedBinNumber; //! type of quantization happening, either FBN/FBS/Equal
  std::string m_QuantizationExtent = "ROI"; //! extent of quantization happening, either ROI-based or Image-based
  std::string m_initializedTimestamp; //! timestamp to append to all results - keeps outputs in sync with current process
  float m_resamplingResolution = 0.0; //! resolution to resample the images and mask to before doing any kind of computation
  std::string m_resamplingInterpolator_Image = "Linear", //! type of interpolator to use if resampling is happening, ignored if m_resamplingResolution = 0
    m_resamplingInterpolator_Mask = "Nearest";

  bool m_SliceComputation = false; //! Controls whether non-Intensity features are calculated along the slice with the largest area along the 3 axes: valid for 3D images only

  float m_gaborFMax = 0.25; //! TBD: what is the description of this?
  float m_gaborGamma = sqrtf(2); //! TBD: what is the description of this?
  int m_gaborLevel = 4; //! TBD: what is the description of this?

  float m_edgesETA = 10; //! TBD: what is the description of this?
  float m_edgesEpsilon = 10; //! TBD: what is the description of this?
};

#include "FeatureExtraction.hxx" // to process templates
