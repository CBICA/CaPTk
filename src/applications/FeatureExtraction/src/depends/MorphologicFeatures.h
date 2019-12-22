/**
\file  MorphologicFeatures.h

\brief Morphologic feature calculation

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include "itkImage.h"
#include "itkLabelGeometryImageFilter.h"

#include "itkHardConnectedComponentImageFilter.h"

#include "FeatureBase.h"

template< typename TImageType, typename TShapeImageType >
class MorphologicFeatures : public FeatureBase < TImageType >
{
public:
  //! Default constructor
  MorphologicFeatures() { };

  //! Default destructor
  ~MorphologicFeatures() { };

  void SetMaskShape(typename TShapeImageType::Pointer maskShape)
  {
    m_maskShape = maskShape;
  }

  void SetRange(int inputRange)
  {
    m_extractionType = inputRange;
  }

  void EnableCalculateFeretDiameter()
  {
    m_calculateFeretDiameter = true;
  }
  
  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      using LabelType = short;
      using OutputImageType = itk::Image< LabelType, 3 >;
      //auto connected = itk::ConnectedComponentImageFilter < TShapeImageType, OutputImageType >::New();
      auto connected = itk::HardConnectedComponentImageFilter < TShapeImageType, TShapeImageType >::New();
      connected->SetInput(m_maskShape);
      //connected->SetFullyConnected(true);
      //connected->SetBackgroundValue(0);
      connected->Update();

      //cbica::WriteImage< TShapeImageType >(connected->GetOutput(), "C:/Projects/CaPTk_myFork/src/applications/FeatureExtraction/data/ibsi_phantom/connected.nii.gz");

      /* TBD
      // using this provides distance threshold https://itk.org/Doxygen/html/classitk_1_1NeighborhoodConnectedImageFilter.html
      // this filter needs an initial seed to be placed in the mask and then the region growing to happen
      */

      auto i2l = itk::LabelImageToShapeLabelMapFilter < TShapeImageType >::New();
      i2l->SetInput(connected->GetOutput());
      if (m_calculateFeretDiameter)
      {
        i2l->ComputeFeretDiameterOn();
      }
      i2l->ComputePerimeterOn();
      //i2l->ComputeOrientedBoundingBoxOn();
      i2l->Update();
      auto labelMap = i2l->GetOutput();

      auto numberOfLabelObjects = labelMap->GetNumberOfLabelObjects();

      //if (numberOfLabelObjects > 2)
      //{
      //  std::cerr << "Number of connected components are more than 1, cannot compute morphologic features.\n";
      //  return;
      //}

      int componentToConsider = 1;
      int largestComponentSize = 0;
      if (m_extractionType == ExtractionType::Largest)
      {
        if (numberOfLabelObjects > 1)
        {
          for (size_t i = 1; i < numberOfLabelObjects; i++) // 0 is always background
          {
            auto currentComponentSize = labelMap->GetNthLabelObject(i)->GetNumberOfPixels();
            if (largestComponentSize < currentComponentSize)
            {
              componentToConsider = i;
              largestComponentSize = currentComponentSize;
            }
          }
          this->m_features["LargestComponentSize"] = largestComponentSize;
        }
        else
        {
          componentToConsider = 0;
          auto currentComponentSize = labelMap->GetNthLabelObject(componentToConsider)->GetNumberOfPixels();
          this->m_features["LargestComponentSize"] = currentComponentSize;
        }
        auto labelObject = labelMap->GetNthLabelObject(componentToConsider);

        auto ellipseDiameter = labelObject->GetEquivalentEllipsoidDiameter();
        auto orientedBoundingBoxSize = labelObject->GetOrientedBoundingBoxSize();

        double numerator = 1;
        if (ellipseDiameter.Size() == 2)
        {
          numerator = std::pow(ellipseDiameter[0], 2);
        }
        else if (ellipseDiameter.Size() == 3)
        {
          numerator = ellipseDiameter[0] * ellipseDiameter[1];
        }

        this->m_features["Eccentricity"] = std::sqrt(1 - numerator / std::pow(ellipseDiameter[ellipseDiameter.Size() - 1], 2));

        for (size_t d = 0; d < TShapeImageType::ImageDimension; d++)
        {
          this->m_features["EllipseDiameter_Axis-" + std::to_string(d)] = ellipseDiameter[d];
          this->m_features["OrientedBoundingBoxSize_Axis-" + std::to_string(d)] = orientedBoundingBoxSize[d];

          //std::cout << "EllipseDiameter" + "_Axis-" + std::to_string(d) << " = " << ellipseDiameter[d] << "\n";
        }

        if (m_calculateFeretDiameter)
        {
          this->m_features["FeretDiameter"] = labelObject->GetFeretDiameter();
        }
        this->m_features["PerimeterOnBorder"] = labelObject->GetPerimeterOnBorder();
        this->m_features["Perimeter"] = labelObject->GetPerimeter();
        this->m_features["PixelsOnBorder"] = labelObject->GetNumberOfPixelsOnBorder();
        this->m_features["PhysicalSize"] = labelObject->GetPhysicalSize();
        this->m_features["NumberOfPixels"] = labelObject->GetNumberOfPixels();
        this->m_features["EquivalentSphericalRadius"] = labelObject->GetEquivalentSphericalRadius();
        this->m_features["EquivalentSphericalPerimeter"] = labelObject->GetEquivalentSphericalPerimeter();
        this->m_features["PerimeterOnBorderRatio"] = labelObject->GetPerimeterOnBorderRatio();
        this->m_features["Roundness"] = labelObject->GetRoundness();
        this->m_features["Flatness"] = labelObject->GetFlatness();
        this->m_features["Elongation"] = labelObject->GetElongation();
        //m_features["OrientedBoundingBoxSize"] = labelObject->GetOrientedBoundingBoxSize(); // vector<double 2> cannot be converted to double

        //std::cout << "Eccentricity" << " = " << this->m_features["Eccentricity"] << "\n";
        //std::cout << "PerimeterOnBorder" << " = " << labelObject->GetPerimeterOnBorder() << "\n";
        //std::cout << "Perimeter" << " = " << labelObject->GetPerimeter() << "\n";

      }
      else
      {
        auto size = this->m_Mask->GetBufferedRegion().GetSize();
        size_t totalSizeThreshold = m_extractionType * 0.000001; // 10^-6 threshold in case all connected components are chosen
        for (size_t d = 1; d < TImageType::ImageDimension; d++)
        {
          totalSizeThreshold *= size[d];
        }

        for (size_t i = 1; i < numberOfLabelObjects; i++) // 0 is always background
        {
          auto labelObject = labelMap->GetNthLabelObject(i);

          std::string labelString = ""; // TBD: redundant
          if (numberOfLabelObjects > 2)
          {
            labelString = "-Label-" + std::to_string(i);
          }
          auto currentComponentSize = labelMap->GetNthLabelObject(i)->GetNumberOfPixels();

          if (currentComponentSize > totalSizeThreshold)
          {
            auto ellipseDiameter = labelObject->GetEquivalentEllipsoidDiameter();
            auto orientedBoundingBoxSize = labelObject->GetOrientedBoundingBoxSize();

            double numerator = 1;
            if (ellipseDiameter.Size() == 2)
            {
              numerator = std::pow(ellipseDiameter[0], 2);
            }
            else if (ellipseDiameter.Size() == 3)
            {
              numerator = ellipseDiameter[0] * ellipseDiameter[1];
            }

            this->m_features["Eccentricity" + labelString] = std::sqrt(1 - numerator / std::pow(ellipseDiameter[ellipseDiameter.Size() - 1], 2));

            for (size_t d = 0; d < TShapeImageType::ImageDimension; d++)
            {
              this->m_features["EllipseDiameter" + labelString + "_Axis-" + std::to_string(d)] = ellipseDiameter[d];
              this->m_features["OrientedBoundingBoxSize" + labelString + "_Axis-" + std::to_string(d)] = orientedBoundingBoxSize[d];

              //std::cout << "EllipseDiameter" + labelString + "_Axis-" + std::to_string(d) << " = " << ellipseDiameter[d] << "\n";
            }

            if (m_calculateFeretDiameter)
            {
              this->m_features["FeretDiameter" + labelString] = labelObject->GetFeretDiameter();
            }
            this->m_features["PerimeterOnBorder" + labelString] = labelObject->GetPerimeterOnBorder();
            this->m_features["Perimeter" + labelString] = labelObject->GetPerimeter();
            this->m_features["PixelsOnBorder" + labelString] = labelObject->GetNumberOfPixelsOnBorder();
            this->m_features["PhysicalSize" + labelString] = labelObject->GetPhysicalSize();
            this->m_features["NumberOfPixels" + labelString] = labelObject->GetNumberOfPixels();
            this->m_features["EquivalentSphericalRadius" + labelString] = labelObject->GetEquivalentSphericalRadius();
            this->m_features["EquivalentSphericalPerimeter" + labelString] = labelObject->GetEquivalentSphericalPerimeter();
            this->m_features["PerimeterOnBorderRatio" + labelString] = labelObject->GetPerimeterOnBorderRatio();
            this->m_features["Roundness" + labelString] = labelObject->GetRoundness();
            this->m_features["Flatness" + labelString] = labelObject->GetFlatness();
            this->m_features["Elongation" + labelString] = labelObject->GetElongation();
            //m_features["OrientedBoundingBoxSize"] = labelObject->GetOrientedBoundingBoxSize(); // vector<double 2> cannot be converted to double
          }
        }

      }

      this->m_algorithmDone = true;
    }
  }

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
  }

private:

  typename TShapeImageType::Pointer m_maskShape;

  enum ExtractionType
  {
    Largest, AnythingElse
  };

  int m_extractionType = ExtractionType::Largest;

  bool m_calculateFeretDiameter = false; //! substantially increases compute time
};
