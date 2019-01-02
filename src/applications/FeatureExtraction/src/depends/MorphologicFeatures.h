/**
\file  MorphologicFeatures.h

\brief Morphologic feature calculation

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#pragma once

#include "itkImage.h"
#include "itkLabelGeometryImageFilter.h"

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
  
  /**
  \brief Update calculate five feature values
  **/
  void Update()
  {
    if (!this->m_algorithmDone)
    {
      auto size = this->m_Mask->GetBufferedRegion().GetSize();
      auto minSize = size[0];
      auto maxSize = size[0];
      for (size_t d = 1; d < TImageType::ImageDimension; d++)
      {
        if (size[d] < minSize)
        {
          minSize = size[d];
        }
        if (size[d] > maxSize)
        {
          maxSize = size[d];
        }
      }

      using LabelType = short;
      using OutputImageType = itk::Image< LabelType, TShapeImageType::ImageDimension >;
      auto connected = itk::ConnectedComponentImageFilter < TShapeImageType, OutputImageType >::New();
      connected->SetInput(m_maskShape);
      connected->SetFullyConnected(true);
      connected->SetBackgroundValue(itk::NumericTraits< LabelType >(0));
      connected->Update();

      /* TBD
      // using this provides distance threshold https://itk.org/Doxygen/html/classitk_1_1NeighborhoodConnectedImageFilter.html
      // this filter needs an initial seed to be placed in the mask and then the region growing to happen
      */

      auto i2l = itk::LabelImageToShapeLabelMapFilter < OutputImageType, itk::LabelMap< itk::ShapeLabelObject< LabelType, TShapeImageType::ImageDimension > > >::New();
      i2l->SetInput(connected->GetOutput());
      //i2l->ComputeFeretDiameterOn();
      i2l->ComputePerimeterOn();
      //i2l->ComputeOrientedBoundingBoxOn();
      i2l->Update();
      auto labelMap = i2l->GetOutput();

      auto numberOfLabelObjects = labelMap->GetNumberOfLabelObjects();

      if (numberOfLabelObjects > 1)
      {
        std::cerr << "Number of connected components are more than 1, cannot compute morphologic features.\n";
        return;
      }

      for (unsigned int i = 0; i < numberOfLabelObjects; i++)
      {
        auto labelObject = labelMap->GetNthLabelObject(i);

        std::string labelString = "";
        if (numberOfLabelObjects > 1)
        {
          labelString = "-Label-" + std::to_string(i);
        }
        /// TBD: convert this into a function of minSize or maxSize instead of the hard-coded '20'
        auto tt = labelObject->GetNumberOfPixels();
        if (labelObject->GetNumberOfPixels() > (0.1 * minSize)) // this is to ensure really small portions of an ROI do not get picked 
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

          this->m_features["Eccentricity" + labelString] = sqrt(1 - numerator / std::pow(ellipseDiameter[ellipseDiameter.Size() - 1], 2));

          for (size_t d = 0; d < TShapeImageType::ImageDimension; d++)
          {
            this->m_features["EllipseDiameter" + labelString + "_Axis-" + std::to_string(d)] = ellipseDiameter[d];
            this->m_features["OrientedBoundingBoxSize" + labelString + "_Axis-" + std::to_string(d)] = orientedBoundingBoxSize[d];
          }
          //m_features["FeretDiameter" + labelString] = labelObject->GetFeretDiameter();
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
        else
        {
          std::cerr << "LabelObject '" << i << "' is too small (total pixels = " << labelObject->GetNumberOfPixels() << ").\n";
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
};
