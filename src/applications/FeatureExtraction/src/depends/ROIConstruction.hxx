#pragma once

template< class TImage >
void ROIConstruction< TImage >::SetInputMask(typename TImage::Pointer input)
{
  m_mask = input;
  m_algorithmDone = false;
}

template< class TImage >
void ROIConstruction< TImage >::SetNewLogFile(const std::string &logFile)
{
  m_logger.UseNewFile(logFile);
}

template< class TImage >
void ROIConstruction< TImage >::SetSelectedROIsAndLabels(std::vector< int > roi, std::vector< std::string > roi_labels)
{
  m_roi = roi;
  m_roiLabels = roi_labels;
  m_algorithmDone = false;
}

template< class TImage >
void ROIConstruction< TImage >::ROIComputation(std::vector< typename TImage::IndexType > &inputLatticeGrid,
  itk::Size< TImage::ImageDimension > &inputLatticeRadius)
{
  TConstIteratorType totalMaskIterator(m_mask, m_mask->GetLargestPossibleRegion());
  totalMaskIterator.GoToBegin();

  if (!m_roi.empty())
  {
    m_output.resize(m_roi.size());

    while (!totalMaskIterator.IsAtEnd())
    {
      auto currentTotalMaskImageValue = totalMaskIterator.Get();
      if (currentTotalMaskImageValue > 0)
      {
        auto currentIndex = totalMaskIterator.GetIndex();

        for (size_t i = 0; i < m_roi.size(); i++)
        {
          if (currentTotalMaskImageValue == static_cast< typename TImage::PixelType >(m_roi[i]))
          {
            //if (m_output[i].maskImage.IsNull() || (m_output[i].maskImage->GetLargestPossibleRegion().GetSize()[0] == 0)) // check if this image has been initialized before
            //{
            //  m_output[i].maskImage = cbica::CreateImage< TImage >(m_mask); // create a blank image taking the properties from total_mask for the current ROI
            //}
            //TIteratorType maskIterator(m_output[i].maskImage, m_output[i].maskImage->GetLargestPossibleRegion());

            //maskIterator.SetIndex(currentIndex);
            //maskIterator.Set(1);

            m_output[i].value = m_roi[i]; // set the value of the current ROI 
            m_output[i].label = m_roiLabels[i];
            m_output[i].nonZeroIndeces.push_back(currentIndex); // store the non-zero pixel/voxel indeces of the image for the current ROI
            typename TImage::IndexType currentIndex;
            currentIndex.Fill(0);
            m_output[i].centerIndex = currentIndex;
          }
        } // end of m_roi iteration loop
      } // end of currentTotalMaskImageValue check
      ++totalMaskIterator;
    }// end of totalMaskIterator while loop

     // construct the lattice patches
    if (!inputLatticeGrid.empty()) // append "_Full" to the ROI name 
    {
      // neighborhood iterator is using ConstantBoundaryCondition, which defaults to using '0' at the boundaries
      using BoundaryCondition = itk::ConstantBoundaryCondition< TImage >;
      //if (m_fluxNeumannEnabled) // TBD: no idea if this actually works
      //{
      //  using BoundaryCondition = itk::ZeroFluxNeumannBoundaryCondition< TImage >;
      //}
      //else if (m_zeroPaddingEnabled)
      //{
      //  using BoundaryCondition = itk::ConstantBoundaryCondition< TImage >;
      //}
      using TNeighborhoodIteratorType = itk::NeighborhoodIterator< TImage, BoundaryCondition >;
      TNeighborhoodIteratorType maskNeighborhoodIt(inputLatticeRadius, m_mask, m_mask->GetBufferedRegion());

      for (size_t j = 0; j < m_roiLabels.size(); j++)
      {
        m_output[j].label = m_roiLabels[j] + "_Full";
      }

      for (size_t i = 0; i < inputLatticeGrid.size(); ++i)
      {
        typename TImage::IndexType currentIndex = inputLatticeGrid[i];
        //m_mask->TransformPhysicalPointToIndex(inputLatticeGrid[i], currentIndex);

        maskNeighborhoodIt.SetLocation(currentIndex);
        auto currentValue = maskNeighborhoodIt.GetCenterValue();
        ROIProperties currentGridNodeProperties;

        // calculate the neighborhood only if the current value is part of the roi calculation 
        auto roiItr = std::find(m_roi.begin(), m_roi.end(), *currentValue);
        // change the condition of this if statement to consider all grid nodes regardless of their center value
        if (roiItr != m_roi.end()) // the value is present in the requested ROIs
        {
          m_latticePoints++;
          auto currentGridIndexString = "[" + std::to_string(inputLatticeGrid[i][0]); // this will be used to distinguish the "patch" ROI
          for (size_t d = 1; d < TImage::ImageDimension; d++)
          {
            currentGridIndexString += ";" + std::to_string(inputLatticeGrid[i][d]);
          }
          currentGridIndexString += "]";

          size_t roiLocation = roiItr - m_roi.begin();
          currentGridNodeProperties.label = m_roiLabels[roiLocation] + "_Lattice" /*+ currentGridIndexString*/; // append current grid node to the original roi label 
          currentGridNodeProperties.value = m_roi[roiLocation];
          currentGridNodeProperties.centerIndex = currentIndex;
          currentGridNodeProperties.latticeGridPoint = true;
          
          size_t roiValCounter = 0; // used to save the number of voxels having the same value as center

          for (size_t n = 0; n < maskNeighborhoodIt.Size(); n++)
          {
            if (maskNeighborhoodIt.IndexInBounds(n))
            {
              auto currentNeighborhoodItIndex = maskNeighborhoodIt.GetIndex(n);
              auto currentNeighborhoodItValue = static_cast< int >(maskNeighborhoodIt.GetPixel(n));

              if (m_patchOnRoiEnabled)// only consider the values that are the same as the current ROI (i.e., the center index) value
              {
                if (currentNeighborhoodItValue == m_roi[roiLocation]) // when the center index is of the selected ROI and everything non-zero is considered
                {
                  currentGridNodeProperties.nonZeroIndeces.push_back(currentNeighborhoodItIndex); // store the non-zero pixel/voxel indeces of the image for the current ROI
                  roiValCounter++;
                }
                if (m_patchBoundaryDisregarded) // in the case of "None", as defined by CBIG texture feature pipeline
                {
                  m_maskBackgroundWeighting = true;
                }
                else
                {
                  m_maskBackgroundWeighting = false;
                }
              }
              else // consider everything that is not part of the background
              {
                // increment the counter when the value is same as the center index
                if (currentNeighborhoodItValue == m_roi[roiLocation])
                {
                  roiValCounter++;
                }
                currentGridNodeProperties.nonZeroIndeces.push_back(currentNeighborhoodItIndex); // store the non-zero pixel/voxel indeces of the image for te mask
                m_maskBackgroundWeighting = true;
              }
            }
          } // maskNeighborhoodIt loop ends

          // check for the background weighting
          if (m_maskBackgroundWeighting)
          {
            // the border portion is never used until the border calculation has been enabled
            //bool tempBorderPixelDetected = false; // checks if the currentIndex is at the border of the image or not
            //auto currentSize = m_mask->GetBufferedRegion().GetSize();
            //size_t numberOfBorders = 0;
            //for (size_t d = 0; d < TImage::ImageDimension; d++)
            //{
            //  if (currentIndex[i] == 0) // at the origin in any of the axes
            //  {
            //    tempBorderPixelDetected = true; // this gets hit even at places not at the border, for some reason
            //    numberOfBorders++;
            //  }
            //  else if (currentIndex[i] == currentSize[i] - 1) // at the end in any of the axes
            //  {
            //    tempBorderPixelDetected = true;
            //    numberOfBorders++;
            //  }
            //}
            //if (tempBorderPixelDetected) // if border pixel has been detected, divide size by the number of detected border regions
            //{
            //  currentGridNodeProperties.weight = static_cast< float >(roiValCounter) / static_cast< float >(maskNeighborhoodIt.Size() / numberOfBorders);
            //}
            //else // do the actual calculation if currentIndex is not at the border
            {
              currentGridNodeProperties.weight = static_cast< float >(roiValCounter) / static_cast< float >(maskNeighborhoodIt.Size());
            }
          }
          else // otherwise, the weight is 1
          {
            currentGridNodeProperties.weight = 1.0;
          }
          m_output.push_back(currentGridNodeProperties);
        } // m_roi value check
      } // inputLatticeGrid iterator
    } // inputLatticeGrid size check
  } // end of m_roi.empty() check
  else
  {
    m_logger.WriteError("m_roi is empty in ROIConstruction::NonLatticeComputation");
  }
}

template< class TImage >
void ROIConstruction< TImage >::Update()
{
  if (!m_algorithmDone) // if this flag has not been set, it means processing has not happened
  {
    auto t1 = std::chrono::high_resolution_clock::now();
    if (m_mask.IsNull())
    {
      m_logger.WriteError("The mask image is not defined or has not been read properly.");
      exit(EXIT_FAILURE);
    }
    // compute the physical coordinate where the image ends 
    TConstIteratorType totalMaskIterator(m_mask, m_mask->GetLargestPossibleRegion());
    auto size = m_mask->GetLargestPossibleRegion().GetSize();

    totalMaskIterator.GoToBegin();
    auto start_imageCoordinates = totalMaskIterator.GetIndex(), end_imageCoordinates = start_imageCoordinates;
    for (size_t i = 0; i < TImage::ImageDimension; i++)
    {
      end_imageCoordinates[i] = size[i] - 1;
    }

    //! use the type identified for one of the coordiantes to declare the lattice - fancy C++11 stuff
    std::vector< decltype(start_imageCoordinates) > latticeGrid;

    auto distances = cbica::GetDistances< TImage >(m_mask);

    // in case the LatticeWindow has not been defined, don't bother checking for LatticeStep do lattice computation
    if (m_latticeWindow != 0)
    {
      m_latticeEnabled = true;
      auto spacing = m_mask->GetSpacing();
      for (size_t i = 0; i < TImage::ImageDimension; i++)
      {
        if (m_latticeWindow >= distances[i])
        {
          m_logger.WriteError("LatticeWindow is greater than or equal to the image size in world coordinates along the axis '" + std::to_string(i) + "', so lattice has been disabled");
          m_latticeEnabled = false;
        }

        // check for latticeStep
        if (m_latticeStep >= distances[i])
        {
          m_logger.WriteError("LatticeStep is greater than or equal to the image size in world coordinates along the axis '" + std::to_string(i) + "', so lattice has been disabled");
          m_latticeEnabled = false;
        }

        // calculate the latticeWindow, i.e., the radius for neighborhoodIterator, in the image space
        auto temp = m_latticeWindow / spacing[i];
        if ((temp < 1) && (temp > 0)) // this is a contingency in cases where the window has been initialized to be less than the pixel spacing
        {
          m_latticeRadius[i] = 1;
        }
        else
        {
          m_latticeRadius[i] = std::round(temp / 2);
        }

        if (m_latticeRadius[i] == 0)
        {
          m_logger.WriteError("LatticeWindow in the image space along the axis '" + std::to_string(i) + "' has been computed to '0', so lattice has been disabled");
          m_latticeEnabled = false;
        }
      }
    }
    else if (m_latticeStep > 0)
    {
      m_logger.WriteError("LatticeStep has been initialized but not LatticeWindow, so lattice has been disabled");
    }

    // if both the latticeWindow and latticeStep are valid sizes, proceed with grid construction
    if (m_latticeEnabled)
    {
      m_logger.Write("Found valid latticeStep and latticeWindowSize; lattice-based computation has been enabled");
      //totalROIsToCompute *= numberOfSteps;

      ///
      // change this entire section to construct the lattice using image coordinates because of NIfTI orientation issues
      ///

      auto latticeHalfStepImage = m_latticeStepImage;
      for (size_t d = 0; d < TImage::ImageDimension; d++)
      {
        auto temp = std::round(m_latticeStepImage[d] / 2);
        if (temp > 0)
        {
          auto tempDiv = m_latticeRadius[d] / temp;
          if ((tempDiv > 1) || std::isnan(tempDiv) || std::isinf(tempDiv))
          {
            latticeHalfStepImage[d] = std::ceil(tempDiv) * temp;
          }
          else
          {
            latticeHalfStepImage[d] = temp;
          }
        }
        else
        {
          latticeHalfStepImage[d] = 1;
        }
      }

      switch (TImage::ImageDimension)
      {
      case 2:
      {
        for (auto x = latticeHalfStepImage[0]; x < end_imageCoordinates[0] - latticeHalfStepImage[0]; x += m_latticeStepImage[0])
        {
          for (auto y = latticeHalfStepImage[1]; y < end_imageCoordinates[1] - latticeHalfStepImage[1]; y += m_latticeStepImage[1])
          {
            auto currentIndex = start_imageCoordinates;
            currentIndex[0] = x;
            currentIndex[1] = y;

            if (m_mask->GetPixel(currentIndex) != 0)
            {
              latticeGrid.push_back(currentIndex);
            }
          } // end of y
        } // end of x
        break;
      }
      case 3:
      {
        for (auto x = latticeHalfStepImage[0]; x < end_imageCoordinates[0] - latticeHalfStepImage[0]; x += m_latticeStepImage[0])
        {
          for (auto y = latticeHalfStepImage[1]; y < end_imageCoordinates[1] - latticeHalfStepImage[1]; y += m_latticeStepImage[1])
          {
            for (auto z = latticeHalfStepImage[2]; z < end_imageCoordinates[2] - latticeHalfStepImage[2]; z += m_latticeStepImage[2])
            {
              auto currentIndex = start_imageCoordinates;
              currentIndex[0] = x;
              currentIndex[1] = y;
              currentIndex[2] = z;

              if (m_mask->GetPixel(currentIndex) != 0)
              {
                latticeGrid.push_back(currentIndex);
              }
            }
          } // end of y
        } // end of x

        break;
      }
      default:
        break;
      }
    }

    totalMaskIterator.GoToBegin();

    if (m_roi.empty())
    {
      // No ROIs have been specified, so calculated for ALL of them.
      // this loop populates m_roi and m_roiLabels so that when the Update() function is called again,
      // the "else" does not get triggered

      std::cout << "Either the selected ROI value or Label is empty, using whatever non-zero value(s) is/are present in the mask for computation.\n";
      auto tempValuesInMask = cbica::GetUniqueValuesInImage< TImage >(m_mask, true);

      // this next step is to ensure that when there are very small differences in the floating point values,
      // the int cast takes care of it and we only want to keep the unique elements
      auto tempROIs = m_roi;
      for (size_t i = 0; i < tempValuesInMask.size(); i++)
      {
        tempROIs.push_back(static_cast<int>(tempValuesInMask[i]));
      }
      tempROIs.erase(
        std::unique(tempROIs.begin(), tempROIs.end()), // only use the unique elements
        tempROIs.end());

      for (size_t i = 0; i < tempROIs.size(); i++)
      {
        if (tempROIs[i] != 0)
        {
          m_roi.push_back(tempROIs[i]);
          m_roiLabels.push_back(std::to_string(tempROIs[i]));
        }
      }
      // at this point, m_roi and m_roiLabels have been populated
    } // end of m_roi.empty() check

    ROIComputation(latticeGrid, m_latticeRadius);

    m_algorithmDone = true;
    auto t2 = std::chrono::high_resolution_clock::now();
    // std::cout << "ROI Construction took " << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << " milliseconds\n";
  } // m_algorithmDone check
}
