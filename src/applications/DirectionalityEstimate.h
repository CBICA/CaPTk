#pragma once

#include "cbicaUtilities.h"
#include "cbicaITKUtilities.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaCmdParser.h"
#include "cbicaITKImageInfo.h"

#include "itkConnectedThresholdImageFilter.h"


/**
\brief This class computes the maximum distance and corresponding coordinate from a seed point in a label map

Requirements: label map or mask and a coordinate defined inside the label map

*/
//using TImageType = itk::Image< float, 3 >;
template< class TImageType = itk::Image< float, 3 > >
class DirectionalityEstimate
{
public:
  //! Default Constructor
  DirectionalityEstimate() {};

  //! Default Destructor
  ~DirectionalityEstimate() {};

  //! Set input Mask
  void SetInputMask(typename TImageType::Pointer inputMask)
  {
    m_mask = inputMask;
    algorithmDone = false;
  }

  /**
  \brief Set the input seed

  This is *ALWAYS* in image coordinates and not in world coordinates
  */
  void SetInputSeed(typename TImageType::IndexType inputSeed)
  {
    m_seed = inputSeed;
    algorithmDone = false;
  }

  //! Run the algorithm
  void Update()
  {
    if (!algorithmDone)
    {
      auto imageSize = m_mask->GetLargestPossibleRegion().GetSize();
      for (size_t i = 0; i < TImageType::ImageDimension; i++)
      {
        if (static_cast<float>(m_seed[i]) >= static_cast<float>(imageSize[i]))
        {
          std::cerr << "Seed for Directionality placed out of bounds.\n";
          return;
        }
      }

      // setup the connected component segmentation
      auto connectedComponentFilter = itk::ConnectedThresholdImageFilter< TImageType, TImageType >::New();
      connectedComponentFilter->SetInput(m_mask);
      connectedComponentFilter->SetSeed(m_seed);
      connectedComponentFilter->SetReplaceValue(1);
      // we only want the selected voxel value to be segmented and nothing else
      auto currentPixelVal = m_mask->GetPixel(m_seed);
      connectedComponentFilter->SetLower(currentPixelVal);
      connectedComponentFilter->SetUpper(currentPixelVal);
      connectedComponentFilter->Update();

      auto inputLabelMask = connectedComponentFilter->GetOutput();

      // m_outputMask contains the quadrant-based result 
      //m_outputMask = typename TImageType::New();
      //m_outputMask->SetRegions(inputLabelMask->GetLargestPossibleRegion());
      //m_outputMask->SetOrigin(inputLabelMask->GetOrigin());
      //m_outputMask->SetDirection(inputLabelMask->GetDirection());
      //m_outputMask->SetSpacing(inputLabelMask->GetSpacing());
      //m_outputMask->Allocate();
      //m_outputMask->FillBuffer(0);

      itk::ImageRegionConstIterator <TImageType> iterator(inputLabelMask, inputLabelMask->GetLargestPossibleRegion());
      
      // non-const iterator because we need write access here
      //itk::ImageRegionIteratorWithIndex< TImageType > iterator_outputMask(m_outputMask, m_outputMask->GetLargestPossibleRegion());

      // iterate through the whole image and find maximum distance
      for (iterator.GoToBegin(); !iterator.IsAtEnd(); ++iterator)
      {
        if (iterator.Get() > 0)
        {
          auto currentIndex = iterator.GetIndex();
          //iterator_outputMask.SetIndex(currentIndex);
          //iterator_outputMask.Set(1);
          float currentDist = 0.0;// TBD  gcc is unable to deduce suitable type. Please fix this -> GetDistanceBetweenIndeces(currentIndex, indexToUse);
          currentDist = cbica::GetDistanceBetweenIndeces<TImageType>(currentIndex, m_seed);

          if (currentDist > m_maxDistance)
          {
            m_maxDistance = currentDist;
            m_outputIndex = currentIndex;
            //iterator_outputMask.Set(2);
          }
        }
      }

      algorithmDone = true;
    }
    
  }

  //! Get output index
  typename TImageType::IndexType GetOutputIndex()
  {
    if (!algorithmDone)
    {
      Update();
    }
    return m_outputIndex;
  }

  //! Get output maximum distance
  float GetOutputMaxDistance()
  {
    if (!algorithmDone)
    {
      Update();
    }
    return m_maxDistance;
  }

  /** 
  \brief Get the output indeces and distances in the 3 planes and the overall 

  Returns a map with corresponding index-distance pairs for 'ALL', 'XY', 'YZ' and 'XZ'
  */
  std::map< std::string, std::pair< typename TImageType::IndexType, float > > GetOutputFull()
  {
    auto output_index = GetOutputIndex();

    auto output_xy_index = output_index, output_yz_index = output_index, output_xz_index = output_index;

    output_xy_index[2] = m_seed[2];
    output_yz_index[0] = m_seed[0];
    output_xz_index[1] = m_seed[1];

    auto dist_00_index = cbica::GetDistanceBetweenIndeces(m_seed, output_index);
    auto dist_xy_index = cbica::GetDistanceBetweenIndeces(m_seed, output_xy_index);
    auto dist_yz_index = cbica::GetDistanceBetweenIndeces(m_seed, output_yz_index);
    auto dist_xz_index = cbica::GetDistanceBetweenIndeces(m_seed, output_xz_index);

    std::map< std::string, std::pair< typename TImageType::IndexType, float > > output;
    output["ALL"] = std::make_pair(output_index, cbica::GetDistanceBetweenIndeces(m_seed, output_index));
    output["XY"] = std::make_pair(output_xy_index, cbica::GetDistanceBetweenIndeces(m_seed, output_xy_index));
    output["YZ"] = std::make_pair(output_yz_index, cbica::GetDistanceBetweenIndeces(m_seed, output_yz_index));
    output["XZ"] = std::make_pair(output_xz_index, cbica::GetDistanceBetweenIndeces(m_seed, output_xz_index));

    return output;
  }

  ////! Get the output label mask with quadrant based result
  //typename TImageType::Pointer GetOutputLabelMask()
  //{
  //  if (!algorithmDone)
  //  {
  //    Update();
  //  }
  //  return m_outputMask;
  //}

private:
  typename TImageType::Pointer m_mask;

  //typename TImageType::Pointer m_outputMask;

  typename TImageType::IndexType m_seed;

  typename TImageType::IndexType m_outputIndex;

  float m_maxDistance = 0;

  bool algorithmDone = false;
};

//
//int main(int argc, char **argv)
//{
//  auto parser = cbica::CmdParser(argc, argv);
//  parser.addRequiredParameter("l", "labelMap", cbica::Parameter::FILE, "NIfTI file", "The label map on which calculations are to be done");
//  parser.addRequiredParameter("o", "outputFile", cbica::Parameter::FILE, "text file", "Where output gets stored");
//  parser.addOptionalParameter("r", "real", cbica::Parameter::FLOAT, "labelMap Dimensions", "The real world coordinates (in mm) of file",
//    "Take example from tissue point file", "Needs to be in same dimensionality as labelMap", "Delineation is done using ','");
//  parser.addOptionalParameter("i", "index", cbica::Parameter::INTEGER, "labelMap Dimensions", "The image indeces (in voxels) of file",
//    "Needs to be in same dimensionality as labelMap", "Delineation is done using ','");
//
//  std::string file_labelMap, file_output;
//
//  parser.getParameterValue("l", file_labelMap);
//  parser.getParameterValue("o", file_output);
//
//  std::string output_message = "Output Point:\nCoordinate = ";
//  std::string seedPoint;
//  bool realCoordinatesPassed = false;
//
//  if ((argc < 2) || (argc > 5))
//  {
//    parser.echoUsage();
//    return EXIT_FAILURE;
//  }
//  if (parser.isPresent("u"))
//  {
//    parser.echoUsage();
//    return EXIT_SUCCESS;
//  }
//
//  if (parser.isPresent("h"))
//  {
//    parser.echoHelp();
//    return EXIT_SUCCESS;
//  }
//
//  if (parser.isPresent("v"))
//  {
//    parser.echoVersion();
//    return EXIT_SUCCESS;
//  }
//
//  if (parser.isPresent("r"))
//  {
//    realCoordinatesPassed = true;
//    parser.getParameterValue("r", seedPoint);
//  }
//
//  if (!realCoordinatesPassed)
//  {
//    if (parser.isPresent("i"))
//    {
//      parser.getParameterValue("i", seedPoint);
//    }
//  }
//
//  auto seedPoint_vec = cbica::stringSplit(seedPoint, ",");
//
//  auto imageInfo = cbica::ImageInfo(file_labelMap);
//
//  switch (imageInfo.GetImageDimensions())
//  {
//  case 2:
//  {
//    using ImageType = itk::Image< float, 2 >;
//
//    auto inputLabel = cbica::ReadImage< ImageType >(file_labelMap);
//
//    ImageType::IndexType seedPoints_index, output_index;
//
//    auto imageOrigin = inputLabel->GetOrigin();
//    auto imageSpacing = inputLabel->GetSpacing();
//
//    if (realCoordinatesPassed)
//    {
//      seedPoints_index[0] = std::abs((std::atof(seedPoint_vec[0].c_str()) - imageOrigin[0]) / imageSpacing[0]);
//      seedPoints_index[1] = std::abs((std::atof(seedPoint_vec[1].c_str()) - imageOrigin[1]) / imageSpacing[1]);
//    }
//    else
//    {
//      seedPoints_index[0] = std::atof(seedPoint_vec[0].c_str());
//      seedPoints_index[1] = std::atof(seedPoint_vec[1].c_str());
//    }
//
//    auto output = cbica::GetMaxDistanceInLabelMap< ImageType >(inputLabel, seedPoints_index);
//    output_index = output.second;
//
//    output_message += "[" + std::to_string(output_index[0]) + "," + std::to_string(output_index[1]) + "]"
//      + "; Distance = " + std::to_string(output.first) + "\n";
//
//    break;
//  }
//  case 3:
//  {
//    using ImageType = itk::Image< float, 3 >;
//
//    auto inputLabel = cbica::ReadImage< ImageType >(file_labelMap);
//
//    ImageType::IndexType seedPoints_index, output_index;
//
//    auto imageOrigin = inputLabel->GetOrigin();
//    auto imageSpacing = inputLabel->GetSpacing();
//
//    if (realCoordinatesPassed)
//    {
//      seedPoints_index[0] = std::abs((std::atof(seedPoint_vec[0].c_str()) - imageOrigin[0]) / imageSpacing[0]);
//      seedPoints_index[1] = std::abs((std::atof(seedPoint_vec[1].c_str()) - imageOrigin[1]) / imageSpacing[1]);
//      seedPoints_index[2] = std::abs((std::atof(seedPoint_vec[2].c_str()) - imageOrigin[2]) / imageSpacing[2]);
//    }
//    else
//    {
//      seedPoints_index[0] = std::atof(seedPoint_vec[0].c_str());
//      seedPoints_index[1] = std::atof(seedPoint_vec[1].c_str());
//      seedPoints_index[2] = std::atof(seedPoint_vec[2].c_str());
//    }
//
//    auto output = cbica::GetMaxDistanceInLabelMap< ImageType >(inputLabel, seedPoints_index);
//    output_index = output.second;
//
//    output_message += "[" + std::to_string(output_index[0]) + "," + std::to_string(output_index[1]) + std::to_string(output_index[2]) + "]"
//      + "; Distance = " + std::to_string(output.first) + "\n";
//
//    break;
//  }
//  default:
//    break;
//  }
//
//  std::ofstream myFile(file_output);
//  myFile << output_message;
//  myFile.close();
//
//  std::cout << "Finished successfully\n";
//
//  return EXIT_SUCCESS;
//}