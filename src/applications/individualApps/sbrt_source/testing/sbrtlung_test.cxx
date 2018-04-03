#include <algorithm>
#include <stdlib.h>
#include <iostream>
#include <string>

#include "itkImageFileReader.h"
#include "cbicaUtilities.h"
#include "cbicaLogging.h"
#include "cbicaCmdParser.h"
#include "cbicaITKSafeImageIO.h"
#include "cbicaPreProcessImage.h"
#include "Joint_Segm.h"
#include "itkTestingComparisonImageFilter.h"

//const unsigned int  Dimension = 3;
typedef  float           PixelType;
typedef itk::Image< PixelType, Dimension > ImageType;
typedef ImageType::Pointer InputImagePointerType;

int main(int argc, char** argv)
{
  std::string datadirectory, ctImageFileName, petImageFileName, maskImageFileName;

  Joint_Segm segmenation;

  bool differenceFailed = false;
  //Joint_Segm segmenation;


  datadirectory = "../../data";

  // ensure dataDir ends with a '/'
  if (datadirectory[datadirectory.length() - 1] != '/')
  {
    datadirectory.append("/");
  }

  if (!cbica::isDir(datadirectory))
  {

    return EXIT_FAILURE;
  }

  std::vector< std::string > detectedSubjects = cbica::subdirectoriesInDirectory(datadirectory);
  for (size_t subject = 0; subject < detectedSubjects.size(); subject++)
  {
    const std::string subjectUnderConsideration = datadirectory + detectedSubjects[subject] + "/";
    std::vector< std::string > filesPerSubject = cbica::filesInDirectory(subjectUnderConsideration);
    InputImagePointerType imgin = ImageType::New();
    InputImagePointerType PET_imgin = ImageType::New();


    for (size_t j = 0; j < filesPerSubject.size(); j++)
    {
      FileNameParts  fileUnderConsideration = FileNameParts(subjectUnderConsideration + filesPerSubject[j]);
      if (!cbica::fileExists(fileUnderConsideration.fullFileName))
      {
        //  std::cout << "Aborting:   Could not read " << fileUnderConsideration.fullFileName << std::endl;
        return EXIT_FAILURE;
      }
      if (fileUnderConsideration.base == "crop_Ct2Pet")
      {
        if (!cbica::fileExists(fileUnderConsideration.fullFileName))
        {
          // std::cout << "Aborting:   Could not read " << fileUnderConsideration.fullFileName << std::endl;
          return EXIT_FAILURE;
        }
        else
          imgin = cbica::ReadImage(fileUnderConsideration.fullFileName);
      }
      if (fileUnderConsideration.base == "crop_Pet")
      {
        if (!cbica::fileExists(fileUnderConsideration.fullFileName))
        {
          // std::cout << "Aborting:   Could not read " << fileUnderConsideration.fullFileName << std::endl;
          return EXIT_FAILURE;
        }
        PET_imgin = cbica::ReadImage(fileUnderConsideration.fullFileName);
      }

    }
    // std::cout << " Input Images read.." << std::endl;
    maskImageFileName = subjectUnderConsideration + "mask.nii";
    auto sbrtseeds = segmenation.generateseeds<itk::Image< float, 3> >(imgin, PET_imgin);


    InputImagePointerType groundtruthseeds = cbica::ReadImage(subjectUnderConsideration + "gt.nii");
    InputImagePointerType test_result = ImageType::New();


    typedef itk::Testing::ComparisonImageFilter<ImageType, ImageType> DiffType;
    DiffType::Pointer diff = DiffType::New();
    diff->SetTestInput(sbrtseeds);
    diff->SetValidInput(groundtruthseeds);
    diff->UpdateLargestPossibleRegion();

    double numberOfPixelsTolerance = 10;

    const double averageIntensityDifference = diff->GetTotalDifference();
    const unsigned long numberOfPixelsWithDifferences =
      diff->GetNumberOfPixelsWithDifferences();
    if (averageIntensityDifference > 0.0)
    {
      if (static_cast<int>(numberOfPixelsWithDifferences) >
        numberOfPixelsTolerance)
      {
        return EXIT_FAILURE;
      }
    }
  }

  return EXIT_SUCCESS;
}