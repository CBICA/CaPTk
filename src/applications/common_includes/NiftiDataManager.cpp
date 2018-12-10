#include "NiftiDataManager.h"
//#include "CAPTk.h"
#include "CaPTkEnums.h"

NiftiDataManager::NiftiDataManager()
{
}

NiftiDataManager::~NiftiDataManager()
{
}
ImageTypeFloat3D::Pointer NiftiDataManager::ReadNiftiImage(std::string filename)
{
  ImageTypeFloat3D::Pointer image;
  itk::ImageIOBase::Pointer reader = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
  if (!reader)
  {
    /*	mLastError = "Unable to read file.";
    return false;*/
  }

  reader->SetFileName(filename);
  reader->ReadImageInformation();
  std::string InputPixelType = reader->GetComponentTypeAsString(reader->GetComponentType());

  image = ReadImageWithDimAndInputPixelType<ImageTypeFloat3D::PixelType, ImageTypeFloat3D::PixelType, ImageTypeFloat3D::ImageDimension>(filename);
  return image;
}

itk::Image<float, 4>::Pointer NiftiDataManager::Read4DNiftiImage(std::string filename)
{
  itk::Image<float, 4>::Pointer image;
  itk::ImageIOBase::Pointer reader = itk::ImageIOFactory::CreateImageIO(filename.c_str(), itk::ImageIOFactory::ReadMode);
  if (!reader)
  {
    /*	mLastError = "Unable to read file.";
    return false;*/
  }
  reader->SetFileName(filename);
  reader->ReadImageInformation();
  std::string InputPixelType = reader->GetComponentTypeAsString(reader->GetComponentType());
  image = ReadImageWithDimAndInputPixelType<ImageTypeFloat3D::PixelType, ImageTypeFloat3D::PixelType, 4>(filename);
  return image;
}

void NiftiDataManager::LoadTrainingData(ImageTypeFloat3D::Pointer labelImagePointer,
  ImageTypeFloat3D::Pointer recurrenceMaskImage,
  ImageTypeFloat3D::Pointer nonrecurrenceMaskImage,
  ImageTypeFloat3D::Pointer t1ceImagePointer,
  ImageTypeFloat3D::Pointer t2flairImagePointer,
  ImageTypeFloat3D::Pointer t1ImagePointer,
  ImageTypeFloat3D::Pointer t2ImagePointer,
  ImageTypeFloat4D::Pointer perfImagePointerNifti,
  ImageTypeFloat3D::Pointer axImagePointer,
  ImageTypeFloat3D::Pointer faImagePointer,
  ImageTypeFloat3D::Pointer radImagePointer,
  ImageTypeFloat3D::Pointer trImagePointer,
  VectorVectorDouble &pNearIntensities,
  VectorVectorDouble &pFarIntensities,
  VectorVectorDouble &mNearIntensities,
  VectorVectorDouble &mFarIntensities,
  VectorDouble &tNearIntensities,
  VectorDouble &tFarIntensities,
  int imagetype,
  bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData)
{
  typedef ImageTypeFloat3D ImageType;
  typedef ImageTypeFloat4D PerfusionImageType;
  std::vector<ImageType::IndexType> recurIndices;
  std::vector<ImageType::IndexType> nonrecurIndices;
  std::vector<ImageType::IndexType> tumorIndices;

  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType RecMaskIt(recurrenceMaskImage, recurrenceMaskImage->GetLargestPossibleRegion());
  IteratorType NonRecMaskIt(nonrecurrenceMaskImage, nonrecurrenceMaskImage->GetLargestPossibleRegion());

  RecMaskIt.GoToBegin();
  NonRecMaskIt.GoToBegin();

  while (!RecMaskIt.IsAtEnd())
  {
    if (RecMaskIt.Get() != 0)
      recurIndices.push_back(RecMaskIt.GetIndex());
    if (NonRecMaskIt.Get() != 0)
      nonrecurIndices.push_back(NonRecMaskIt.GetIndex());

    ++RecMaskIt;
    ++NonRecMaskIt;
  }
  int recurSize = recurIndices.size();
  int nonrecurSize = nonrecurIndices.size();
  //--------------------perfusion intensities---------------------
  if (usePerfData)
  {
    auto timeStamps = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3];

    for (int i = 0; i < recurSize; i++)
    {
      VectorDouble perfusionIntensitiesPerVoxel;
      PerfusionImageType::IndexType perfVoxelIndex;

      perfVoxelIndex[0] = recurIndices[i][0];
      perfVoxelIndex[1] = recurIndices[i][1];
      perfVoxelIndex[2] = recurIndices[i][2];

      for (unsigned int j = 0; j < timeStamps; j++)
      {
        perfVoxelIndex[3] = j;
        perfusionIntensitiesPerVoxel.push_back(std::round(perfImagePointerNifti.GetPointer()->GetPixel(perfVoxelIndex)));
      }
      pNearIntensities.push_back(perfusionIntensitiesPerVoxel);
    }
    for (int i = 0; i < nonrecurSize; i++)
    {
      VectorDouble perfusionIntensitiesPerVoxel;
      PerfusionImageType::IndexType perfVoxelIndex;

      perfVoxelIndex[0] = nonrecurIndices[i][0];
      perfVoxelIndex[1] = nonrecurIndices[i][1];
      perfVoxelIndex[2] = nonrecurIndices[i][2];
      for (unsigned int j = 0; j < timeStamps; j++)
      {
        perfVoxelIndex[3] = j;
        perfusionIntensitiesPerVoxel.push_back(std::round(perfImagePointerNifti.GetPointer()->GetPixel(perfVoxelIndex)));
      }
      pFarIntensities.push_back(perfusionIntensitiesPerVoxel);
    }
  }
  //--------------------remaining intensities---------------------
  VectorVectorDouble mIntensitiesNearPoints;
  VectorVectorDouble mIntensitiesFarPoints;
  for (int i = 0; i < recurSize; i++)
  {
    VectorDouble otherIntensitiesPerVoxel;
	if (useConventionalData)
	{
		otherIntensitiesPerVoxel.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(recurIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(recurIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(recurIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(recurIndices[i])));
	}
    if (useDTIData)
    {
      otherIntensitiesPerVoxel.push_back(std::round(faImagePointer.GetPointer()->GetPixel(recurIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(radImagePointer.GetPointer()->GetPixel(recurIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(trImagePointer.GetPointer()->GetPixel(recurIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(axImagePointer.GetPointer()->GetPixel(recurIndices[i])));
    }
    mNearIntensities.push_back(otherIntensitiesPerVoxel);
  }
  for (int i = 0; i < nonrecurSize; i++)
  {
    VectorDouble otherIntensitiesPerVoxel;
	if (useConventionalData)
	{
		otherIntensitiesPerVoxel.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
	}
    if (useDTIData)
    {
      otherIntensitiesPerVoxel.push_back(std::round(faImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(radImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(trImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(axImagePointer.GetPointer()->GetPixel(nonrecurIndices[i])));
    }
    mFarIntensities.push_back(otherIntensitiesPerVoxel);
  }

  //---------------------------distance metric------------------------------------
  if (useDistData)
  {
	  ImageType::Pointer DistMap = GetDistanceMap<ImageType>(labelImagePointer);
	  for (int i = 0; i < recurSize; i++)
		  tNearIntensities.push_back(DistMap.GetPointer()->GetPixel(recurIndices[i]));
	  for (int i = 0; i < nonrecurSize; i++)
		  tFarIntensities.push_back(DistMap.GetPointer()->GetPixel(nonrecurIndices[i]));
  }
}


std::vector<ImageTypeFloat3D::IndexType> NiftiDataManager::LoadTestData(ImageTypeFloat3D::Pointer t1ceImagePointer,
  ImageTypeFloat3D::Pointer t2flairImagePointer,
  ImageTypeFloat3D::Pointer t1ImagePointer,
  ImageTypeFloat3D::Pointer t2ImagePointer,
  ImageTypeFloat4D::Pointer perfImagePointerNifti,
  ImageTypeFloat3D::Pointer axImagePointer,
  ImageTypeFloat3D::Pointer faImagePointer,
  ImageTypeFloat3D::Pointer radImagePointer,
  ImageTypeFloat3D::Pointer trImagePointer,
  ImageTypeFloat3D::Pointer labelImagePointer,
  ImageTypeFloat3D::Pointer dilatedEdemaPointer,
  VectorVectorDouble & pIntensities,
  VectorVectorDouble & mIntensities,
  VectorDouble & tIntensities,
  int imagetype,
  bool useConventionalData,bool useDTIData, bool usePerfData, bool useDistData)
{
  typedef ImageTypeFloat3D ImageType;
  typedef ImageTypeFloat4D PerfusionImageType;
  //-----------------------get edema and tumor pixels-----------------------
  std::vector<ImageType::IndexType> testIndices;
  std::vector<ImageType::IndexType> tumorIndices;
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  //IteratorType EdeIt(dilatedEdemaPointer, dilatedEdemaPointer->GetLargestPossibleRegion());

  //EdeIt.GoToBegin();

  //while (!EdeIt.IsAtEnd())
  //{
	 // if (EdeIt.Get() == GLISTR_OUTPUT_LABELS::EDEMA_1 || EdeIt.Get() == GLISTR_OUTPUT_LABELS::EDEMA_2)
  //    testIndices.push_back(EdeIt.GetIndex());

  //  ++EdeIt;
  //}
  if (useDistData)
  {
    IteratorType TumIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
    TumIt.GoToBegin();
	while (!TumIt.IsAtEnd())
	{
		if (TumIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::TUMOR || TumIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::NONENHANCING)
			tumorIndices.push_back(TumIt.GetIndex());
		else if (TumIt.Get() == CAPTK::GLISTR_OUTPUT_LABELS::EDEMA)
			testIndices.push_back(TumIt.GetIndex());
		++TumIt;
	}
  }
  //--------------------perfusion intensities---------------------
  if (usePerfData)
  {
    PerfusionImageType::RegionType regionperf = perfImagePointerNifti->GetLargestPossibleRegion();
    auto timeStamps = perfImagePointerNifti->GetLargestPossibleRegion().GetSize()[3];
    for (unsigned int i = 0; i < testIndices.size(); i++)
    {
      VectorDouble perfusionIntensitiesPerVoxel;
      PerfusionImageType::IndexType perfVoxelIndex;

      perfVoxelIndex[0] = testIndices[i][0];
      perfVoxelIndex[1] = testIndices[i][1];
      perfVoxelIndex[2] = testIndices[i][2];

      for (unsigned int j = 0; j < timeStamps; j++)
      {
        perfVoxelIndex[3] = j;
		perfusionIntensitiesPerVoxel.push_back(std::round(perfImagePointerNifti.GetPointer()->GetPixel(perfVoxelIndex)));
      }
      pIntensities.push_back(perfusionIntensitiesPerVoxel);
    }
  }
  //--------------------remaining intensities---------------------
  for (unsigned int i = 0; i < testIndices.size(); i++)
  {
    VectorDouble otherIntensitiesPerVoxel;
	if (useConventionalData)
	{
		otherIntensitiesPerVoxel.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(testIndices[i])));
	}
    if (useDTIData)
    {
		otherIntensitiesPerVoxel.push_back(std::round(faImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(radImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(trImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(axImagePointer.GetPointer()->GetPixel(testIndices[i])));
    }
    if (otherIntensitiesPerVoxel.size()>0)
      mIntensities.push_back(otherIntensitiesPerVoxel);
  }
  //---------------------------distance metric------------------------------------
  if (useDistData)
  {
    for (unsigned int i = 0; i < testIndices.size(); i++)
    {
      ImageType::IndexType testIndex = testIndices[i];
      double min = 0;

      for (unsigned int j = 0; j < tumorIndices.size(); j++)
      {
        ImageType::IndexType tumorIndex = tumorIndices[j];
        double distance = std::sqrt(std::pow(testIndex[0] - tumorIndex[0], 2) + std::pow(testIndex[1] - tumorIndex[1], 2) + std::pow(testIndex[2] - tumorIndex[2], 2));
        if (j == 0)
          min = distance;
        else
        {
          if (distance < min)
            min = distance;
        }
      }
      tIntensities.push_back(min);
    }
  }
  return testIndices;
}

void NiftiDataManager::LoadTrainingData(ImageTypeFloat3D::Pointer labelImagePointer,
	VectorVectorDouble nearPoints,
	VectorVectorDouble farPoints,
	ImageTypeFloat3D::Pointer t1ceImagePointer,
	ImageTypeFloat3D::Pointer t2flairImagePointer,
	ImageTypeFloat3D::Pointer t1ImagePointer,
	ImageTypeFloat3D::Pointer t2ImagePointer,
	std::vector<ImageTypeFloat3D::Pointer> perfImagePointer,
	std::vector<ImageTypeFloat3D::Pointer> dtiImagePointer,
	VectorVectorDouble & pNearIntensities,
	VectorVectorDouble & pFarIntensities,
	VectorVectorDouble & mNearIntensities,
	VectorVectorDouble & mFarIntensities,
	VectorDouble & tNearIntensities,
	VectorDouble & tFarIntensities,
	int imagetype,
	bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData)
{
	typedef ImageTypeFloat3D ImageType;
	std::vector<ImageType::IndexType> recurIndices;
	std::vector<ImageType::IndexType> nonrecurIndices;
	//std::vector<ImageType::IndexType> tumorIndices;
	//typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
	//if (useDistData)
	//{
	//	IteratorType TumIt(labelImagePointer, labelImagePointer->GetLargestPossibleRegion());
	//	TumIt.GoToBegin();
	//	while (!TumIt.IsAtEnd())
	//	{
	//		if (TumIt.Get() == 175 || TumIt.Get() == 200)
	//			tumorIndices.push_back(TumIt.GetIndex());
	//		++TumIt;
	//	}
	//}
	std::vector<ImageType::IndexType> nearIndices;
	std::vector<ImageType::IndexType> farIndices;
	for (unsigned int i = 0; i < nearPoints.size(); i++)
	{
		ImageType::IndexType voxelIndex;
		voxelIndex[0] = nearPoints[i][0];
		voxelIndex[1] = nearPoints[i][1];
		voxelIndex[2] = nearPoints[i][2];
		nearIndices.push_back(voxelIndex);
	}
	for (unsigned int i = 0; i < farPoints.size(); i++)
	{
		ImageType::IndexType voxelIndex;
		voxelIndex[0] = farPoints[i][0];
		voxelIndex[1] = farPoints[i][1];
		voxelIndex[2] = farPoints[i][2];
		farIndices.push_back(voxelIndex);
	}
	//--------------------perfusion intensities---------------------
	if (usePerfData)
	{
		const unsigned int timeStamps = perfImagePointer.size();
		for (unsigned int i = 0; i < nearIndices.size(); i++)
		{
			VectorDouble perfusionIntensitiesPerVoxel;
			for (unsigned int j = 0; j < timeStamps; j++)
				perfusionIntensitiesPerVoxel.push_back(std::round(perfImagePointer[j].GetPointer()->GetPixel(nearIndices[i])));
			pNearIntensities.push_back(perfusionIntensitiesPerVoxel);
		}
		for (unsigned int i = 0; i < farIndices.size(); i++)
		{
			VectorDouble perfusionIntensitiesPerVoxel;
			for (unsigned int j = 0; j < timeStamps; j++)
				perfusionIntensitiesPerVoxel.push_back(std::round(perfImagePointer[j].GetPointer()->GetPixel(farIndices[i])));
			pFarIntensities.push_back(perfusionIntensitiesPerVoxel);
		}

	}
	//--------------------remaining intensities---------------------
	VectorVectorDouble mIntensitiesNearPoints,
		mIntensitiesFarPoints;
	for (unsigned int i = 0; i < nearIndices.size(); i++)
	{
		VectorDouble otherIntensitiesPerVoxel;

		if (useConventionalData)
		{
			otherIntensitiesPerVoxel.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(nearIndices[i])));
				otherIntensitiesPerVoxel.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(nearIndices[i])));
				otherIntensitiesPerVoxel.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(nearIndices[i])));
				otherIntensitiesPerVoxel.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(nearIndices[i])));
		}
		if (useDTIData)
		{
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[0].GetPointer()->GetPixel(nearIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[1].GetPointer()->GetPixel(nearIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[2].GetPointer()->GetPixel(nearIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[3].GetPointer()->GetPixel(nearIndices[i])));
		}
		mNearIntensities.push_back(otherIntensitiesPerVoxel);
	}
	for (unsigned int i = 0; i < farIndices.size(); i++)
	{
		VectorDouble otherIntensitiesPerVoxel;
		if (useConventionalData)
		{
			otherIntensitiesPerVoxel.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(farIndices[i])));
				otherIntensitiesPerVoxel.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(farIndices[i])));
				otherIntensitiesPerVoxel.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(farIndices[i])));
				otherIntensitiesPerVoxel.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(farIndices[i])));
		}
		if (useDTIData)
		{
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[0].GetPointer()->GetPixel(farIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[1].GetPointer()->GetPixel(farIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[2].GetPointer()->GetPixel(farIndices[i])));
			otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[3].GetPointer()->GetPixel(farIndices[i])));
		}
		mFarIntensities.push_back(otherIntensitiesPerVoxel);
	}

	//---------------------------distance metric------------------------------------
	if (useDistData)
	{
		ImageType::Pointer DistMap = GetDistanceMap<ImageType>(labelImagePointer);
		for (unsigned int i = 0; i < nearIndices.size(); i++)
			tNearIntensities.push_back(DistMap.GetPointer()->GetPixel(nearIndices[i]));
		for (unsigned int i = 0; i < farIndices.size(); i++)
			tFarIntensities.push_back(DistMap.GetPointer()->GetPixel(farIndices[i]));


		//  VariableLengthVectorType distances;
		//  distances.SetSize(tumorIndices.size());

		//  for (unsigned int j = 0; j < distances.Size(); j++)
		//    distances[j] = std::sqrt(std::pow(nearIndices[i][0] - tumorIndices[j][0], 2) + std::pow(nearIndices[i][1] - tumorIndices[j][1], 2) + std::pow(nearIndices[i][2] - tumorIndices[j][2], 2));

		//  double min = distances[0];
		//  for (unsigned int j = 1; j < distances.Size(); j++)
		//  {
		//    if (distances[j] < min)
		//      min = distances[j];
		//  }
		//  tNearIntensities.push_back(min);
		//}
		//for (unsigned int i = 0; i < farIndices.size(); i++)
		//{
		//  VariableLengthVectorType distances;
		//  distances.SetSize(tumorIndices.size());

		//  for (unsigned int j = 0; j < distances.Size(); j++)
		//    distances[j] = std::sqrt(std::pow(farIndices[i][0] - tumorIndices[j][0], 2) + std::pow(farIndices[i][1] - tumorIndices[j][1], 2) + std::pow(farIndices[i][2] - tumorIndices[j][2], 2));

		//  double min = distances[0];
		//  for (unsigned int j = 1; j < distances.Size(); j++)
		//  {
		//    if (distances[j] < min)
		//      min = distances[j];
		//  }
		//  tFarIntensities.push_back(min);
	}
}

std::vector<ImageTypeFloat3D::IndexType> NiftiDataManager::LoadTestData(ImageTypeFloat3D::Pointer t1ceImagePointer,
  ImageTypeFloat3D::Pointer t2flairImagePointer,
  ImageTypeFloat3D::Pointer t1ImagePointer,
  ImageTypeFloat3D::Pointer t2ImagePointer,
  std::vector<ImageTypeFloat3D::Pointer> perfImagePointer,
  std::vector<ImageTypeFloat3D::Pointer> dtiImagePointer,
  ImageTypeFloat3D::Pointer labelImagePointer,
  ImageTypeFloat3D::Pointer dilatedEdemaPointer,
  VectorVectorDouble & pIntensities,
  VectorVectorDouble & mIntensities,
  VectorDouble & tIntensities,
  int imagetype,
  bool useConventionalData, bool useDTIData, bool usePerfData, bool useDistData)
{
  typedef ImageTypeFloat3D ImageType;
  ////-----------------------get edema and tumor pixels-----------------------
  std::vector<ImageType::IndexType> testIndices;
  typedef itk::ImageRegionIteratorWithIndex <ImageType> IteratorType;
  IteratorType EdeIt(dilatedEdemaPointer, dilatedEdemaPointer->GetLargestPossibleRegion());
  EdeIt.GoToBegin();
  while (!EdeIt.IsAtEnd())
  {
    if (EdeIt.Get() == CAPTK::VOXEL_STATUS::ON)
      testIndices.push_back(EdeIt.GetIndex());
    ++EdeIt;
  }
  //--------------------perfusion intensities---------------------
  if (usePerfData)
  {
    for (unsigned int i = 0; i < testIndices.size(); i++)
    {
      VectorDouble perfusionIntensitiesPerVoxel;
      for (unsigned int j = 0; j < perfImagePointer.size(); j++)
        perfusionIntensitiesPerVoxel.push_back(perfImagePointer[j].GetPointer()->GetPixel(testIndices[i]));
      pIntensities.push_back(perfusionIntensitiesPerVoxel);
    }
  }
  //--------------------remaining intensities---------------------
  for (unsigned int i = 0; i < testIndices.size(); i++)
  {
    VectorDouble otherIntensitiesPerVoxel;
	if (useConventionalData)
	{
		otherIntensitiesPerVoxel.push_back(std::round(t1ImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t1ceImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t2ImagePointer.GetPointer()->GetPixel(testIndices[i])));
		otherIntensitiesPerVoxel.push_back(std::round(t2flairImagePointer.GetPointer()->GetPixel(testIndices[i])));
	}
    if (useDTIData)
    {
      otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[0].GetPointer()->GetPixel(testIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[1].GetPointer()->GetPixel(testIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[2].GetPointer()->GetPixel(testIndices[i])));
      otherIntensitiesPerVoxel.push_back(std::round(dtiImagePointer[3].GetPointer()->GetPixel(testIndices[i])));
    }
    if (otherIntensitiesPerVoxel.size() > 0)
      mIntensities.push_back(otherIntensitiesPerVoxel);
  }
  //---------------------------distance metric------------------------------------
  if (useDistData)
  {
	  ImageType::Pointer DistMap = GetDistanceMap<ImageType>(labelImagePointer);
	  for (unsigned int i = 0; i < testIndices.size(); i++)
		  tIntensities.push_back(DistMap.GetPointer()->GetPixel(testIndices[i]));


	  //for (unsigned int i = 0; i < testIndices.size(); i++)
	  //{
	  //  ImageType::IndexType testIndex = testIndices[i];
	  //  double min = 0;

	  //  for (unsigned int j = 0; j < tumorIndices.size(); j++)
	  //  {
	  //    ImageType::IndexType tumorIndex = tumorIndices[j];
	  //    double distance = std::sqrt(std::pow(testIndex[0] - tumorIndex[0], 2) + std::pow(testIndex[1] - tumorIndex[1], 2) + std::pow(testIndex[2] - tumorIndex[2], 2));
	  //    if (j == 0)
	  //      min = distance;
	  //    else
	  //    {
	  //      if (distance < min)
	  //        min = distance;
	  //    }
	  //  }
	  //  tIntensities.push_back(min);

  }
  return testIndices;
}