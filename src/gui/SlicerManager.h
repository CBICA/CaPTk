/**
\file  SlicerManager.h

\brief Declaration of SlicerManager class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _SlicerManager_h_
#define _SlicerManager_h_


#include "CAPTk.h"
#include "Slicer.h"
#include "InteractorStyleNavigator.h"
#include "PreprocessingPipelineClass.h"
#include "itkImageToVTKImageFilter.h"
#include "itkExtractImageFilter.h"

#define PRESET_AUTO 0
#define PRESET_USER 1
#define PRESET_LABEL 2
#define PRESET_LABEL2 3
#define PRESET_THRESHOLD 4
#define PRESET_PROB 5
#define PRESET_GEODESIC 6

/**
\class SlicerManager

\brief This class handles the 3 different slices being visualized from an abstract level

*/
#if __GNUC__
#pragma GCC visibility push(hidden)
#endif
class SlicerManager : public QObject
{
  Q_OBJECT

public:
  SlicerManager(int numberOfSlicers, Landmarks* landmarks, Landmarks* seedPoints, Landmarks* tissuePoints);
  ~SlicerManager();

  std::string GetLastError()
  {
    return mLastError;
  }

  void UpdateVtkImage()
  {
    typedef itk::Image<float, 3> OutputImageType3D;
    // Convert from ITK object to VTK object
    typedef itk::ImageToVTKImageFilter<OutputImageType3D> ConverterType;
    ConverterType::Pointer converter = ConverterType::New();
    mItkToVtkConverters = dynamic_cast< itk::ProcessObject *>(converter.GetPointer());
    converter->SetInput(mITKImage);
    converter->Update();
    mImage = converter->GetOutput();

    vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
    matrix->Identity();
    // for (unsigned int i = 0; i < mITKImage->GetImageDimension(); i++)
    // {
    //   for (unsigned int j = 0; j < mITKImage->GetImageDimension(); j++)
    //   {
    //     (*matrix)[i][j] = mITKImage->GetDirection()[i][j];
    //     (*matrix)[i][3] -= (*matrix)[i][j] * mITKImage->GetOrigin()[j];
    //   }
    //   (*matrix)[i][3] += mITKImage->GetOrigin()[i];
    // }
    for (unsigned int i = 0; i < mITKImage->GetImageDimension(); i++)
    {
      float temp = 1;
      for (unsigned int j = 0; j < mITKImage->GetImageDimension(); j++)
      {
        //(*matrix)[i][j] = mITKImage->GetDirection()[i][j];
        //(*matrix)[i][3] -= (*matrix)[i][j] * mITKImage->GetOrigin()[j];
        matrix->SetElement(i, j, mITKImage->GetDirection()[i][j]);
        temp -= matrix->GetElement(i, j) * mITKImage->GetOrigin()[j];
      }
      //(*matrix)[i][3] += mITKImage->GetOrigin()[i];
      temp += mITKImage->GetOrigin()[i];
      matrix->SetElement(i, 3, temp);
    }

    matrix->Invert();
    mTransform->SetMatrix(matrix);
    for (int i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetImage(mImage, mTransform);
      mSlicers[i]->GetImageActor()->SetInterpolate(0);
    }
  }

  void SetImage(ImageTypeFloat3D::Pointer t1ceimage)
  {
    mITKImage = t1ceimage;
    UpdateVtkImage();

  }

  void SetOriginalOrigin(std::vector< double > origin)
  {
    mOrigin[0] = origin[0];
    mOrigin[1] = origin[1];
    if (origin.size() == 2)
    {
      mOrigin[2] = 0;
    }
    else
    {
      mOrigin[2] = origin[2];
    }
  }

  void SetOriginalDirection(ImageTypeFloat3D::DirectionType direction)
  {
    mDirection = direction;
  }

  void SetPerfImage(ImageTypeFloat4D::Pointer image)
  {
    mPerfusionImagePointer = image;
    //crop one 3D image from the 4D series
    ImageTypeFloat4D::RegionType region = image->GetLargestPossibleRegion();
    ImageTypeFloat4D::IndexType regionIndex;
    ImageTypeFloat4D::SizeType regionSize;
    regionSize[0] = region.GetSize()[0];
    regionSize[1] = region.GetSize()[1];
    regionSize[2] = region.GetSize()[2];
    regionSize[3] = 0;
    regionIndex[0] = 0;
    regionIndex[1] = 0;
    regionIndex[2] = 0;
    regionIndex[3] = 0;
    ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);

    typedef itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D > FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetExtractionRegion(desiredRegion);
    filter->SetInput(mPerfusionImagePointer);

    filter->SetDirectionCollapseToIdentity(); // This is required.
    filter->Update();
    mITKImage = filter->GetOutput();
    UpdateVtkImage();
  }


  void SetMask(vtkSmartPointer< vtkImageData >  mask);
  vtkSmartPointer< vtkImageData >  GetMask()
  {
    return mMask;
  }

  std::string GetPathFileName()
  {
    return mPathFileName;
  }
  std::string GetFileName()
  {
    return mFileName;
  }
  std::string GetBaseFileName()
  {
    return mBaseFileName;
  }

  void ToggleInterpolation();
  Slicer* GetSlicer(int i);
  void UpdateSlicer(int num, bool state);
  void SetSlicerWindow(int i, vtkRenderWindow* RW);
  void SetInteractorStyleNavigator(int i, vtkInteractorStyle* style);

  int GetNumberOfSlicers()
  {
    return mSlicers.size();
  }
  vtkSmartPointer< vtkImageData >  GetImage()
  {
    return mImage;
  }

  void SetId(const std::string &id)
  {
    mId = id;
  }
  std::string GetId()
  {
    return mId;
  }

  int GetDimension()
  {
    if (mImage) {
      int dim = 0;
#if VTK_MAJOR_VERSION <= 5
      int* extent = mImage->GetWholeExtent();
#else
      int* extent = mImage->GetExtent();
#endif
      if (extent[4] != extent[5]) {
        dim = 3;
      }
      else if (extent[3] != extent[4]) {
        dim = 2;
      }
      else if (extent[0] != extent[1]) {
        dim = 1;
      }
      return dim;
    }
    else {
      return -1;
    }
  }
  void SetOrder(int order)
  {
    mOrder = order;
  }
  int GetOrder()
  {
    return mOrder;
  }

  void SetFilename(const std::string &f);
  void SetSliceOrientation(int slicer, int orientation);
  void SetTSlice(int slice);
  void SetNextTSlice(int originating_slicer);
  void SetPreviousTSlice(int originating_slicer);
  void SetTSliceInSlicer(int tslice, int slicer);

  void GenerateDefaultLookupTable();
  void SetColorWindow(double s);
  void SetColorLevel(double s);
  void SetOpacity(int i, double factor);
  void SetPreset(int preset);
  void SetThresholdIndex(double i)
  {
    mThresholdIndex = i;
  }

  double GetColorWindow();
  double GetColorLevel();
  int GetPreset()
  {
    return mPreset;
  }
  double GetThresholdIndex()
  {
    return mThresholdIndex;
  }

  void UpdateViews(int slicer);
  void UpdateLinked(int slicer);
  void updateToRefCam(Slicer *refSlicer);
  void ResetTransformationToIdentity();
  void Render();

  void AddLink(const std::string &newId)
  {
    mLinkedId.push_back(newId);
    mLinkedId.unique();
  }
  void RemoveLink(const std::string &oldId)
  {
    mLinkedId.remove(oldId);
  }
  bool IsLinked()
  {
    return mLinkedId.size() > 0;
  }

  void RemoveActors();
  void Activated();
  void Picked();
  void UpdateInfoOnCursorPosition(int slicer);
  void UpdateInfoOnCurrentPosition(int slicer);
  void UpdateWindowLevel();
  void UpdateSlice(int slicer);
  void UpdateSliceRange(int slicer);

  void AddLandmark(float x, float y, float z);
  //void AddLoadedLandmark(double x, double y,double z,int landmarktype);
  void RemoveLandmark(float x, float y, float z);
  void AddLandmarkRadius(float x, float y, float z);
  void AddLandmarkShift(float x, float y, float z);

  void NextImageWithOrder(int order);
  void LeftButtonReleaseEvent(int slicer);
  void VerticalSliderHasChanged(int slicer, int slice);
  double GetScalarComponentAsDouble(vtkSmartPointer< vtkImageData > image, double X, double Y, double Z, int component = 0);

  void SetCurrentLandmarksType(int type, int row, int col);

  void EraseCompleteNearDrawing();
  void EraseCompleteFarDrawing();

  Landmarks* mLandmarks;
  Landmarks* mSeedPoints;
  Landmarks* mTissuePoints;
  int mCurrentLandmarksType;
  int mCurrentLandmarksRow;
  int mCurrentLandmarksCol;
  int mTissueSelectedEntry;

  std::string m_tempFolderLocation;
  void setTempFolderLocation(const std::string &tempFolder)
  {
    m_tempFolderLocation = tempFolder;
  }

  std::string GetFileNameInADicomDirectory(QString directoryname);
  void UpdateBorderCoordinates(double startX, double startY, double endX, double endY);
  void UpdateBorderCoordinates(double startZ, double endZ);
  void ActionAdded(std::vector<PointVal>& points);
  void Get3DImageAtCurrentPerfusionIndex(int index);




signals:
  void currentImageChanged(std::string &id);
  void currentPickedImageChanged(std::string id);
  void UpdatePosition(int visibility, double x, double y, double z, double X, double Y, double Z, double value);
  void UpdateOrientation(int slicer, int orientation);
  void UpdateSlice(int slicer, int slice);
  void UpdateSliceRange(int slice, int min, int max);
  void WindowLevelChanged();
  void UpdateLinkManager(std::string, int slicer, double x, double y, double z);
  void LandmarkAdded();
  void UpdateNumberOfPoints(int drawn_points, VectorDouble tissueTypeCounter);
  void SeedPointsAdded();
  void SeedPointsAdded(int index, bool update);
  void TissuePointsAdded(int index);
  void HighlightTableIndexOnDeletion(double x, double y, double z, double x1, double y1, double z1, double value);
  void TissuePointsRemoved(int index);
  void ChangeImageWithOrder(SlicerManager *sm, int order);
  void LeftButtonReleaseSignal(int slicer);
  void AVerticalSliderHasChanged(int slicer, int slice);



  void UpdateNumberOfPoints(int nearpoints, int farpoints);
  void UpdateBorderWidgetInMain(double startX, double startY, double endX, double endY);
  void UpdateBorderWidgetInMain(double startZ, double endZ);
  void UpdateActionInMain(const QVariantList& points);




public:
  vtkSmartPointer<vtkImageData> mImage;
  ImageTypeFloat3D::Pointer mITKImage;
  ImageTypeFloat3D::DirectionType mDirection;
  ImageTypeFloat3D::PointType mOrigin;
  ImageTypeFloat4D::Pointer mPerfusionImagePointer;
  std::vector< ImageTypeFloat3D::Pointer > mPerfusionImagePointerDicom;
  vtkSmartPointer<vtkTransform> mTransform;
  vtkSmartPointer<vtkImageData> mMask;
  itk::ImageSeriesReader< ImageTypeFloat3D >::DictionaryArrayRawPointer mDicomDictionaryArray;
  itk::GDCMSeriesFileNames::FileNamesContainerType mFileNames;
  itk::ImageSeriesReader< ImageTypeShort3D >::Pointer mSeriesReader;
  typedef itk::ProcessObject::Pointer ConverterPointer;
  ConverterPointer mItkToVtkConverters;
  std::vector<vtkSmartPointer<Slicer> > mSlicers;

  int mPreset;
  std::string mPathFileName;
  int mImageType;
  int mImageSubType;
  std::string mFileName;
  std::string mBaseFileName;
  int mBaseFileNameNumber;
  std::string mId;
  std::string mLastError;
  std::list<std::string> mLinkedId;

  std::vector<int> mPreviousSlice;
  std::vector<int> mPreviousTSlice;

  int mOrder;
  double mThresholdIndex;
  double mThreshold;

};
#if __GNUC__
#pragma GCC visibility pop
#endif

#endif
