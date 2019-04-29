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

//#include "CAPTk.h"
//#include "CaPTkUtils.h"
#include "CaPTkDefines.h"
#include "CaPTkUtils.h"
#include "vtkSmartPointer.h"
#include "itkImageSeriesReader.h"
#include "itkGDCMSeriesFileNames.h"
#include "vtkTransform.h"
#include <QObject>

class Slicer;
class Landmarks;
class vtkInteractorStyle;
class vtkImageData;
//class itkImageSeriesReader;
class vtkRenderWindow;
class vtkTransform;

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

  //! Get/Set Comparison mode
  void SetComparisonMode(bool mode);
  bool GetComparisonMode();

  //! set image on a specific slicer
  void SetImageOnSlicer(vtkSmartPointer< vtkImageData >, vtkTransform* transform, int slicer);

  //! Get image as ITK image
  ImageTypeFloat3D::Pointer GetITKImage();

  void UpdateVtkImage();

  void SetImage(ImageTypeFloat3D::Pointer t1ceimage);

  void SetOriginalOrigin(std::vector< double > origin);

  inline void SetOriginalDirection(ImageTypeFloat3D::DirectionType direction)
  {
    mDirection = direction;
  }

  void SetPerfImage(ImageTypeFloat4D::Pointer image);

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

  inline void SetId(const std::string &id)
  {
    mId = id;
  }
  inline std::string GetId()
  {
    return mId;
  }

  int GetDimension();

  inline void SetOrder(int order)
  {
    mOrder = order;
  }
  inline int GetOrder()
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
  inline void SetThresholdIndex(double i)
  {
    mThresholdIndex = i;
  }

  double GetColorWindow();
  double GetColorLevel();
  inline int GetPreset()
  {
    return mPreset;
  }
  inline double GetThresholdIndex()
  {
    return mThresholdIndex;
  }

  void UpdateViews(int slicer);
  void UpdateLinked(int slicer);
  void updateToRefCam(Slicer *refSlicer);
  void ResetTransformationToIdentity();
  void Render();

  void AddLink(const std::string &newId);

  void RemoveLink(const std::string &oldId);

  inline bool IsLinked()
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
  inline void setTempFolderLocation(const std::string &tempFolder)
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
  bool m_ComparisonMode; //! comparison mode

};
#if __GNUC__
#pragma GCC visibility pop
#endif

#endif
