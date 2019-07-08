/**
\file  SlicerManager.h

\brief Declaration of SlicerManager class

https://www.med.upenn.edu/sbia/software/ <br>
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

  //! Default contructor
  SlicerManager(int numberOfSlicers, Landmarks* landmarks, Landmarks* seedPoints, Landmarks* tissuePoints);

  //! Default destructor
  ~SlicerManager();

  //! Get last error
  std::string GetLastError()
  {
    return mLastError;
  }

  //! Get/Set Comparison mode
  void SetComparisonMode(bool mode);

  //! Get/Set Comparison mode
  bool GetComparisonMode();

  //! set image on a specific slicer
  void SetImageOnSlicer(vtkSmartPointer< vtkImageData >, vtkTransform* transform, int slicer);

  //! Get image as ITK image
  ImageTypeFloat3D::Pointer GetITKImage();

  //! Update VTK image
  void UpdateVtkImage();

  //! Set new ITK image
  void SetImage(ImageTypeFloat3D::Pointer t1ceimage);

  //! Set the original Origin (since we are setting directions to identity) using a std::vector
  void SetOriginalOrigin(std::vector< double > origin);

  //! Set the original Origin (since we are setting directions to identity) using itk::Point
  void SetOriginalOrigin(ImageTypeFloat3D::PointType origin);

  //! Set the original direction (since we are setting directions to identity) 
  inline void SetOriginalDirection(ImageTypeFloat3D::DirectionType direction)
  {
    mDirection = direction;
  }

  //! Set new 4D image 
  void SetPerfImage(ImageTypeFloat4D::Pointer image);

  //! Set the vtk mask
  void SetMask(vtkSmartPointer< vtkImageData >  mask);

  //! Get the mask
  vtkSmartPointer< vtkImageData > GetMask()
  {
    return mMask;
  }

  //! Get the full file path of input image
  std::string GetPathFileName()
  {
    return mPathFileName;
  }
  //! Get the full file name of input image
  std::string GetFileName()
  {
    return mFileName;
  }
  //! Get the base file name of input image
  std::string GetBaseFileName()
  {
    return mBaseFileName;
  }

  //! Enable/Disable interpolation for all vktImage actors
  void ToggleInterpolation();

  //! Get the Slicer from the index
  Slicer* GetSlicer(int i);

  //! Update selected Slicer
  void UpdateSlicer(int num, bool state);

  //! Set Slicer Window for selected Slicer with the defined vtkRenderWindow
  void SetSlicerWindow(int i, vtkRenderWindow* RW);

  //! Set Slicer style navigator for selected Slicer with the defined vtkInteractorStyle
  void SetInteractorStyleNavigator(int i, vtkInteractorStyle* style);

  //! Get the total number of Slicers
  int GetNumberOfSlicers()
  {
    return mSlicers.size();
  }

  //! Get the vtkImage
  vtkSmartPointer< vtkImageData >  GetImage()
  {
    return mImage;
  }

  //! Set ID of the SlicerManager
  inline void SetId(const std::string &id)
  {
    mId = id;
  }

  //! Get ID of the SlicerManager
  inline std::string GetId()
  {
    return mId;
  }

  //! Get Dimension (always returns 3 for CaPTk)
  int GetDimension();

  //! Set order for the Slicer
  inline void SetOrder(int order)
  {
    mOrder = order;
  }

  //! Get order for the Slicer
  inline int GetOrder()
  {
    return mOrder;
  }

  //! Set the file name of the input image
  void SetFilename(const std::string &f);

  //! Set the orientation of the selected slicer
  void SetSliceOrientation(int slicer, int orientation);
  //void SetTSlice(int slice);
  //void SetNextTSlice(int originating_slicer);
  //void SetPreviousTSlice(int originating_slicer);
  //void SetTSliceInSlicer(int tslice, int slicer);

  //! Generate default LUT
  void GenerateDefaultLookupTable();

  //! Set the color window
  void SetColorWindow(double s);

  //! Set the color level
  void SetColorLevel(double s);

  //! Get the color level
  double GetColorLevel();

  //! Get the color level
  double GetColorWindow();

  //! Set the opacity for the selected slicer
  void SetOpacity(int i, double factor);

  //! Set the rendering preset
  void SetPreset(int preset);

  //! Get the rendering preset
  inline int GetPreset()
  {
    return mPreset;
  }

  //! Set the threshold index 
  inline void SetThresholdIndex(double i)
  {
    mThresholdIndex = i;
  }

  //! Get the threshold index 
  inline double GetThresholdIndex()
  {
    return mThresholdIndex;
  }

  //! Update views for selected Slicer
  void UpdateViews(int slicer);

  //! Update linked views for selected Slicer
  void UpdateLinked(int slicer);

  //! Update views to the reference camera
  void updateToRefCam(Slicer *refSlicer);

  //! Update all transforms to identity
  void ResetTransformationToIdentity();

  //! Render all Slicers
  void Render();

  //! Add Link ID
  void AddLink(const std::string &newId);

  //! Remove Link ID
  void RemoveLink(const std::string &oldId);

  //! Check link status
  inline bool IsLinked()
  {
    return mLinkedId.size() > 0;
  }

  //! Remove actors from all Slicers
  void RemoveActors();

  //! Check if Slicer is activated on image panel - used to check when an image focus get changed
  void Activated();

  //! Check if Slicer is picked on image panel
  void Picked();

  //! Update image information on current cursor position
  void UpdateInfoOnCursorPosition(int slicer);

  //! Update Window Level
  void UpdateWindowLevel();

  //! Update visualized slice for selected Slicer
  void UpdateSlice(int slicer);

  //! Update slice range for selected Slicer
  void UpdateSliceRange(int slicer);

  //! Add tissue point (no radius information)
  void AddLandmark(float x, float y, float z);
  //void AddLoadedLandmark(double x, double y,double z,int landmarktype);

  //! Remove tissue point (no radius information)
  void RemoveLandmark(float x, float y, float z);

  //! Add tumor point (with radius information taken from the second point initialized with Ctrl+Space)
  void AddLandmarkRadius(float x, float y, float z);

  //! Shift tumor point 
  void AddLandmarkShift(float x, float y, float z);

  //! Change selected image to the one defined with order
  void NextImageWithOrder(int order);

  //! Left button release event
  void LeftButtonReleaseEvent(int slicer);

  //! Function to control change in the vertical sliders
  void VerticalSliderHasChanged(int slicer, int slice);

  //! Get scalar component as double
  double GetScalarComponentAsDouble(vtkSmartPointer< vtkImageData > image, double X, double Y, double Z, int component = 0);

  //! Set current landmarks type
  void SetCurrentLandmarksType(int type, int row, int col);

  //! Set temporary folder location
  inline void setTempFolderLocation(const std::string &tempFolder)
  {
    m_tempFolderLocation = tempFolder;
  }

  //! Update border coordinates for x and y axes
  void UpdateBorderCoordinates(double startX, double startY, double endX, double endY);

  //! Update border coordinates for the z axis
  void UpdateBorderCoordinates(double startZ, double endZ);

  //! Add an action - used for undo
  void ActionAdded(std::vector<PointVal>& points);

  //! Get the 3D image at the current index for 4D image
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

  Landmarks* mLandmarks;
  Landmarks* mSeedPoints;
  Landmarks* mTissuePoints;
  int mCurrentLandmarksType;
  int mCurrentLandmarksRow;
  int mCurrentLandmarksCol;
  int mTissueSelectedEntry;

  std::string m_tempFolderLocation;
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
