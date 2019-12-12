/**
\file  Slicer.h

\brief Implementation of Slicer class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _Slicer_h_
#define _Slicer_h_


//#include "CAPTk.h"
//#include "Landmarks.h"
#include "QBorderWidget.h"
#include "vtkImageViewer2.h"

class Landmarks;
class vtkImageReslice;
class vtkCursor2D;
class vtkImageMapToColors;
class vtkBox;
class vtkCursor3D;
class vtkGlyph3D;
class vtkImageMapToColors;
class vtkClipPolyData;
class vtkPolyDataMapper;
class vtkVertexGlyphFilter;
class vtkLabeledDataMapper;
class vtkRegularPolygonSource;
class vtkInteractorStyle;
class vtkCornerAnnotation;

class Slicer : public vtkImageViewer2
{
public:
  Slicer();
  ~Slicer();

  static Slicer *New();
#if VTK_MAJOR_VERSION >= 6
  vtkTypeMacro(Slicer, vtkImageViewer2);
#else
  vtkTypeRevisionMacro(Slicer, vtkImageViewer2);
#endif

  //! set interactor style
  void SetInteractorStyle(vtkInteractorStyle* style);

  void SetImage(vtkImageData* image, vtkTransform* transform);
  vtkImageData* GetImage() {
    return mImage;
  }
  vtkTransform* GetTransform() {
    return mTransform;
  }
  void SetOverlay(vtkImageData* overlay);
  vtkImageData* GetOverlay()
  {
    return mOverlay;
  }
  void SetOverlayOpacity(double opacity);
  void RemoveOverlay();
  void SetMask(vtkImageData* mask);
  vtkImageData* GetMask()
  {
    return mMask;
  }
  void SetMaskOpacity(double opacity);
  double GetMaskOpacity();
  void RemoveMask();

  void ResetMap();

  void AddLabelToMap(int label);
  void ShowLabelOnMap(int label);

  void SetLandmarks(Landmarks* landmarks, int type);
  void SetSliceOrientation(int orientation);
  void SetSlice(int s);

  void SetOpacity(double s);
  void SetRenderWindow(int orientation, vtkRenderWindow * rw);
  void SetDisplayMode(bool i);
  void FlipHorizontalView();
  void FlipVerticalView();
  static double GetScalarComponentAsDouble(vtkSmartPointer< vtkImageData > image, double X, double Y, double Z, int &ix, int &iy, int &iz, int component = 0);
  void Render();
  void ResetCamera();

  //! Get/Set Comparison mode
  void SetComparisonMode(bool mode);
  bool GetComparisonMode();

  void SetInitPosition();
  void SetCurrentPosition(double x, double y, double z);
  double* GetCurrentPosition() {
	  return mCursor;
  }
  double* GetCursorPosition() {
    return mCursor;
  }
  void UpdateCursorPosition();
  void SetCursorVisibility(bool s);
  bool GetCursorVisibility();
  void SetActive(bool active);
  bool GetActive();
  void SetCursorColor(double r, double g, double b);
  bool GetCornerAnnotationVisibility();
  void SetLandmarksVisibility(bool s);
  bool GetLandmarksVisibility();

  void UpdateLandmarks();
  void ForceUpdateDisplayExtent();

  virtual void SetColorWindow(double s);
  virtual void SetColorLevel(double s);

  int* GetDisplayExtent();
  int GetOrientation();
  void UpdateOrientation();
  void UpdateDisplayExtent();

  void AdjustResliceToSliceOrientation(vtkImageReslice *reslice);
  void ConvertImageToImageDisplayExtent(vtkSmartPointer< vtkImageData > sourceImage, const int sourceExtent[6], vtkSmartPointer< vtkImageData > targetImage, int targetExtent[6]);
  void ClipDisplayedExtent(int extent[6], int refExtent[6]);

  //! set image series description to be displayed on slicer
  void SetImageSeriesDescription(std::string description);

  bool mActive;


  int mCornerAnnotationVisibility;
  int mLandmarksVisibility;

  vtkSmartPointer<vtkImageData> mImage;
  vtkSmartPointer<vtkImageData> mOverlay;
  vtkSmartPointer<vtkImageData> mMask, mMaskOriginal; // mMaskOriginal is used only when the mask manipulation needs to happen
  vtkSmartPointer<vtkTransform> mTransform;
  Landmarks* mLandmarks;
  int mLandmarksType;
  bool mOriginalMaskSaved = false;

  vtkSmartPointer<vtkImageReslice> mImageReslice;
  vtkSmartPointer<vtkCursor2D> crossCursor;
  vtkSmartPointer<vtkPolyDataMapper2D> pdm;
  vtkSmartPointer<vtkActor2D> pdmA;
  //
  vtkSmartPointer<vtkImageReslice> mOverlayReslice;
  vtkSmartPointer<vtkImageActor> mOverlayActor;
  vtkSmartPointer<vtkImageMapToWindowLevelColors> mOverlayMapper;
  double mOverlayOpacity;
  //
  vtkSmartPointer<vtkImageReslice> mMaskReslice;
  vtkSmartPointer<vtkImageActor> mMaskActor;
  vtkSmartPointer<vtkImageMapToColors> mMaskMapper;
  double mMaskOpacity;
  //
  vtkSmartPointer<vtkBox> mClipBox;
  //
  vtkSmartPointer<vtkCursor3D> mCross;
  vtkSmartPointer<vtkGlyph3D> mLandGlyph;
  vtkSmartPointer<vtkClipPolyData> mLandClipper;
  vtkSmartPointer<vtkPolyDataMapper> mLandMapper;
  vtkSmartPointer<vtkActor> mLandActor;
  //
  vtkSmartPointer<vtkVertexGlyphFilter> mLandLabelGlyph;
  vtkSmartPointer<vtkClipPolyData> mLandLabelClipper;
  vtkSmartPointer<vtkLabeledDataMapper> mLandLabelMapper;
  vtkSmartPointer<vtkActor2D> mLandLabelActor;
  //
  vtkSmartPointer<vtkRegularPolygonSource> mCircle;
  vtkSmartPointer<vtkGlyph3D> mLandRadiusGlyph;
  vtkSmartPointer<vtkClipPolyData> mLandRadiusClipper;
  vtkSmartPointer<vtkPolyDataMapper> mLandRadiusMapper;
  vtkSmartPointer<vtkActor> mLandRadiusActor;


  double mCursor[3];

  vtkSmartPointer<vtkBorderWidget> borderWidget;
  vtkBorderCallback * borderCallback;

  //!slicer knows if the viewing mode is comparison mode or not
  bool m_ComparisonMode; 

  vtkSmartPointer<vtkCornerAnnotation> mCornerAnnotation;
};


#endif
