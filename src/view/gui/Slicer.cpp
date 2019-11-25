/**
\file  Slicer.cpp

\brief Implementation of Slicer class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "Slicer.h"
#include "SlicerManagerCommand.h"
#include "Landmarks.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkInformation.h"
#include "vtkCursor2D.h"
#include "vtkPolyDataMapper2D.h"
#include "vtkActor2D.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkImageMapper3D.h"
#include "vtkImageReslice.h"
#include "vtkProperty2D.h"
#include "vtkImageData.h"
#include "vtkTransform.h"
#include "vtkAbstractTransform.h"
#include "vtkLookupTable.h"
#include "vtkScalarsToColors.h"
#include "vtkCursor3D.h"
#include "vtkRenderWindow.h"
#include "vtkGlyph3D.h"
#include "vtkClipPolyData.h"
#include "vtkBox.h"
#include "vtkImplicitFunction.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkVertexGlyphFilter.h"
#include "vtkLabeledDataMapper.h"
#include "vtkRegularPolygonSource.h"
#include "vtkTextProperty.h"
#include "vtkCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkRenderWindowInteractor.h"
#include <vtkCornerAnnotation.h>

///// debug
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
///// debug

vtkStandardNewMacro(Slicer);
#if VTK_MAJOR_VERSION < 6
  vtkCxxRevisionMacro(Slicer, "$Revision: 1.70 $");
#endif
Slicer::Slicer()
{
  this->UnInstallPipeline();
  mImage = NULL;



  mCursor[0] = -VTK_DOUBLE_MAX;
  mCursor[1] = -VTK_DOUBLE_MAX;
  mCursor[2] = -VTK_DOUBLE_MAX;

  crossCursor = vtkSmartPointer<vtkCursor2D>::New();
  crossCursor->AllOff();
  crossCursor->AxesOn();
  crossCursor->SetTranslationMode(0);
  crossCursor->SetRadius(0);


  mCornerAnnotationVisibility = 1;
  mLandmarksVisibility = 1;

  pdm = vtkSmartPointer<vtkPolyDataMapper2D>::New();
#if VTK_MAJOR_VERSION <= 5
  pdm->SetInput(crossCursor->GetOutput());
#else
  pdm->SetInputConnection(crossCursor->GetOutputPort());
  pdm->Update();
#endif

  pdmA = vtkSmartPointer<vtkActor2D>::New();
  pdmA->SetMapper(pdm);
  // red
  pdmA->GetProperty()->SetColor(1.0, 0.1, 0.1);
  // green
  //pdmA->GetProperty()->SetColor(0.1, 1.0, 0.1);
  pdmA->GetProperty()->SetOpacity(0.5);
  pdmA->SetVisibility(1);
  pdmA->SetPickable(0);



  mLandmarks = NULL;
  mLandmarksType = LANDMARK_TYPE::NONE;

  mActive = false;

  mOverlay = NULL;
  mMask = NULL;
  mMaskOriginal = NULL;

  mOverlayOpacity = 0.5;
  mMaskOpacity = 0.5;

  this->WindowLevel = vtkImageMapToWindowLevelColors::New();

  this->InstallPipeline();

#if VTK_MAJOR_VERSION >= 6 || (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 10)
  this->GetImageActor()->GetMapper()->BorderOn();
#endif



  borderWidget = vtkSmartPointer<vtkBorderWidget>::New();
  borderCallback = new vtkBorderCallback();

  this->mCornerAnnotation = vtkSmartPointer<vtkCornerAnnotation>::New();
  this->mCornerAnnotation->GetTextProperty()->SetColor(1, 1, 1);
  this->Renderer->AddViewProp(mCornerAnnotation);
}

void Slicer::SetActive(bool active)
{
  mActive = active;
  if (mActive) {
    SetCursorColor(0.1, 1.0, 0.1);
  }
  else {
    SetCursorColor(1.0, 0.1, 0.1);
  }
}
bool Slicer::GetActive()
{
  return mActive;
}

void Slicer::SetCursorColor(double r, double g, double b)
{
  pdmA->GetProperty()->SetColor(r, g, b);
}

void Slicer::SetCursorVisibility(bool bVis)
{
  pdmA->SetVisibility(bVis);
}
bool Slicer::GetCursorVisibility()
{
  bool bVis = (bool)pdmA->GetVisibility();
  return bVis;
}

bool Slicer::GetCornerAnnotationVisibility()
{
  return mCornerAnnotationVisibility;
}
void Slicer::SetLandmarksVisibility(bool s)
{
  mLandmarksVisibility = s;
  if (mLandActor)
    mLandActor->SetVisibility(mLandmarksVisibility);
}
bool Slicer::GetLandmarksVisibility()
{
  return mLandmarksVisibility;
}

Slicer::~Slicer()
{

}

void Slicer::SetInitPosition()
{
  double bX0 = this->GetImageActor()->GetDisplayExtent()[0];
  double bX1 = this->GetImageActor()->GetDisplayExtent()[1];
  double bY0 = this->GetImageActor()->GetDisplayExtent()[2];
  double bY1 = this->GetImageActor()->GetDisplayExtent()[3];
  double bZ0 = this->GetImageActor()->GetDisplayExtent()[4];
  double bZ1 = this->GetImageActor()->GetDisplayExtent()[5];
  double cX, cY, cZ;

  cX = (int)((bX0 + bX1) / 2);
  cY = (int)((bY0 + bY1) / 2);
  cZ = (int)((bZ0 + bZ1) / 2);

  mCursor[0] = cX * this->GetInput()->GetSpacing()[0] + this->GetInput()->GetOrigin()[0];
  mCursor[1] = cY * this->GetInput()->GetSpacing()[1] + this->GetInput()->GetOrigin()[1];
  mCursor[2] = cZ * this->GetInput()->GetSpacing()[2] + this->GetInput()->GetOrigin()[2];
}

void Slicer::SetCurrentPosition(double x, double y, double z)
{
	mCursor[0] = x;
	mCursor[1] = y;
	mCursor[2] = z;
}

void Slicer::SetInteractorStyle(vtkInteractorStyle * style)
{
  this->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
}

void Slicer::SetImage(vtkImageData* image, vtkTransform* transform)
{
  if (image != NULL && transform != NULL) {
    mImage = image;
    mTransform = transform;

    if (!mImageReslice) {
      mImageReslice = vtkSmartPointer<vtkImageReslice>::New();
      //mImageReslice->SetInterpolationModeToLinear();
      mImageReslice->SetInterpolationModeToNearestNeighbor();
      mImageReslice->AutoCropOutputOn();
      mImageReslice->SetBackgroundColor(0, 0, 0, 1);
    }
    mImageReslice->SetResliceTransform(transform);
#if VTK_MAJOR_VERSION <= 5
    mImageReslice->SetInput(0, image);
    mImageReslice->UpdateInformation();
#else
    mImageReslice->SetInputData(0, image);
    mImageReslice->Update();
#endif

#if VTK_MAJOR_VERSION <= 5
    this->Superclass::SetInput(mImageReslice->GetOutput());
#else
    this->Superclass::SetInputConnection(mImageReslice->GetOutputPort());
#endif

    int extent[6];
#if VTK_MAJOR_VERSION <= 5
    this->GetInput()->GetWholeExtent(extent);
#else
    this->GetInput()->GetExtent(extent);
#endif

    if (Slice < extent[SliceOrientation * 2] || Slice >= extent[SliceOrientation * 2 + 1]) {
      Slice = (int)((extent[SliceOrientation * 2 + 1] + extent[SliceOrientation * 2]) / 2.0);
    }

    extent[SliceOrientation * 2] = Slice;
    extent[SliceOrientation * 2 + 1] = Slice;
#if VTK_MAJOR_VERSION <= 5
    mImageReslice->GetOutput()->SetUpdateExtent(extent);
    mImageReslice->GetOutput()->Update();
#else
    //mImageReslice->SetUpdateExtent(extent);
    mImageReslice->UpdateExtent(extent);
    mImageReslice->Update();
#endif

    this->UpdateDisplayExtent();
  }


}

void Slicer::SetMask(vtkImageData* mask)
{
  if (mask != NULL)
  {
    mMask = mask;

    if (!mMaskReslice)
    {
      mMaskReslice = vtkSmartPointer<vtkImageReslice>::New();
      //mMaskReslice->SetInterpolationModeToLinear();
      mMaskReslice->SetInterpolationModeToNearestNeighbor();
      mMaskReslice->AutoCropOutputOn();
      mMaskReslice->SetBackgroundColor(0, 0, 0, 0);
    }
    mMaskReslice->SetResliceTransform(mTransform);
#if VTK_MAJOR_VERSION <= 5
    mMaskReslice->SetInput(0, mask);
    //mMaskReslice->UpdateInformation();
#else
    mMaskReslice->SetInputData(0, mask);
    mMaskReslice->Update();
#endif

    vtkSmartPointer<vtkLookupTable> LUT = vtkSmartPointer<vtkLookupTable>::New();
    LUT->SetNumberOfTableValues(10); //TBD avoid hard coding
    LUT->SetTableRange(0, 9);
    
    // background - kept for consistency in masks instead of using the same for erase
    LUT->SetTableValue(0, 0, 0, 0, 0);
    // label_1 -- RED
    LUT->SetTableValue(1, 1, 0, 0, mMaskOpacity);
    // label_2 -- GREEN
    LUT->SetTableValue(2, 0, 1, 0, mMaskOpacity);
    // label_3 -- YELLOW
    LUT->SetTableValue(3, 1, 1, 0, mMaskOpacity);
    // label_4 -- BLUE
    LUT->SetTableValue(4, 0, 0, 1, mMaskOpacity);
    // label_5 -- MAGENTA
    LUT->SetTableValue(5, 1, 0, 1, mMaskOpacity);
    // label_5 -- CYAN
    LUT->SetTableValue(6, 0, 1, 1, mMaskOpacity);
    // label_5 -- RED_2
    LUT->SetTableValue(7, 1, 0.5, 0.5, mMaskOpacity);
    // label_5 -- GREEN_2
    LUT->SetTableValue(8, 0.5, 1, 0.5, mMaskOpacity);
    // label_5 -- BLUE_2
    LUT->SetTableValue(9, 0.5, 0.5, 1, mMaskOpacity);



    if (!mMaskMapper)
      mMaskMapper = vtkSmartPointer<vtkImageMapToColors>::New();

    mMaskMapper->SetLookupTable(LUT);
    mMaskMapper->PassAlphaToOutputOn();
#if VTK_MAJOR_VERSION <= 5
    mMaskMapper->SetInput(mMaskReslice->GetOutput());
#else
    mMaskMapper->SetInputConnection(mMaskReslice->GetOutputPort());
    mMaskMapper->Update();
#endif

    if (!mMaskActor)
    {
      mMaskActor = vtkSmartPointer<vtkImageActor>::New();
#if VTK_MAJOR_VERSION <= 5
      mMaskActor->SetInput(mMaskMapper->GetOutput());
#else
      mMaskActor->GetMapper()->SetInputConnection(mMaskMapper->GetOutputPort());
      mMaskActor->Update();
#endif
      mMaskActor->SetPickable(0);
      mMaskActor->SetVisibility(true);
      ImageActor->SetOpacity(0.75);
      mMaskActor->SetInterpolate(0);
#if VTK_MAJOR_VERSION >= 6 || (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 10)
      mMaskActor->GetMapper()->BorderOn();
#endif
    }

    this->GetRenderer()->AddActor(mMaskActor);

    AdjustResliceToSliceOrientation(mMaskReslice);
    this->UpdateDisplayExtent();
   }
}
void Slicer::RemoveMask()
{
  if (mMaskActor) {
    Renderer->RemoveActor(mMaskActor);
    mMask = NULL;
    mMaskActor = NULL;
    mMaskMapper = NULL;
  }
}
double Slicer::GetMaskOpacity()
{
  return mMaskOpacity;
}
void Slicer::SetMaskOpacity(double opacity)
{
  if (opacity > 1)
  {
    mMaskOpacity = opacity / 10;
  }
  else
  {
    mMaskOpacity = opacity;
  }

  if (mMaskActor)
  {
    mMaskActor->SetOpacity(opacity);
    this->GetRenderWindow()->Render();
  }
}

void Slicer::SetLandmarks(Landmarks* landmarks, int type)
{
  if (landmarks == NULL || type == LANDMARK_TYPE::NONE) {
    return;
  }

  mLandmarks = landmarks;
  mLandmarksType = type;

  double r, g, b;

  r = g = b = 0;
  if (type == LANDMARK_TYPE::DEFAULT)
  {
    r = 0.1; g = 1; b = 0.1;
  }
  else if (type == LANDMARK_TYPE::TUMOR_POINTS)
  {
    r = 1; g = 0.1; b = 1;
  }
  else if (type == LANDMARK_TYPE::TISSUE_POINTS)
  {
    r = 1; g = 1; b = 0.1;
  }

  ///////////////////////////////////////////////////////////////////////////
  if (!mClipBox) 
  {
    mClipBox = vtkSmartPointer<vtkBox>::New();
  }
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  if (!mCross) 
  {
    mCross = vtkSmartPointer<vtkCursor3D>::New();
  }
  mCross->SetFocalPoint(0, 0, 0);
  mCross->SetModelBounds(-3, 3, -3, 3, -3, 3);
  mCross->AllOff();
  mCross->AxesOn();

  if (!mLandGlyph) {
    mLandGlyph = vtkSmartPointer<vtkGlyph3D>::New();
  }
#if VTK_MAJOR_VERSION <= 5
  mLandGlyph->SetSource(mCross->GetOutput());
  mLandGlyph->SetInput(landmarks->mLandDataWithID);
#else
  mLandGlyph->SetSourceConnection(mCross->GetOutputPort());
  mLandGlyph->SetInputData(landmarks->mLandDataWithID);
  mLandGlyph->Update();
#endif
  mLandGlyph->SetScaleModeToDataScalingOff();

  if (!mLandClipper) {
    mLandClipper = vtkSmartPointer<vtkClipPolyData>::New();
  }
  mLandClipper->InsideOutOn();
#if VTK_MAJOR_VERSION <= 5
  mLandClipper->SetInput(mLandGlyph->GetOutput());
#else
  mLandClipper->SetInputConnection(mLandGlyph->GetOutputPort());
  mLandClipper->Update();
#endif
  mLandClipper->SetClipFunction(mClipBox);
  if (!mLandMapper) {
    mLandMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  }
  mLandMapper->SetInputConnection(mLandClipper->GetOutputPort());
  mLandMapper->ScalarVisibilityOff();
  if (!mLandActor) {
    mLandActor = vtkSmartPointer<vtkActor>::New();
  }
  mLandActor->SetMapper(mLandMapper);
  mLandActor->GetProperty()->SetColor(r, g, b);
  mLandActor->GetProperty()->SetOpacity(0.7);
  mLandActor->SetPickable(0);
  mLandActor->SetVisibility(mLandmarksVisibility);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  if (!mLandLabelGlyph) {
    mLandLabelGlyph = vtkSmartPointer<vtkVertexGlyphFilter>::New();
  }
#if VTK_MAJOR_VERSION <= 5
  mLandLabelGlyph->SetInput(landmarks->mLandDataWithID);
#else
  mLandLabelGlyph->SetInputData(landmarks->mLandDataWithID);
#endif
  if (!mLandLabelClipper) {
    mLandLabelClipper = vtkSmartPointer<vtkClipPolyData>::New();
  }
  mLandLabelClipper->InsideOutOn();
#if VTK_MAJOR_VERSION <= 5
  mLandLabelClipper->SetInput(mLandLabelGlyph->GetOutput());
#else
  mLandLabelClipper->SetInputConnection(mLandLabelGlyph->GetOutputPort());
  mLandLabelClipper->Update();
#endif
  mLandLabelClipper->SetClipFunction(mClipBox);
  if (!mLandLabelMapper) 
  {
    mLandLabelMapper = vtkSmartPointer<vtkLabeledDataMapper>::New();
  }
  mLandLabelMapper->SetFieldDataName("index");
#if VTK_MAJOR_VERSION <= 5
  mLandLabelMapper->SetInput(mLandLabelClipper->GetOutput());
#else
  mLandLabelMapper->SetInputConnection(mLandLabelClipper->GetOutputPort());
  mLandLabelMapper->Update();
#endif
  mLandLabelMapper->SetLabelModeToLabelScalars();
  mLandLabelMapper->SetLabelFormat("%.0f");
  mLandLabelMapper->GetLabelTextProperty()->SetColor(r, g, b);
  mLandLabelMapper->GetLabelTextProperty()->SetFontSize(10);
  mLandLabelMapper->GetLabelTextProperty()->SetShadow(0);
  mLandLabelMapper->GetLabelTextProperty()->SetOpacity(0.7);
  {
    vtkSmartPointer<vtkTransform> translate = vtkSmartPointer<vtkTransform>::New();
    switch (SliceOrientation)
    {
    case SLICE_ORIENTATION_XY:
      translate->Translate(3, -3, 0);
      break;
    case SLICE_ORIENTATION_XZ:
      translate->Translate(3, 0, 3);
      break;
    case SLICE_ORIENTATION_YZ:
      translate->Translate(0, 3, 3);
      break;
    default:
      break;
    }
    mLandLabelMapper->SetTransform(translate);
  }
  if (!mLandLabelActor) 
  {
    mLandLabelActor = vtkSmartPointer<vtkActor2D>::New();
  }
  mLandLabelActor->SetMapper(mLandLabelMapper);
  mLandLabelActor->SetVisibility(mLandmarksVisibility);
  ///////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////
  if (!mCircle) {
    mCircle = vtkSmartPointer<vtkRegularPolygonSource>::New();
  }
  mCircle->SetNumberOfSides(50);
  mCircle->SetRadius(1);
  mCircle->SetCenter(0, 0, 0);
  mCircle->SetGeneratePolygon(0);
  switch (SliceOrientation)
  {
  case SLICE_ORIENTATION_XY:
    mCircle->SetNormal(0, 0, 1);
    break;
  case SLICE_ORIENTATION_XZ:
    mCircle->SetNormal(0, 1, 0);
    break;
  case SLICE_ORIENTATION_YZ:
    mCircle->SetNormal(1, 0, 0);
    break;
  default:
    break;
  }
  mCircle->Update();

  if (!mLandRadiusGlyph) {
    mLandRadiusGlyph = vtkSmartPointer<vtkGlyph3D>::New();
  }
#if VTK_MAJOR_VERSION <= 5
  mLandRadiusGlyph->SetSource(mCircle->GetOutput());
  mLandRadiusGlyph->SetInput(landmarks->mLandDataWithRadius);
#else
  mLandRadiusGlyph->SetSourceConnection(mCircle->GetOutputPort());
  mLandRadiusGlyph->SetInputData(landmarks->mLandDataWithRadius);
  mLandRadiusGlyph->Update();
#endif
  mLandRadiusGlyph->SetScaleModeToScaleByScalar();

  if (!mLandRadiusClipper) {
    mLandRadiusClipper = vtkSmartPointer<vtkClipPolyData>::New();
  }
  mLandRadiusClipper->InsideOutOn();
#if VTK_MAJOR_VERSION <= 5
  mLandRadiusClipper->SetInput(mLandRadiusGlyph->GetOutput());
#else
  mLandRadiusClipper->SetInputConnection(mLandRadiusGlyph->GetOutputPort());
  mLandRadiusClipper->Update();
#endif
  mLandRadiusClipper->SetClipFunction(mClipBox);
  if (!mLandRadiusMapper) 
  {
    mLandRadiusMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  }
  mLandRadiusMapper->SetInputConnection(mLandRadiusClipper->GetOutputPort());
  mLandRadiusMapper->ScalarVisibilityOff();
  if (!mLandRadiusActor) 
  {
    mLandRadiusActor = vtkSmartPointer<vtkActor>::New();
  }
  mLandRadiusActor->SetMapper(mLandRadiusMapper);
  mLandRadiusActor->GetProperty()->SetColor(r, g, b);
  //mLandRadiusActor->GetProperty()->SetOpacity(0.7);
  mLandRadiusActor->SetPickable(0);
  mLandRadiusActor->SetVisibility(mLandmarksVisibility);

  ///////////////////////////////////////////////////////////////////////////
  this->UpdateDisplayExtent();
  this->GetRenderer()->AddActor(mLandActor);
  this->GetRenderer()->AddActor(mLandLabelActor);
  this->GetRenderer()->AddActor(mLandRadiusActor);
  ///////////////////////////////////////////////////////////////////////////
}

void Slicer::SetSliceOrientation(int orientation)
{
  int extent[6];
#if VTK_MAJOR_VERSION <= 5
  this->GetInput()->GetWholeExtent(extent);
#else
  this->GetInput()->GetExtent(extent);
#endif
  if (extent[5] - extent[4] <= 2) {
    orientation = vtkImageViewer2::SLICE_ORIENTATION_XY;
  }

  if (orientation < vtkImageViewer2::SLICE_ORIENTATION_YZ ||
    orientation > vtkImageViewer2::SLICE_ORIENTATION_XY) {
    vtkErrorMacro("Error - invalid slice orientation " << orientation);
    return;
  }

  this->SliceOrientation = orientation;

  int *range = this->GetSliceRange();
  if (range) {
    this->Slice = static_cast<int>((range[0] + range[1]) * 0.5);
  }

  this->UpdateOrientation();
  this->UpdateDisplayExtent();

  if (this->Renderer && this->GetInput()) {
    double scale = this->Renderer->GetActiveCamera()->GetParallelScale();
    this->Renderer->ResetCamera();
    this->Renderer->GetActiveCamera()->SetParallelScale(scale);
  }
}
int Slicer::GetOrientation()
{
  return this->SliceOrientation;
}

void Slicer::UpdateDisplayExtent()
{
  vtkSmartPointer< vtkImageData > input = this->GetInput();
  if (!input || !this->ImageActor) {
    return;
  }
#if VTK_MAJOR_VERSION <= 5
  input->UpdateInformation();
#endif

  int w_ext[6];
#if VTK_MAJOR_VERSION <= 5
  int* ext = GetInput()->GetWholeExtent();
#else
  int* ext = GetInput()->GetExtent();
#endif
  for (int i = 0; i < 6; i++) {
    w_ext[i] = ext[i];
  }

  int s = this->Slice > ext[this->SliceOrientation * 2 + 1] ? ext[this->SliceOrientation * 2 + 1] : this->Slice;
  w_ext[this->SliceOrientation * 2] = s;
  w_ext[this->SliceOrientation * 2 + 1] = s;

  this->ImageActor->SetDisplayExtent(w_ext);

  if (mOverlay && mOverlayActor->GetVisibility()) {
    AdjustResliceToSliceOrientation(mOverlayReslice);
    int overExtent[6];
    this->ConvertImageToImageDisplayExtent(input, w_ext, mOverlayReslice->GetOutput(), overExtent);
#if VTK_MAJOR_VERSION <= 5
    ClipDisplayedExtent(overExtent, mOverlayMapper->GetInput()->GetWholeExtent());
#else
    ClipDisplayedExtent(overExtent, mOverlayMapper->GetImageDataInput(0)->GetExtent());
#endif
    mOverlayActor->SetDisplayExtent(overExtent);
  }

  if (mMask && mMaskActor->GetVisibility()) {
    AdjustResliceToSliceOrientation(mMaskReslice);
    int overExtent[6];
    this->ConvertImageToImageDisplayExtent(input, w_ext, mMaskReslice->GetOutput(), overExtent);
#if VTK_MAJOR_VERSION <= 5
    ClipDisplayedExtent(overExtent, mMaskMapper->GetInput()->GetWholeExtent());
#else
    ClipDisplayedExtent(overExtent, mMaskMapper->GetImageDataInput(0)->GetExtent());
#endif
    mMaskActor->SetDisplayExtent(overExtent);
  }

  double* camera = Renderer->GetActiveCamera()->GetPosition();
  double* image_bounds = ImageActor->GetBounds();
  double position[3] = { 0, 0, 0 };
  position[this->SliceOrientation] = image_bounds[this->SliceOrientation * 2];

  double offset = 1;
  if (camera[this->SliceOrientation] < image_bounds[this->SliceOrientation * 2]) {
    offset = -1;
  }

  if (mLandActor) {
    if (mClipBox) {
      double bounds[6];
      for (unsigned int i = 0; i < 6; i++) {
        bounds[i] = ImageActor->GetBounds()[i];
      }
      bounds[this->SliceOrientation * 2] = ImageActor->GetBounds()[this->SliceOrientation * 2] - fabs(this->GetInput()->GetSpacing()[this->SliceOrientation]) / 2;
      bounds[this->SliceOrientation * 2 + 1] = ImageActor->GetBounds()[this->SliceOrientation * 2 + 1] + fabs(this->GetInput()->GetSpacing()[this->SliceOrientation]) / 2;
      mClipBox->SetBounds(bounds);
      UpdateLandmarks();
    }

    position[this->SliceOrientation] = offset;
    mLandActor->SetPosition(position);
  }

  if (this->Renderer) {
    if (this->InteractorStyle &&
      this->InteractorStyle->GetAutoAdjustCameraClippingRange()) {
      this->Renderer->ResetCameraClippingRange();
    }
    else {
      vtkCamera *cam = this->Renderer->GetActiveCamera();
      if (cam) {
        double bounds[6];
        this->ImageActor->GetBounds(bounds);
        double spos = (double)bounds[this->SliceOrientation * 2];
        double cpos = (double)cam->GetPosition()[this->SliceOrientation];
        double range = fabs(spos - cpos);
        double *spacing = input->GetSpacing();
        double avg_spacing = ((double)spacing[0] + (double)spacing[1] + (double)spacing[2]) / 3.0;
        cam->SetClippingRange(range - avg_spacing * 3.0, range + avg_spacing * 3.0);
      }
    }
  }
}

void Slicer::UpdateOrientation()
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (cam)
  {
    switch (this->SliceOrientation)
    {
    case vtkImageViewer2::SLICE_ORIENTATION_XY:
		this->mCornerAnnotation->SetText(0, "Axial");
      cam->SetFocalPoint(0, 0, 0);
      cam->SetPosition(0, 0, -1);
      cam->SetViewUp(0, -1, 0);
      break;
    case vtkImageViewer2::SLICE_ORIENTATION_XZ:
		this->mCornerAnnotation->SetText(0, "Coronal");
      cam->SetFocalPoint(0, 0, 0);
      cam->SetPosition(0, -1, 0);
      cam->SetViewUp(0, 0, 1);
      break;
    case vtkImageViewer2::SLICE_ORIENTATION_YZ:
		this->mCornerAnnotation->SetText(0, "Sagittal");
      cam->SetFocalPoint(0, 0, 0);
      cam->SetPosition(1, 0, 0);
      cam->SetViewUp(0, 0, 1);
      break;
    }
  }
}

void Slicer::SetOpacity(double s)
{
  this->GetImageActor()->SetOpacity(s);
}

void Slicer::SetRenderWindow(int orientation, vtkRenderWindow * rw)
{
  if (rw == nullptr)
    this->Superclass::SetRenderWindow(rw);
  else
  {
    this->Superclass::SetRenderWindow(rw);
    this->SetupInteractor(rw->GetInteractor());
    this->GetRenderer()->AddActor(pdmA);
    this->GetRenderer()->ResetCamera();
    SetSliceOrientation(2 - (orientation % 3));
    ResetCamera();
    vtkCamera *camera = this->GetRenderer()->GetActiveCamera();
    camera->SetParallelScale(camera->GetParallelScale() * 0.8);
  }
}

void Slicer::ResetCamera()
{
  this->GetRenderer()->ResetCamera();
}

void Slicer::SetDisplayMode(bool i)
{
  this->GetRenderer()->SetDraw(i);
  if (i) {
    UpdateDisplayExtent();
  }
}

void Slicer::FlipHorizontalView()
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (cam) {
    double *position = cam->GetPosition();
    double factor[3] = { 1, 1, 1 };
    factor[this->SliceOrientation] = -1;
    cam->SetPosition(factor[0] * position[0], factor[1] * position[1], factor[2] * position[2]);

    this->Renderer->ResetCameraClippingRange();
    this->UpdateDisplayExtent();
  }
}

void Slicer::FlipVerticalView()
{
  vtkCamera *cam = this->Renderer ? this->Renderer->GetActiveCamera() : NULL;
  if (cam) {
    FlipHorizontalView();
    double *viewup = cam->GetViewUp();
    cam->SetViewUp(-viewup[0], -viewup[1], -viewup[2]);
    this->UpdateDisplayExtent();
  }
}

void Slicer::SetColorWindow(double window)
{
  vtkLookupTable* LUT = static_cast<vtkLookupTable*>(this->GetWindowLevel()->GetLookupTable());
  if (LUT) {
    double level = this->GetWindowLevel()->GetLevel();
    LUT->SetTableRange(level - fabs(window) / 2, level + fabs(window) / 2);
    LUT->Build();
  }
  this->vtkImageViewer2::SetColorWindow(window);
}

void Slicer::SetColorLevel(double level)
{
  vtkLookupTable* LUT = static_cast<vtkLookupTable*>(this->GetWindowLevel()->GetLookupTable());
  if (LUT) {
    double window = this->GetWindowLevel()->GetWindow();
    LUT->SetTableRange(level - fabs(window) / 2, level + fabs(window) / 2);
    LUT->Build();
  }
  this->vtkImageViewer2::SetColorLevel(level);
}

double Slicer::GetScalarComponentAsDouble(vtkSmartPointer< vtkImageData > image, double X, double Y, double Z, int &ix, int &iy, int &iz, int component)
{
  ix = ROUND(X);
  iy = ROUND(Y);
  iz = ROUND(Z);
#if VTK_MAJOR_VERSION <= 5
  if (ix < image->GetWholeExtent()[0] || ix > image->GetWholeExtent()[1] ||
    iy < image->GetWholeExtent()[2] || iy > image->GetWholeExtent()[3] ||
    iz < image->GetWholeExtent()[4] || iz > image->GetWholeExtent()[5]) {
#else
  if (ix < image->GetExtent()[0] || ix > image->GetExtent()[1] ||
    iy < image->GetExtent()[2] || iy > image->GetExtent()[3] ||
    iz < image->GetExtent()[4] || iz > image->GetExtent()[5]) {
#endif
    return std::numeric_limits<double>::quiet_NaN();
  }

#if VTK_MAJOR_VERSION <= 5
  image->SetUpdateExtent(ix, ix, iy, iy, iz, iz);
  image->Update();
#else
  int extent[6] = { ix, ix, iy, iy, iz, iz };
  //vtkStreamingDemandDrivenPipeline::SetUpdateExtent(image->GetInformation(), extent);
  vtkStreamingDemandDrivenPipeline::SetWholeExtent(image->GetInformation(), extent);
  //image->SetExtent(ix, ix, iy, iy, iz, iz);
#endif
  return image->GetScalarComponentAsDouble(ix, iy, iz, component);
  }

void Slicer::Render()
{
	if (this->GetInput() == nullptr)
		return;

  if (pdmA->GetVisibility()) {
    double bX0 = this->GetImageActor()->GetDisplayExtent()[0];
    double bX1 = this->GetImageActor()->GetDisplayExtent()[1];
    double bY0 = this->GetImageActor()->GetDisplayExtent()[2];
    double bY1 = this->GetImageActor()->GetDisplayExtent()[3];
    double bZ0 = this->GetImageActor()->GetDisplayExtent()[4];
    double bZ1 = this->GetImageActor()->GetDisplayExtent()[5];
    double bx0 = bX0 * this->GetInput()->GetSpacing()[0] + this->GetInput()->GetOrigin()[0];
    double bx1 = bX1 * this->GetInput()->GetSpacing()[0] + this->GetInput()->GetOrigin()[0];
    double by0 = bY0 * this->GetInput()->GetSpacing()[1] + this->GetInput()->GetOrigin()[1];
    double by1 = bY1 * this->GetInput()->GetSpacing()[1] + this->GetInput()->GetOrigin()[1];
    double bz0 = bZ0 * this->GetInput()->GetSpacing()[2] + this->GetInput()->GetOrigin()[2];
    double bz1 = bZ1 * this->GetInput()->GetSpacing()[2] + this->GetInput()->GetOrigin()[2];
    double bounds[6];
    double x = mCursor[0];
    double y = mCursor[1];
    double z = mCursor[2];
    double xCursor = (x - this->GetInput()->GetOrigin()[0]) / this->GetInput()->GetSpacing()[0];
    double yCursor = (y - this->GetInput()->GetOrigin()[1]) / this->GetInput()->GetSpacing()[1];
    double zCursor = (z - this->GetInput()->GetOrigin()[2]) / this->GetInput()->GetSpacing()[2];
    //
    // round up pixel values
    xCursor = ROUND(xCursor);
    yCursor = ROUND(yCursor);
    zCursor = ROUND(zCursor);
    x = xCursor * this->GetInput()->GetSpacing()[0] + this->GetInput()->GetOrigin()[0];
    y = yCursor * this->GetInput()->GetSpacing()[1] + this->GetInput()->GetOrigin()[1];
    z = zCursor * this->GetInput()->GetSpacing()[2] + this->GetInput()->GetOrigin()[2];
    //
    if (xCursor >= bX0 && xCursor <= bX1 &&
      yCursor >= bY0 && yCursor <= bY1 &&
      zCursor >= bZ0 && zCursor <= bZ1)
    {
      vtkRenderer * renderer = this->Renderer;

      renderer->WorldToView(x, y, z);
      renderer->ViewToNormalizedViewport(x, y, z);
      renderer->NormalizedViewportToViewport(x, y);
      renderer->ViewportToNormalizedDisplay(x, y);
      renderer->NormalizedDisplayToDisplay(x, y);

      this->GetRenderer()->WorldToView(bx0, by0, bz0);
      this->GetRenderer()->ViewToNormalizedViewport(bx0, by0, bz0);
      this->GetRenderer()->NormalizedViewportToViewport(bx0, by0);
      this->GetRenderer()->ViewportToNormalizedDisplay(bx0, by0);
      this->GetRenderer()->NormalizedDisplayToDisplay(bx0, by0);
      //
      this->GetRenderer()->WorldToView(bx1, by1, bz1);
      this->GetRenderer()->ViewToNormalizedViewport(bx1, by1, bz1);
      this->GetRenderer()->NormalizedViewportToViewport(bx1, by1);
      this->GetRenderer()->ViewportToNormalizedDisplay(bx1, by1);
      this->GetRenderer()->NormalizedDisplayToDisplay(bx1, by1);

      bounds[0] = __min(bx0, bx1);
      bounds[1] = __max(bx0, bx1);
      bounds[2] = __min(by0, by1);
      bounds[3] = __max(by0, by1);
      bounds[4] = __min(bz0, bz1);
      bounds[5] = __max(bz0, bz1);
      crossCursor->SetModelBounds(bounds);
      crossCursor->SetOutline(1);

      crossCursor->SetFocalPoint(x, y, z);
    }
    else {
      crossCursor->SetFocalPoint(-1, -1, z);
    }
  }

  if (mOverlay && mOverlayActor->GetVisibility()) {
#if VTK_MAJOR_VERSION <= 5
    mOverlayMapper->GetOutput()->SetUpdateExtent(mOverlayActor->GetDisplayExtent());
    mOverlayMapper->GetOutput()->Update();
#else
    //vtkStreamingDemandDrivenPipeline::SetUpdateExtent(mOverlayMapper->GetInformation(), mOverlayActor->GetDisplayExtent());
    vtkStreamingDemandDrivenPipeline::SetWholeExtent(mOverlayMapper->GetInformation(), mOverlayActor->GetDisplayExtent());
#endif
    mOverlayMapper->Update();
  }

  if (mMask && mMaskActor->GetVisibility()) {
#if VTK_MAJOR_VERSION <= 5
    mMaskMapper->GetOutput()->SetUpdateExtent(mMaskActor->GetDisplayExtent());
    mMaskMapper->GetOutput()->Update();
#else
    //vtkStreamingDemandDrivenPipeline::SetUpdateExtent(mMaskMapper->GetInformation(), mMaskActor->GetDisplayExtent());
    vtkStreamingDemandDrivenPipeline::SetWholeExtent(mMaskMapper->GetInformation(), mMaskActor->GetDisplayExtent());
#endif
    mMaskMapper->Update();
  }

  if (mLandMapper) {
    UpdateLandmarks();
  }

  this->GetRenderWindow()->Render();
}

void Slicer::UpdateCursorPosition()
{

}

void Slicer::UpdateLandmarks()
{
  if (!mLandClipper) {
    return;
  }

  vtkPolyData *pd = static_cast<vtkPolyData*>(mLandClipper->GetInput());
  if (pd->GetPoints())
  {
    mLandGlyph->Modified();
    mLandGlyph->Update();
    mLandClipper->Update();
    mLandMapper->Update();
  }
}

void Slicer::SetSlice(int slice)
{
  int *range = this->GetSliceRange();
  if (range) {
    if (slice < range[0]) {
      slice = range[0];
    }
    else if (slice > range[1]) {
      slice = range[1];
    }
  }

  if (this->Slice == slice) {
    return;
  }

  this->Slice = slice;
  this->Modified();
  this->UpdateDisplayExtent();
}

void Slicer::ForceUpdateDisplayExtent()
{
  this->UpdateDisplayExtent();
}

int* Slicer::GetDisplayExtent()
{
  return this->GetImageActor()->GetDisplayExtent();
}

void Slicer::AdjustResliceToSliceOrientation(vtkImageReslice *reslice)
{
  reslice->SetOutputOriginToDefault();
  reslice->SetOutputSpacingToDefault();
#if VTK_MAJOR_VERSION <= 5
  reslice->GetOutput()->UpdateInformation();
#else
  //reslice->UpdateInformation();
  reslice->Update();
#endif

  double origin[3];
  double spacing[3];
  reslice->GetOutput()->GetOrigin(origin);
  reslice->GetOutput()->GetSpacing(spacing);

  spacing[this->SliceOrientation] = mImageReslice->GetOutput()->GetSpacing()[this->SliceOrientation];

  origin[this->SliceOrientation] -= mImageReslice->GetOutput()->GetOrigin()[this->SliceOrientation];
  origin[this->SliceOrientation] /= mImageReslice->GetOutput()->GetSpacing()[this->SliceOrientation];

  origin[this->SliceOrientation] = ceil(origin[this->SliceOrientation]);

  origin[this->SliceOrientation] *= mImageReslice->GetOutput()->GetSpacing()[this->SliceOrientation];
  origin[this->SliceOrientation] += mImageReslice->GetOutput()->GetOrigin()[this->SliceOrientation];

  reslice->SetOutputOrigin(origin);
  reslice->SetOutputSpacing(spacing);
#if VTK_MAJOR_VERSION <= 5
  reslice->UpdateInformation();
  reslice->GetOutput()->UpdateInformation();
#else
  //reslice->UpdateInformation();
  reslice->Update();
#endif
}

void Slicer::ConvertImageToImageDisplayExtent(vtkSmartPointer< vtkImageData > sourceImage, const int sourceExtent[6], vtkSmartPointer< vtkImageData > targetImage, int targetExtent[6])
{
  double dExtents[6];
  for (unsigned int i = 0; i < 6; i++) {
    dExtents[i] = sourceImage->GetOrigin()[i / 2] + sourceImage->GetSpacing()[i / 2] * sourceExtent[i];

    dExtents[i] = (dExtents[i] - targetImage->GetOrigin()[i / 2]) / targetImage->GetSpacing()[i / 2];

    targetExtent[i] = (int)floor(dExtents[i]);
  }
}

void Slicer::ClipDisplayedExtent(int extent[6], int refExtent[6])
{
  bool out = false;
  int maxBound = 6;

  if (refExtent[4] == refExtent[5]) {
    maxBound = 4;
    extent[4] = refExtent[4];
    extent[5] = refExtent[5];
  }

  for (int i = 0; i < maxBound; i = i + 2) {
    if (extent[i] > refExtent[i + 1] || extent[i + 1] < refExtent[i]) {
      out = true;
      break;
    }
    extent[i] = (extent[i] > refExtent[i]) ? extent[i] : refExtent[i];
    extent[i] = (extent[i] < refExtent[i + 1]) ? extent[i] : refExtent[i + 1];
    extent[i + 1] = (extent[i + 1] > refExtent[i]) ? extent[i + 1] : refExtent[i];
    extent[i + 1] = (extent[i + 1] < refExtent[i + 1]) ? extent[i + 1] : refExtent[i + 1];
  }
  if (out) {
    for (int i = 0; i < maxBound; i = i + 2) {
      extent[i] = refExtent[i];
      extent[i + 1] = refExtent[i];
    }
  }
}
void Slicer::SetImageSeriesDescription(std::string description)
{
	this->mCornerAnnotation->SetText(vtkCornerAnnotation::UpperLeft, description.c_str());
}
void Slicer::ResetMap()
{
  mMask->DeepCopy(mMaskOriginal);
  mOriginalMaskSaved = false;
}

void Slicer::AddLabelToMap(int label)
{
  if (!mOriginalMaskSaved)
  {
    mMaskOriginal = vtkSmartPointer<vtkImageData>::New();
    mMaskOriginal->DeepCopy(mMask);
    mOriginalMaskSaved = true;
  }

  auto dims = mMaskOriginal->GetDimensions();
  for (int z = 0; z < dims[2]; z++)
  {
    for (int y = 0; y < dims[1]; y++)
    {
      for (int x = 0; x < dims[0]; x++)
      {
        auto source = static_cast< float* >(mMaskOriginal->GetScalarPointer(x, y, z));
        auto destination = static_cast< float* >(mMask->GetScalarPointer(x, y, z));

        if (static_cast< int >(*source) == label)
        {
          *destination = label;
        }
      }
    }
  }

}

void Slicer::ShowLabelOnMap(int label)
{
  if (!mOriginalMaskSaved)
  {
    mMaskOriginal = vtkSmartPointer<vtkImageData>::New();
    mMaskOriginal->DeepCopy(mMask);
    mOriginalMaskSaved = true;
  }

  auto dims = mMaskOriginal->GetDimensions();
  for (int z = 0; z < dims[2]; z++)
  {
    for (int y = 0; y < dims[1]; y++)
    {
      for (int x = 0; x < dims[0]; x++)
      {
        auto source = static_cast<float*>(mMaskOriginal->GetScalarPointer(x, y, z));
        auto destination = static_cast<float*>(mMask->GetScalarPointer(x, y, z));

        if (static_cast<int>(*source) == label)
        {
          *destination = label;
        }
        else
        {
          *destination = 0;
        }
      }
    }
  }
}


void Slicer::SetOverlayOpacity(double opacity)
{
  mOverlayOpacity = opacity;

  if (mOverlayActor)
  {
    mOverlayActor->SetOpacity(mOverlayOpacity);
    this->GetRenderWindow()->Render();
  }
}
void Slicer::SetOverlay(vtkImageData* overlay)
{
  if (overlay != NULL) {
    mOverlay = overlay;

    if (!mOverlayReslice) {
      mOverlayReslice = vtkSmartPointer<vtkImageReslice>::New();
      //mOverlayReslice->SetInterpolationModeToLinear();
      mOverlayReslice->SetInterpolationModeToNearestNeighbor();
      mOverlayReslice->AutoCropOutputOn();
      mOverlayReslice->SetBackgroundColor(0, 0, 0, 1);
    }
    mOverlayReslice->SetResliceTransform(mTransform);
#if VTK_MAJOR_VERSION <= 5
    mOverlayReslice->SetInput(0, overlay);
    //mOverlayReslice->UpdateInformation();
#else
    mOverlayReslice->SetInputData(0, overlay);
    mOverlayReslice->Update();
#endif

    if (!mOverlayMapper) {
      mOverlayMapper = vtkSmartPointer<vtkImageMapToWindowLevelColors>::New();
    }
#if VTK_MAJOR_VERSION <= 5
    mOverlayMapper->SetInput(mOverlayReslice->GetOutput());
#else
    mOverlayMapper->SetInputConnection(mOverlayReslice->GetOutputPort());
    mOverlayMapper->Update();
#endif

    if (!mOverlayActor) {
      mOverlayActor = vtkSmartPointer<vtkImageActor>::New();
#if VTK_MAJOR_VERSION <= 5
      mOverlayActor->SetInput(mOverlayMapper->GetOutput());
#else
      mOverlayActor->GetMapper()->SetInputConnection(mOverlayMapper->GetOutputPort());
      mOverlayActor->Update();
#endif
      mOverlayActor->SetPickable(0);
      mOverlayActor->SetVisibility(true);
      mOverlayActor->SetOpacity(mOverlayOpacity);
      //
      mOverlayActor->SetInterpolate(0);
      //
#if VTK_MAJOR_VERSION >= 6 || (VTK_MAJOR_VERSION >= 5 && VTK_MINOR_VERSION >= 10)
      mOverlayActor->GetMapper()->BorderOn();
#endif
    }

    this->GetRenderer()->AddActor(mOverlayActor);

    AdjustResliceToSliceOrientation(mOverlayReslice);
    this->UpdateDisplayExtent();
  }
}
void Slicer::RemoveOverlay()
{
  if (mOverlayActor) {
    Renderer->RemoveActor(mOverlayActor);
    mOverlay = NULL;
    mOverlayActor = NULL;
    mOverlayMapper = NULL;
  }
}

void Slicer::SetComparisonMode(bool mode)
{
  this->m_ComparisonMode = mode;
}

bool Slicer::GetComparisonMode()
{
  return this->m_ComparisonMode;
}
