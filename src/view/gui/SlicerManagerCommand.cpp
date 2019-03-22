///////////////////////////////////////////////////////////////////////////////////////
// SlicerManagerCommand.cxx
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#include "SlicerManagerCommand.h"
#include "Slicer.h"
#include "fMainWindow.h"
#include "InteractorStyleNavigator.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include "vtkCamera.h"

SlicerManagerCommand::SlicerManagerCommand()
{
  mStartSlicer = -1;
  mSlicerNumber = -1;
  mw = NULL;
}

int SlicerManagerCommand::FindSlicerNumber(vtkRenderWindow* renwin)
{
  int rvalue;
  if (renwin != SM->GetSlicer(mSlicerNumber)->GetRenderWindow() || !SM->GetSlicer(mSlicerNumber)->GetRenderer()->GetDraw())
  {
    rvalue = -1;
  }
  else
  {
    rvalue = mSlicerNumber;
  }
  return rvalue;
}

void SlicerManagerCommand::Execute(vtkObject *caller, unsigned long event, void *vtkNotUsed(callData))
{
  InteractorStyleNavigator *isi = dynamic_cast<InteractorStyleNavigator *>(caller);
  if (isi)
  {
    int VisibleInWindow = this->FindSlicerNumber(isi->GetInteractor()->GetRenderWindow());
    fMainWindow* w = qobject_cast<fMainWindow*>(qApp->activeWindow());
    if (w != NULL)
    {
      mw = w;
    }
    if (mw == NULL)
      return;

    double x = isi->GetInteractor()->GetEventPosition()[0];
    double y = isi->GetInteractor()->GetEventPosition()[1];

    bool bCtrlKey = isi->GetInteractor()->GetControlKey();
    bool bShiftKey = isi->GetInteractor()->GetShiftKey();
    bool bAltKey = isi->GetInteractor()->GetAltKey();
    std::string KeyPress;

    if (event == vtkCommand::KeyPressEvent) 
    {
      KeyPress = isi->GetInteractor()->GetKeySym();
    }

    if (event == vtkCommand::StartPickEvent && VisibleInWindow < 0)
    {
      for (int i = 0; i < this->SM->GetNumberOfSlicers(); i++) 
      {
        if (this->SM->GetSlicer(i)->GetCursorVisibility() && !this->SM->IsLinked()) 
        {
          this->SM->GetSlicer(i)->SetCursorVisibility(0);
          this->SM->GetSlicer(i)->Render();
        }
      }
    }
    if (VisibleInWindow >= 0)
    {
      // vtkRenderer* renderer = NULL;
      // renderer = this->SM->GetSlicer(VisibleInWindow)->GetRenderer();
      //int order = SM->GetOrder();
      std::string tmp = SM->GetId();

      if (event == vtkCommand::KeyPressEvent)
      {
        if (bShiftKey || bAltKey)
        {
          if (KeyPress == "0")
          {
            this->SM->GetSlicer(0)->ResetMap();
          }
          if (KeyPress == "1")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(1);
          }
          if (KeyPress == "2")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(2);
          }
          if (KeyPress == "3")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(3);
          }
          if (KeyPress == "4")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(4);
          }
          if (KeyPress == "5")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(5);
          }
          if (KeyPress == "6")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(6);
          }
          if (KeyPress == "7")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(7);
          }
          if (KeyPress == "8")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(8);
          }
          if (KeyPress == "9")
          {
            this->SM->GetSlicer(0)->AddLabelToMap(9);
          }
          this->SM->GetSlicer(0)->mMask->Modified();
          //this->SM->GetSlicer(0)->Render();
          this->SM->Render();
        }
        if (bCtrlKey)
        {
          if (KeyPress == "0")
          {
            this->SM->GetSlicer(0)->ResetMap();
          }
          if (KeyPress == "1")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(1);
          }
          if (KeyPress == "2")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(2);
          }
          if (KeyPress == "3")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(3);
          }
          if (KeyPress == "4")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(4);
          }
          if (KeyPress == "5")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(5);
          }
          if (KeyPress == "6")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(6);
          }
          if (KeyPress == "7")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(7);
          }
          if (KeyPress == "8")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(8);
          }
          if (KeyPress == "9")
          {
            this->SM->GetSlicer(0)->ShowLabelOnMap(9);
          }
          this->SM->GetSlicer(0)->mMask->Modified();
          //this->SM->GetSlicer(0)->Render();
          this->SM->Render();
        }
        if (KeyPress == "l") 
        {
          this->SM->ToggleInterpolation();
          this->SM->Render();
          return;
        }
        if (KeyPress == "h") 
        {
          for (int i = 0; i < this->SM->GetNumberOfSlicers(); i++)
          {
            int s;
            s = this->SM->GetSlicer(i)->GetCursorVisibility();
            this->SM->GetSlicer(i)->SetCursorVisibility(1 - s);
            //this->SM->GetSlicer(i)->SetLandmarksVisibility(1 - s);
          }
          this->SM->Render();
          return;
        }
        if (KeyPress == "m")
        {
          mw->SetOpacity();
        }
        if (KeyPress == "r") 
        {
          this->SM->GetSlicer(VisibleInWindow)->ResetCamera();
          //
          // adjust scale
          vtkCamera *camera = this->SM->GetSlicer(VisibleInWindow)->GetRenderer()->GetActiveCamera();
          camera->SetParallelScale(camera->GetParallelScale() * 0.8);
          //
          this->SM->GetSlicer(VisibleInWindow)->Render();
          mw->UpdateLinkedNavigation(this->SM->GetSlicer(VisibleInWindow));
          this->SM->Render();
          return;
        }
        if (KeyPress == "a")
        {
          this->SM->SetPreset(PRESET_AUTO);
          this->SM->UpdateWindowLevel();
          this->SM->Render();
          return;
        }
        if (mw->m_drawShapeMode == SHAPE_MODE_NONE)
        {
          if (KeyPress == "1")
          {
            this->SM->NextImageWithOrder(0);
            return;
          }
          if (KeyPress == "2")
          {
            this->SM->NextImageWithOrder(1);
            return;
          }
          if (KeyPress == "3")
          {
            this->SM->NextImageWithOrder(2);
            return;
          }
          if (KeyPress == "4")
          {
            this->SM->NextImageWithOrder(3);
            return;
          }
          if (KeyPress == "5")
          {
            this->SM->NextImageWithOrder(4);
            return;
          }
          if (KeyPress == "6")
          {
            this->SM->NextImageWithOrder(5);
            return;
          }
          if (KeyPress == "7")
          {
            this->SM->NextImageWithOrder(6);
            return;
          }
          if (KeyPress == "8")
          {
            this->SM->NextImageWithOrder(7);
            return;
          }
          if (KeyPress == "9")
          {
            this->SM->NextImageWithOrder(8);
            return;
          }
          if (KeyPress == "space")
          {
            double* current = this->SM->GetSlicer(VisibleInWindow)->GetCurrentPosition();
            if (bShiftKey)
            {
              this->SM->AddLandmarkShift(current[0], current[1], current[2]);
            }
            else if (bCtrlKey)
            {
              this->SM->AddLandmarkRadius(current[0], current[1], current[2]);
            }
            else
            {
              this->SM->AddLandmark(current[0], current[1], current[2]);
            }
            this->SM->GetSlicer(VisibleInWindow)->UpdateLandmarks();
            this->SM->Render();
            return;
          }
        }
      }

      if (event == vtkCommand::EndPickEvent)
      {
        if (!m_shapeBuffer.empty())
        {
          this->SM->ActionAdded(m_shapeBuffer);
          m_shapeBuffer.clear();
        }
        if (!m_undoBuffer.empty())
        {
          this->SM->ActionAdded(m_undoBuffer);
          m_undoBuffer.clear();
        }
        if (VisibleInWindow >= 0)
        {
          this->SM->LeftButtonReleaseEvent(VisibleInWindow);
        }

        return;
      }

      if (event == vtkCommand::StartWindowLevelEvent)
      {
        mStartSlicer = -1;
        this->InitialWindow = this->SM->GetColorWindow();
        this->InitialLevel = this->SM->GetColorLevel();

        if (VisibleInWindow >= 0)
        {
          mStartSlicer = VisibleInWindow;
        }
        return;
      }

      if (event == vtkCommand::EndWindowLevelEvent)
      {
        mStartSlicer = -1;
      }

      if (event == vtkCommand::EndInteractionEvent)
      {
        this->SM->Picked();
        mw->UpdateLinkedNavigation(this->SM->GetSlicer(VisibleInWindow));
        this->SM->Render();
        return;
      }
    }
    if (VisibleInWindow >= 0)
    {
      double xWorld = 0;
      double yWorld = 0;
      double zWorld = 0;

      this->SM->Activated();

      if (event == vtkCommand::LeftButtonReleaseEvent)
      {
        isi->OnLeftButtonUp();
        // Do this only on double click and when active Tab is tumor panel and control is not pressed
        if (isi->isDoubleClick() && !mw->tumorPointSelected() && (mw->getActiveTabId() == 1))
        {        
          double* current = this->SM->GetSlicer(VisibleInWindow)->GetCurrentPosition();
          this->SM->AddLandmarkShift(current[0], current[1], current[2]);
          this->SM->GetSlicer(VisibleInWindow)->UpdateLandmarks();
          this->SM->Render();
        }

      }
      if (event == vtkCommand::MouseWheelForwardEvent && bCtrlKey) {
        double factor = 2;
        this->Dolly(pow((double)1.1, factor), isi->GetInteractor());
        Execute(caller, vtkCommand::EndInteractionEvent, NULL);
      }
      else if (event == vtkCommand::MouseWheelBackwardEvent && bCtrlKey) {
        double factor = -2;
        this->Dolly(pow((double)1.1, factor), isi->GetInteractor());
        Execute(caller, vtkCommand::EndInteractionEvent, NULL);
      }

      if ((event == vtkCommand::MouseWheelForwardEvent && !bCtrlKey) || (event == vtkCommand::KeyPressEvent && KeyPress == "Up")) {
        double* current = this->SM->GetSlicer(VisibleInWindow)->GetCurrentPosition();
        xWorld = current[0];
        yWorld = current[1];
        zWorld = current[2];
        switch (this->SM->GetSlicer(VisibleInWindow)->GetSliceOrientation()) {
        case vtkImageViewer2::SLICE_ORIENTATION_XY:
          zWorld = (this->SM->GetSlicer(VisibleInWindow)->GetSlice() + 1)*this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing()[2] + this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin()[2];
          break;
        case vtkImageViewer2::SLICE_ORIENTATION_XZ:
          yWorld = (this->SM->GetSlicer(VisibleInWindow)->GetSlice() + 1)*this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing()[1] + this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin()[1];
          break;
        case vtkImageViewer2::SLICE_ORIENTATION_YZ:
          xWorld = (this->SM->GetSlicer(VisibleInWindow)->GetSlice() + 1)*this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing()[0] + this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin()[0];
          break;
        }
        this->SM->GetSlicer(VisibleInWindow)->SetCurrentPosition(xWorld, yWorld, zWorld);
        this->SM->GetSlicer(VisibleInWindow)->UpdateCursorPosition();
        //
        this->SM->Picked();
        this->SM->UpdateViews(VisibleInWindow);
        this->SM->UpdateLinked(VisibleInWindow);
        this->SM->UpdateInfoOnCursorPosition(VisibleInWindow);
      }
      else if ((event == vtkCommand::MouseWheelBackwardEvent && !bCtrlKey) || (event == vtkCommand::KeyPressEvent && KeyPress == "Down")) {
        double* current = this->SM->GetSlicer(VisibleInWindow)->GetCurrentPosition();
        xWorld = current[0];
        yWorld = current[1];
        zWorld = current[2];
        switch (this->SM->GetSlicer(VisibleInWindow)->GetSliceOrientation()) {
        case vtkImageViewer2::SLICE_ORIENTATION_XY:
          zWorld = (this->SM->GetSlicer(VisibleInWindow)->GetSlice() - 1)*this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing()[2] + this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin()[2];
          break;
        case vtkImageViewer2::SLICE_ORIENTATION_XZ:
          yWorld = (this->SM->GetSlicer(VisibleInWindow)->GetSlice() - 1)*this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing()[1] + this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin()[1];
          break;
        case vtkImageViewer2::SLICE_ORIENTATION_YZ:
          xWorld = (this->SM->GetSlicer(VisibleInWindow)->GetSlice() - 1)*this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing()[0] + this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin()[0];
          break;
        }
        this->SM->GetSlicer(VisibleInWindow)->SetCurrentPosition(xWorld, yWorld, zWorld);
        this->SM->GetSlicer(VisibleInWindow)->UpdateCursorPosition();
        //
        this->SM->Picked();
        this->SM->UpdateViews(VisibleInWindow);
        this->SM->UpdateLinked(VisibleInWindow);
        this->SM->UpdateInfoOnCursorPosition(VisibleInWindow);
      }
      if (event == vtkCommand::PickEvent || event == vtkCommand::StartPickEvent)
      {
        if (mw->m_drawShapeMode == SHAPE_MODE_NONE || mw->getActiveTabId() != TAB_DRAW) 
        {
          moveCursor(VisibleInWindow, x, y);
        }
        else
        {
          makeStroke(VisibleInWindow, x, y);
        }
      }
    }

    if (event == vtkCommand::WindowLevelEvent && mStartSlicer >= 0 && this->SM->GetPreset() <= 1)
    {

      // Adjust the window level here
      int *size = isi->GetInteractor()->GetRenderWindow()->GetSize();
      double window = this->InitialWindow;
      double level = this->InitialLevel;
      double range[2];
      this->SM->GetImage()->GetScalarRange(range);

      // Compute normalized delta
      double dx = static_cast<double>(isi->GetWindowLevelCurrentPosition()[0] - isi->GetWindowLevelStartPosition()[0]) / size[0];
      double dy = static_cast<double>(isi->GetWindowLevelStartPosition()[1] - isi->GetWindowLevelCurrentPosition()[1]) / size[1];
      //Window is exponential in nature, use exponential to avoid falling into negative numbers
      dx = std::exp(1.0 * (dx*fabs(dx) + dx)); //Quadratic behavior for more reactive interface
      dy = 0.15 * (dy*fabs(dy) + dy) * (range[1] - range[0]);//Quadratic behavior for more reactive interface

      this->SM->SetColorWindow(window*dx);
      this->SM->SetColorLevel(level - dy);
      this->SM->SetPreset(PRESET_USER);
      this->SM->Render();
      this->SM->UpdateWindowLevel();
      return;
    }
  }
}
void SlicerManagerCommand::AddActions()
{
}

void SlicerManagerCommand::Dolly(double factor, vtkRenderWindowInteractor *interactor)
{
  int VisibleInWindow = this->FindSlicerNumber(interactor->GetRenderWindow());
  vtkRenderer* renderer;
  if (VisibleInWindow > -1)
  {
    renderer = this->SM->GetSlicer(VisibleInWindow)->GetRenderer();
  }
  else
  {
    return;
  }
  vtkCamera *camera = renderer->GetActiveCamera();
  if (camera->GetParallelProjection())
  {
    camera->SetParallelScale(camera->GetParallelScale() / factor);
  }
  else
  {
    camera->Dolly(factor);
  }
  if (interactor->GetLightFollowCamera())
  {
    renderer->UpdateLightsGeometryToFollowCamera();
  }
  renderer->ResetCameraClippingRange();
}

void SlicerManagerCommand::moveCursor(int VisibleInWindow, double x, double y)
{

  vtkRenderer* renderer = this->SM->GetSlicer(VisibleInWindow)->GetRenderer();
  renderer->DisplayToNormalizedDisplay(x, y);
  renderer->NormalizedDisplayToViewport(x, y);
  renderer->ViewportToNormalizedViewport(x, y);
  double z = 0;
  renderer->NormalizedViewportToView(x, y, z);
  renderer->ViewToWorld(x, y, z);

  auto origin = this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin();
  auto spacing = this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing();

  double X = (x - origin[0]) / spacing[0];
  double Y = (y - origin[1]) / spacing[1];
  double Z = (z - origin[2]) / spacing[2];

  switch (this->SM->GetSlicer(VisibleInWindow)->GetSliceOrientation()) {
  case vtkImageViewer2::SLICE_ORIENTATION_XY:
    X = ROUND(X);
    Y = ROUND(Y);
    Z = this->SM->GetSlicer(VisibleInWindow)->GetSlice();
    break;
  case vtkImageViewer2::SLICE_ORIENTATION_XZ:
    X = ROUND(X);
    Y = this->SM->GetSlicer(VisibleInWindow)->GetSlice();
    Z = ROUND(Z);
    break;
  case vtkImageViewer2::SLICE_ORIENTATION_YZ:
    X = this->SM->GetSlicer(VisibleInWindow)->GetSlice();
    Y = ROUND(Y);
    Z = ROUND(Z);
    break;
  }
  //
  double xWorld = X * spacing[0] + origin[0];
  double yWorld = Y * spacing[1] + origin[1];
  double zWorld = Z * spacing[2] + origin[2];

  //Update the cursur position 
  this->SM->GetSlicer(VisibleInWindow)->SetCurrentPosition(xWorld, yWorld, zWorld);
  this->SM->Picked();
  this->SM->UpdateViews(VisibleInWindow);
  this->SM->UpdateLinked(VisibleInWindow);
  this->SM->UpdateInfoOnCursorPosition(VisibleInWindow);

}

std::pair<int, int > SlicerManagerCommand::point3Dto2D(const PointVal& pt3D, const int orientation)
{
  std::pair<int, int> pt;
  if (orientation == 2)//SLICE_ORIENTATION_XY
  {
    pt.first = pt3D.x;
    pt.second = pt3D.y;
  }
  else if (orientation == 1)//SLICE_ORIENTATION_XZ
  {
    pt.first = pt3D.x;
    pt.second = pt3D.z;
  }
  else//SLICE_ORIENTATION_YZ ==0
  {
    pt.first = pt3D.y;
    pt.second = pt3D.z;
  }
  return pt;
}

PointVal SlicerManagerCommand::point2Dto3D(const std::pair<int, int >& pt, const int orientation, const int slice, const int value)
{
  PointVal pt3D;
  if (orientation == 2)//SLICE_ORIENTATION_XY
  {
    pt3D.x = pt.first;
    pt3D.y = pt.second;
    pt3D.z = value;
  }
  else if (orientation == 1)//SLICE_ORIENTATION_XZ
  {
    pt3D.x = pt.first;
    pt3D.y = value;
    pt3D.z = pt.second;
  }
  else//SLICE_ORIENTATION_YZ ==0
  {
    pt3D.x = value;
    pt3D.y = pt.first;
    pt3D.z = pt.second;
  }
  pt3D.value = value;
  return pt3D;
}
std::vector< std::pair<int, int> > SlicerManagerCommand::points3Dto2D(const std::vector<PointVal>& points3D, const int orientation)
{
  std::vector< std::pair<int, int> > points2D;
  for (size_t i = 0; i < points3D.size(); i++)
  {
    points2D.push_back(point3Dto2D(points3D[i], orientation));
  }
  return points2D;
}
std::vector<PointVal> SlicerManagerCommand::points2Dto3D(const std::vector< std::pair<int, int> >& points2D, const int orientation, const int slice, const int value)
{
  std::vector<PointVal> points3D;
  for (size_t i = 0; i < points2D.size(); i++)
  {
    points3D.push_back(point2Dto3D(points2D[i], orientation, slice, value));
  }
  return points3D;
}
std::vector< std::pair<int, int> > SlicerManagerCommand::getCirclePoints(int x, int y, int r)
{
  if (r <= 0) r = 1;

  std::vector< std::pair<int, int> > points;
  double dtheta = 2 * M_PI / 8 / r;
  int n = 2 * M_PI / dtheta;
  points.push_back(std::make_pair(x + r, y));
  for (int i = 1; i <= n; i++)
  {
    double theta = i * dtheta;
    int x1 = int(x + r * cos(theta) + 0.5);
    int y1 = int(y + r * sin(theta) + 0.5);
    points.push_back(std::make_pair(x1, y1));
  }
  points.push_back(points[0]);
  return points;
}
std::vector< std::pair<int, int> > SlicerManagerCommand::getLinePoints(int x0, int y0, int x1, int y1, float wd)
{
  std::vector< std::pair<int, int> > points;
  int dx = abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx - dy, e2, x2, y2;                          /* error value e_xy */
  float ed = dx + dy == 0 ? 1 : sqrt((float)dx*dx + (float)dy*dy);

  for (wd = (wd + 1) / 2;;)
  {                                   /* pixel loop */
    points.push_back(std::make_pair(x0, y0));
    e2 = err; x2 = x0;
    if (2 * e2 >= -dx)
    {                                           /* x step */
      for (e2 += dy, y2 = y0; e2 < ed*wd && (y1 != y2 || dx > dy); e2 += dx)
      {
        points.push_back(std::make_pair(x0, y2 += sy));
      }
      if (x0 == x1) break;
      e2 = err; err -= dy; x0 += sx;
    }
    if (2 * e2 <= dy)
    {                                            /* y step */
      for (e2 = dx - e2; e2 < ed*wd && (x1 != x2 || dx < dy); e2 += dy)
      {
        points.push_back(std::make_pair(x2 += sx, y0));
      }

      if (y0 == y1) break;
      err += dx; y0 += sy;
    }
  }
  return points;
}
std::vector<PointVal> SlicerManagerCommand::drawLine(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int width)
{
  std::vector<PointVal> undoBuffer;
  std::pair<int, int> p1 = point3Dto2D(startPt, orientation);
  std::pair<int, int> p2 = point3Dto2D(endPt, orientation);
  std::vector< std::pair<int, int> > points = getLinePoints(p1.first, p1.second, p2.first, p2.second, width);
  std::vector<PointVal> points3D = points2Dto3D(points, orientation, slice, startPt.value);

  for (size_t i = 0; i < points3D.size(); i++)
  {
    PointVal undoPt = drawPoint(points3D[i], image, orientation, slice);
    if (undoPt.isValid())
    {
      undoBuffer.push_back(undoPt);
    }
  }
  return undoBuffer;
}

std::vector<PointVal> SlicerManagerCommand::fillShape(PointVal centerPt, int labelToFill, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice)
{
  std::vector<PointVal> undoBuffer;
  using MaskPixelType = unsigned char;
  using MaskImageType = itk::Image< MaskPixelType, 3 >;
  using MaskSliceType = itk::Image< MaskPixelType, 2 >;
  using ConnectedFilterType = itk::ConfidenceConnectedImageFilter< MaskSliceType, MaskSliceType >;
  using ExtractFilterType = itk::ExtractImageFilter< MaskImageType, MaskSliceType >;
  using RescaleIntensityImageFilter = itk::RescaleIntensityImageFilter< MaskImageType >;

  MaskImageType::Pointer maskItk = MaskImageType::New(); // putting the itk::Image version of 'image' here
  maskItk = convertVtkToItk< float, MaskPixelType, MaskImageType::ImageDimension >(image);
  maskItk->DisconnectPipeline();

#ifndef CAPTK_PACKAGE_PROJECT
  cbica::WriteImage< MaskImageType >(maskItk, cbica::getExecutablePath() + "/tempMask_1.nii.gz");
#endif

  auto pointToConsider = point3Dto2D(centerPt, orientation);
  auto maskFullSize = maskItk->GetLargestPossibleRegion().GetSize();

  itk::ImageRegionIterator< MaskImageType > maskIt(maskItk, maskItk->GetLargestPossibleRegion());
  MaskImageType::IndexType indexOnMaskToConsider;
  //for (; !maskIt.IsAtEnd(); ++maskIt)
  //{
  //  if (maskIt.Get() == startPt.value)
  //  {
  //    indexOnMaskToConsider = maskIt.GetIndex();
  //    break;
  //  }
  //}

  //auto minMaxFilter_full = itk::MinimumMaximumImageCalculator< MaskImageType >::New();
  //minMaxFilter_full->SetImage(maskItk);
  //minMaxFilter_full->Compute();
  //auto maxValue = minMaxFilter_full->GetMaximum();

  auto normalizeFilter = RescaleIntensityImageFilter::New();
  normalizeFilter->SetInput(maskItk);
  normalizeFilter->SetOutputMaximum(255);
  normalizeFilter->SetOutputMinimum(0);
  normalizeFilter->Update();

  MaskImageType::SizeType desiredSize; // what is the size of the slice that's needed
  MaskImageType::IndexType desiredIndex; // the index of maskItk from where extraction needs to happen
  desiredSize.Fill(0); // done because we want a 2D slice of a particular index
  desiredIndex.Fill(0); // done because we want the entire expanse of a slice
  MaskSliceType::IndexType seedPoint; // seed point for the ConnectedFilterType to work on

  seedPoint[0] = pointToConsider.first; // this is done to use the next voxel location as a seed, can work with '1' as well but using '3' for safety
  seedPoint[1] = pointToConsider.second;

  //const unsigned int incrementForSeed = 2;
  switch (orientation)
  {
  case 1: // SLICE_ORIENTATION_XZ
  {
    // get slice expanse
    desiredSize[0] = maskFullSize[0];
    desiredSize[2] = maskFullSize[2];

    indexOnMaskToConsider[0] = seedPoint[0];
    indexOnMaskToConsider[1] = slice;
    indexOnMaskToConsider[2] = seedPoint[1];

    // get index of the slice for extraction
    desiredIndex[1] = /*indexOnMaskToConsider[1]*/slice;
    break;
  }
  case 2: // SLICE_ORIENTATION_XY
  {
    // get slice expanse
    desiredSize[0] = maskFullSize[0];
    desiredSize[1] = maskFullSize[1];

    indexOnMaskToConsider[0] = seedPoint[0];
    indexOnMaskToConsider[1] = seedPoint[1];
    indexOnMaskToConsider[2] = slice;

    // get index of the slice for extraction
    desiredIndex[2] = /*indexOnMaskToConsider[2]*/slice;
    break;
  }
  default: // SLICE_ORIENTATION_YZ
  {
    // get slice expanse
    desiredSize[1] = maskFullSize[1];
    desiredSize[2] = maskFullSize[2];

    indexOnMaskToConsider[0] = slice;
    indexOnMaskToConsider[1] = seedPoint[0];
    indexOnMaskToConsider[2] = seedPoint[1];

    // get index of the slice for extraction
    desiredIndex[0] = /*indexOnMaskToConsider[0]*/slice;
    break;
  }
  }

  // extraction of the slice
  auto extractFilter = ExtractFilterType::New();
  extractFilter->SetExtractionRegion(MaskImageType::RegionType(desiredIndex, desiredSize));
  extractFilter->SetInput(normalizeFilter->GetOutput());
  extractFilter->SetDirectionCollapseToSubmatrix(); // This is required; see documentation for details
  extractFilter->Update(); // at this point, we have extracted the desired slice; it no longer has any connection to base image at this point

  // initialize and run the connected component filter; tried with a bunch of different filters and this seemed to work the best
  //auto rescaledPixelValueToConsider = normalizeFilter->GetOutput()->GetPixel(indexOnMaskToConsider);
  auto connectedFilter = ConnectedFilterType::New();
  connectedFilter->SetInput(extractFilter->GetOutput());
  connectedFilter->SetInitialNeighborhoodRadius(1); // only the first connected region to be considered
  connectedFilter->SetNumberOfIterations(0); // want to finish after the first run itself
  connectedFilter->SetMultiplier(2.5); // confidence interval is the mean plus or minus the "Multiplier" times the standard deviation; 2.5 => 99% of samples are included
  connectedFilter->SetSeed(seedPoint);
  connectedFilter->SetReplaceValue(/*startPt.value*//*rescaledPixelValueToConsider*/labelToFill); // replace with same value as the drawing
  connectedFilter->Update();

  //auto minMaxFilter_slice = itk::MinimumMaximumImageCalculator< MaskSliceType >::New();
  //minMaxFilter_slice->SetImage(connectedFilter->GetOutput());
  //minMaxFilter_slice->Compute();
  //auto maxIndex = minMaxFilter_slice->GetIndexOfMaximum();

#ifndef CAPTK_PACKAGE_PROJECT
    //cbica::WriteImage< MaskImageType >(maskItk, cbica::getExecutablePath() + "/tempMask_1.nii.gz");
    //cbica::WriteImage<MaskSliceType>(extractFilter->GetOutput(), cbica::getExecutablePath() + "/captk_outputMask_" + std::to_string(orientation) + "_extracted.nii.gz");
    //cbica::WriteImage<MaskSliceType>(connectedFilter->GetOutput(), cbica::getExecutablePath() + "/captk_outputMask_" + std::to_string(orientation) + "_connected.nii.gz");
#endif

    // iterate through the 2D slice and 3D slice in tandem to fill the latter up. There is a PastImageFilter but it wasn't working
  itk::ImageRegionIterator< MaskSliceType > maskSliceIt(connectedFilter->GetOutput(), connectedFilter->GetOutput()->GetLargestPossibleRegion());
  MaskImageType::IndexType indexToPopulate;

  // this is done so that the iterator doesn't start from zero but starts from the location where the seed was set
  //seedPoint[0] = seedPoint[0] - (incrementForSeed /*+ 1*/); // the '1' here is for difference between VTK and ITK image spaces
  //seedPoint[1] = seedPoint[1] - (incrementForSeed /*+ 1*/);
  //maskSliceIt.SetIndex(seedPoint);
  maskSliceIt.GoToBegin();
  //if (((maxIndex[0] == seedPoint[0] + 1) && ((maxIndex[1] == seedPoint[1] + 1))) || // proceed if maxIndex = seedPoint
  //  ((maxIndex[0] == seedPoint[0]) && (maxIndex[1] == seedPoint[1])))
  {
    // cannot combine maskIt increment here since the increment can happen over either x, y or z directions
    for (; !maskSliceIt.IsAtEnd(); ++maskSliceIt)
    {
      PointVal tempBuff; // temporary undo buffer
      if (maskSliceIt.Get() == /*startPt.value*//*rescaledPixelValueToConsider*/labelToFill) // only do index processing for values that aren't '0'
      {
        auto index = maskSliceIt.GetIndex(); // index of the 2D slice
        switch (orientation)
        {
        case 1:
        {
          // populate information from 2D slice and startPt to the main image
          indexToPopulate[0] = index[0];
          indexToPopulate[1] = indexOnMaskToConsider[1];
          indexToPopulate[2] = index[1];
          break;
        }
        case 2:
        {
          // populate information from 2D slice and startPt to the main image
          indexToPopulate[0] = index[0];
          indexToPopulate[1] = index[1];
          indexToPopulate[2] = indexOnMaskToConsider[2];
          break;
        }
        default:
        {
          // populate information from 2D slice and startPt to the main image
          indexToPopulate[0] = indexOnMaskToConsider[0];
          indexToPopulate[1] = index[0];
          indexToPopulate[2] = index[1];
          break;
        }
        }

        // populate temporary undo buffer - this thing doesn't work
        tempBuff.x = indexToPopulate[0];
        tempBuff.y = indexToPopulate[1];
        tempBuff.z = indexToPopulate[2];
        tempBuff.value = /*centerPt.value*/labelToFill;
        maskIt.SetIndex(indexToPopulate);
        auto tempVal = maskIt.Get();

        drawPoint(tempBuff, image, orientation, slice);

        tempBuff.value = tempVal/*centerPt.value*/; // this is done to ensure that the previous label value gets picked up for undo functionality
        undoBuffer.push_back(tempBuff);
        //maskIt.SetIndex(indexToPopulate);
        //maskIt.Set(/*centerPt.value*/labelToFill);
      }
    }
  }

  //cbica::WriteImage< MaskImageType >(maskItk, saveDir + "/filledMask.nii.gz"); 

  return undoBuffer;
}

std::vector<PointVal> SlicerManagerCommand::drawRectangle(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int thickness)
{
  PointVal corner1 = startPt;
  PointVal corner2 = startPt;
  switch (orientation)
  {
  case 1: // SLICE_ORIENTATION_XZ
  {
    corner1.x = endPt.x;
    corner2.z = endPt.z;
    break;
  }
  case 2: // SLICE_ORIENTATION_XY
  {
    corner1.x = endPt.x;
    corner2.y = endPt.y;
    break;
  }
  default: // SLICE_ORIENTATION_YZ
  {
    corner1.y = endPt.y;
    corner2.z = endPt.z;
    break;
  }
  }

  std::vector<PointVal> undoBuffer;
  /*std::vector<PointVal>*/ undoBuffer = drawLine(startPt, corner1, image, orientation, slice, thickness);
  std::vector<PointVal> vec2 = drawLine(corner1, endPt, image, orientation, slice, thickness);
  std::vector<PointVal> vec3 = drawLine(endPt, corner2, image, orientation, slice, thickness);
  std::vector<PointVal> vec4 = drawLine(corner2, startPt, image, orientation, slice, thickness);

  undoBuffer.insert(undoBuffer.end(), vec2.begin(), vec2.end());
  undoBuffer.insert(undoBuffer.end(), vec3.begin(), vec3.end());
  undoBuffer.insert(undoBuffer.end(), vec4.begin(), vec4.end());

  std::vector<PointVal> vec_fill;
  switch (orientation)
  {
  case 1:
  {
    for (int x = std::min(startPt.x, endPt.x); x < std::max(startPt.x, endPt.x); x++)
    {
      for (int z = std::min(startPt.z, endPt.z); z < std::max(startPt.z, endPt.z); z++)
      {
        vec_fill = drawLine(PointVal(x, startPt.y, z, startPt.value), PointVal(x, endPt.y, z, endPt.value), image, orientation, slice, thickness);
        undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
      }
    }
    break;
  }
  case 2:
  {
    for (int x = std::min(startPt.x, endPt.x); x < std::max(startPt.x, endPt.x); x++)
    {
      for (int y = std::min(startPt.y, endPt.y); y < std::max(startPt.y, endPt.y); y++)
      {
        vec_fill = drawLine(PointVal(x, y, startPt.z, startPt.value), PointVal(x, y, endPt.z, endPt.value), image, orientation, slice, thickness);
        undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
      }
    }
    break;
    break;
  }
  default:
  {
    for (int y = std::min(startPt.y, endPt.y); y < std::max(startPt.y, endPt.y); y++)
    {
      for (int z = std::min(startPt.z, endPt.z); z < std::max(startPt.z, endPt.z); z++)
      {
        vec_fill = drawLine(PointVal(startPt.x, y, z, startPt.value), PointVal(endPt.x, y, z, endPt.value), image, orientation, slice, thickness);
        undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
      }
    }
    break;
    break;
  }
  }

  //std::vector<PointVal> vec_fill = fillShape(image, orientation, slice, saveDir, startPt);

  //undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());

  return undoBuffer;
}
std::vector<PointVal> SlicerManagerCommand::drawSphere(PointVal startPt, vtkSmartPointer<vtkImageData> image, int radius)
{
  std::vector<PointVal> undoBuffer;
  auto dims = image->GetDimensions();

  auto actualRadius = (radius * 2 + 1);
  // Iterate through phi, theta then convert r,theta,phi to  XYZ
  int x_max = 0, x_min = dims[0], y_max = 0, y_min = dims[1], z_max = 0, z_min = dims[2];
  //double dtheta = 2 * M_PI / 8 / actualRadius;
  //int n = 2 * M_PI / dtheta;
  //auto deltaTheta = M_PI / 12;
  //auto deltaPhi = 2 * M_PI / 10;
  for (double phi = 0.; phi < 2 * M_PI; phi += M_PI / 10.) // Azimuth [0, 2PI]
  {
    for (double theta = 0.; theta < 2 * M_PI; theta += M_PI / 10.) // Elevation [0, PI]
    {
      PointVal tempPoint;
      tempPoint.x = std::round(actualRadius * cos(phi) * sin(theta)) + startPt.x;
      tempPoint.y = std::round(actualRadius * sin(phi) * sin(theta)) + startPt.y;
      tempPoint.z = std::round(actualRadius * cos(theta)) + startPt.z;
      tempPoint.value = startPt.value;
      auto undoPt = drawPoint(tempPoint, image);
      if (undoPt.isValid())
      {
        undoBuffer.push_back(tempPoint);
        if (x_max < tempPoint.x)
        {
          x_max = tempPoint.x;
        }
        if (y_max < tempPoint.y)
        {
          y_max = tempPoint.y;
        }
        if (z_max < tempPoint.z)
        {
          z_max = tempPoint.z;
        }
        if (x_min > tempPoint.x)
        {
          x_min = tempPoint.x;
        }
        if (y_min > tempPoint.y)
        {
          y_min = tempPoint.y;
        }
        if (z_min > tempPoint.z)
        {
          z_min = tempPoint.z;
        }
      }
    }
  }

  for (int x = x_min; x <= x_max; x++)
  {
    for (int y = y_min; y <= y_max; y++)
    {
      for (int z = z_min; z <= z_max; z++)
      {
        PointVal tempPoint;
        tempPoint.x = x;
        tempPoint.y = y;
        tempPoint.z = z;
        tempPoint.value = startPt.value;
        if (tempPoint.getDistanceFrom(startPt) < actualRadius)
        {
          auto undoPt = drawPoint(tempPoint, image);
          if (undoPt.isValid())
          {
            undoBuffer.push_back(tempPoint);
          }
        }
      }
    }
  }
  //for (int x = startPt.x - actualRadius; x <= startPt.x + actualRadius; x++)
  //{
  //  for (int y = startPt.y - actualRadius; y <= startPt.y + actualRadius; y++)
  //  {
  //    for (int z = startPt.z - actualRadius; z <= startPt.z + actualRadius; z++)
  //    {
  //      if ((x < dims[0]) && (y < dims[1]) && (z < dims[2]))
  //      {
  //        PointVal tempPoint;
  //        tempPoint.x = x;
  //        tempPoint.y = y;
  //        tempPoint.z = z;
  //        tempPoint.value = startPt.value;
  //        if (tempPoint.getDistanceFrom(startPt) <= actualRadius) // remove this if to get a cube
  //        {
  //          auto undoPt = drawPoint(tempPoint, image);
  //          if (undoPt.isValid())
  //          {
  //            undoBuffer.push_back(tempPoint);
  //          }
  //        }
  //      }
  //    }
  //  }
  //}

  return undoBuffer;
}
std::vector<PointVal> SlicerManagerCommand::drawCircle(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int width)
{
  std::vector<PointVal> undoBuffer;
  std::pair<int, int> p1 = point3Dto2D(startPt, orientation);
  std::pair<int, int> p2 = point3Dto2D(endPt, orientation);
  std::pair<int, int> center;
  center.first = (p1.first + p2.first) / 2;
  center.second = (p1.second + p2.second) / 2;
  int radius = sqrt((p1.first - p2.first)*(p1.first - p2.first) + (p1.second - p2.second)*(p1.second - p2.second)) / 2;
  std::vector< std::pair<int, int> > circlePoints = getCirclePoints(center.first, center.second, radius);
  std::vector< std::pair<int, int> > thickCirclePoints;
  for (size_t i = 1; i < circlePoints.size(); i++)
  {
    std::vector< std::pair<int, int> > points = getLinePoints(circlePoints[i - 1].first, circlePoints[i - 1].second, circlePoints[i].first, circlePoints[i].second, width);
    thickCirclePoints.insert(thickCirclePoints.end(), points.begin(), points.end());
  }
  //Lets remove duplicates 
  sort(thickCirclePoints.begin(), thickCirclePoints.end());
  thickCirclePoints.erase(std::unique(thickCirclePoints.begin(), thickCirclePoints.end()), thickCirclePoints.end());
  std::vector<PointVal> points3D = points2Dto3D(thickCirclePoints, orientation, slice, startPt.value);
  //std::vector<PointVal> vec_fill = fillShape(image, orientation, slice, saveDir, startPt);
  std::vector<PointVal> vec_fill;

  for (size_t i = 0, j = points3D.size() - 1; i < std::floor(points3D.size() / 2); i++, j--)
  {
    switch (orientation)
    {
    case 1:
    {
      vec_fill = drawLine(PointVal(points3D[i].x, startPt.y, points3D[i].z, startPt.value), PointVal(points3D[j].x, endPt.y, points3D[j].z, endPt.value), image, orientation, slice, 3);
      undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
      break;
    }
    case 2:
    {
      vec_fill = drawLine(PointVal(points3D[i].x, points3D[i].y, startPt.z, startPt.value), PointVal(points3D[j].x, points3D[j].y, endPt.z, endPt.value), image, orientation, slice, 3);
      undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
      break;
    }
    default:
    {
      vec_fill = drawLine(PointVal(startPt.x, points3D[i].y, points3D[i].z, startPt.value), PointVal(endPt.x, points3D[j].y, points3D[j].z, endPt.value), image, orientation, slice, 3);
      undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
      break;
    }
    }
  }
  //for (size_t i = 0; i < points3D.size(); i++)
  //{
  //  PointVal undoPt = drawPoint(points3D[i], image, orientation, slice);
  //  if (undoPt.isValid())
  //  {
  //    undoBuffer.push_back(undoPt);
  //  }
  //}
  undoBuffer.insert(undoBuffer.end(), vec_fill.begin(), vec_fill.end());
  return undoBuffer;
}
PointVal SlicerManagerCommand::drawPoint(PointVal pt, vtkSmartPointer<vtkImageData> image)
{
  if (image == NULL)
  {
    return pt.getInvalidPt();
  }
  //Fix for crash check range before accessing pixel 
  if (pt.isWithinRange(image->GetDimensions()))
  {
    float* pData = (float*)image->GetScalarPointer(pt.x, pt.y, pt.z);
    int oldVal = *pData;
    *pData = (float)pt.value;
    pt.value = oldVal;
    return pt;
  }
  return pt.getInvalidPt();
}
PointVal SlicerManagerCommand::drawPoint(PointVal pt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, const int i, const int j)
{
  if (image == NULL)
  {
    return pt.getInvalidPt();
  }
  if (orientation == 2)//SLICE_ORIENTATION_XY
  {
    pt.x = pt.x + i;
    pt.y = pt.y + j;
    pt.z = slice;
  }
  else if (orientation == 1)//SLICE_ORIENTATION_XZ
  {
    pt.x = pt.x + i;
    pt.y = slice;
    pt.z = pt.z + j;
  }
  else//SLICE_ORIENTATION_YZ
  {
    pt.x = slice;
    pt.y = pt.y + i;
    pt.z = pt.z + j;
  }
  //Fix for crash check range before accessing pixel 
  if (pt.isWithinRange(image->GetDimensions()))
  {
    float* pData = (float*)image->GetScalarPointer(pt.x, pt.y, pt.z);
    int oldVal = *pData;
    *pData = (float)pt.value;
    pt.value = oldVal;
    return pt;
  }
  return pt.getInvalidPt();

}
void SlicerManagerCommand::makeStroke(int VisibleInWindow, double x, double y)
{
  vtkRenderer* renderer = this->SM->GetSlicer(VisibleInWindow)->GetRenderer();
  renderer->DisplayToNormalizedDisplay(x, y);
  renderer->NormalizedDisplayToViewport(x, y);
  renderer->ViewportToNormalizedViewport(x, y);
  double z = 0;
  renderer->NormalizedViewportToView(x, y, z);
  renderer->ViewToWorld(x, y, z);

  auto origin = this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetOrigin();
  auto spacing = this->SM->GetSlicer(VisibleInWindow)->GetInput()->GetSpacing();

  double X = (x - origin[0]) / spacing[0];
  double Y = (y - origin[1]) / spacing[1];
  double Z = (z - origin[2]) / spacing[2];

  int orientation = this->SM->GetSlicer(VisibleInWindow)->GetSliceOrientation();
  int slice = this->SM->GetSlicer(VisibleInWindow)->GetSlice();
  //int* dims = this->SM->mMask->GetDimensions();

  int color = mw->getSelectedDrawLabel();
  PointVal centerPt(ROUND(X), ROUND(Y), ROUND(Z), color);
  if (mw->m_drawShapeMode == SHAPE_MODE_FREE_HAND || mw->m_drawShapeMode == SHAPE_MODE_ERASER)
  {
    if (mw->m_drawShapeMode == SHAPE_MODE_ERASER)
    {
      centerPt.value = 0;
    }
    int size = mw->getSelectedDrawSize();
    for (int j = -size; j <= size; j++)
    {
      for (int i = -size; i <= size; i++)
      {
        PointVal pt = drawPoint(centerPt, this->SM->mMask, orientation, slice, i, j);
        m_undoBuffer.push_back(pt);
      }
    }

  }
  else // for all 'shape' routines
  {
    //Undo previous shape  draw
    if (m_shapeBuffer.empty())
    {
      m_startPoint = centerPt;
    }
    else
    {
      for (std::vector<PointVal>::iterator it = m_shapeBuffer.end(); it != m_shapeBuffer.begin();)
      {
        --it;
        drawPoint(*it, this->SM->mMask, orientation, slice);
      }
    }
    switch (mw->m_drawShapeMode)
    {
      //case SHAPE_MODE_FREE_HAND:
      // already handled previously
    case SHAPE_MODE_LINE:
    {
      m_shapeBuffer = drawLine(m_startPoint, centerPt, this->SM->mMask, orientation, slice, mw->getSelectedDrawSize());
      break;
    }
    case SHAPE_MODE_CIRCLE:
    {
      m_shapeBuffer = drawCircle(m_startPoint, centerPt, this->SM->mMask, orientation, slice, mw->getSelectedDrawSize());
      break;
    }
    case SHAPE_MODE_RECTANGLE:
    {
      m_shapeBuffer = drawRectangle(m_startPoint, centerPt, this->SM->mMask, orientation, slice, mw->getSelectedDrawSize());
      break;
    }
    case SHAPE_MODE_FILL:
    {
      m_shapeBuffer = fillShape(centerPt, mw->getSelectedDrawLabel(), this->SM->mMask, orientation, slice);
      break;
    }
    case SHAPE_MODE_SPHERE:
    {
      if (orientation == 2)//SLICE_ORIENTATION_XY
      {
        centerPt.z = slice;
      }
      else if (orientation == 1)//SLICE_ORIENTATION_XZ
      {
        centerPt.y = slice;
      }
      else//SLICE_ORIENTATION_YZ
      {
        centerPt.x = slice;
      }
      m_shapeBuffer = drawSphere(centerPt, this->SM->mMask, mw->getSelectedDrawSize());
      break;
    }
    default:
    {
      break;
    }
    }
    //if (mw->m_drawShapeMode == SHAPE_MODE_CIRCLE)
    //{
    //  m_shapeBuffer = drawCircle(m_startPoint, centerPt, this->SM->mMask, orientation, slice, mw->getSelectedDrawSize());
    //}
    //else if (mw->m_drawShapeMode == SHAPE_MODE_RECTANGLE)
    //{
    //  m_shapeBuffer = drawRectangle(m_startPoint, centerPt, this->SM->mMask, orientation, slice, mw->getSelectedDrawSize());
    //}
    //else if (mw->m_drawShapeMode == SHAPE_MODE_LINE)
    //{
    //  m_shapeBuffer = drawLine(m_startPoint, centerPt, this->SM->mMask, orientation, slice, mw->getSelectedDrawSize());
    //}
    //else if (mw->m_drawShapeMode == SHAPE_MODE_LINE)
    //{
    //  int size = mw->getSelectedDrawSize();
    //  for (int j = -size; j <= size; j++)
    //  {
    //    for (int i = -size; i <= size; i++)
    //    {
    //      PointVal pt = drawPoint(centerPt, this->SM->mMask, orientation, slice, i, j);
    //      m_undoBuffer.push_back(pt);
    //    }
    //  }
    //}
    m_undoBuffer.insert(m_undoBuffer.end(), m_shapeBuffer.begin(), m_shapeBuffer.end());
  }

  //if (!m_undoBuffer.empty())
  //{
  //  this->SM->ActionAdded(m_undoBuffer);
  //}

  this->SM->mMask->Modified();
  this->SM->GetSlicer(VisibleInWindow)->Render();

}