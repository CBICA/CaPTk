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
        if (KeyPress == "l") {
          this->SM->ToggleInterpolation();
          this->SM->Render();
          return;
        }
        if (KeyPress == "h") {
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
        if (KeyPress == "r") {
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