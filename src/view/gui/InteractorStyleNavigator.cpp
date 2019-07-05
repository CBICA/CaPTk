/**
\file  InteractorStyleNavigator.cpp

\brief Implementation of InteractorStyleNavigator class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "InteractorStyleNavigator.h"
#include <vtkWidgetEvent.h>
#include <vtkObjectFactory.h>
#include <vtkActor.h>
#include <vtkBorderRepresentation.h>
#include <vtkBorderWidget.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSphereSource.h>
#include <vtkRendererCollection.h>
#include <vtkCamera.h>
#include <vtkCallbackCommand.h>


#if VTK_MAJOR_VERSION < 6
vtkCxxRevisionMacro(InteractorStyleNavigator, "$Revision: 1.70 $");
#endif
vtkStandardNewMacro(InteractorStyleNavigator);

InteractorStyleNavigator::InteractorStyleNavigator()
{
  this->WindowLevelStartPosition[0] = 0;
  this->WindowLevelStartPosition[1] = 0;

  this->WindowLevelCurrentPosition[0] = 0;
  this->WindowLevelCurrentPosition[1] = 0;

  this->MotionFactor = 10.0;
  this->NumberOfClicks = 0;
  this->ResetPixelDistance = 5;
  this->PreviousPosition[0] = 0;
  this->PreviousPosition[1] = 0;
}

InteractorStyleNavigator::~InteractorStyleNavigator()
{
  CurrentRenderer = NULL;
}

void InteractorStyleNavigator::FindPokedRenderer(int dummy1, int dummy2)
{
  vtkRenderWindow * renwin = this->GetInteractor()->GetRenderWindow();
  renwin->GetRenderers()->InitTraversal();
  while (true) 
  {
    vtkRenderer* current = renwin->GetRenderers()->GetNextItem();
    if (current == NULL || current->GetDraw()) 
    {
      CurrentRenderer = current;
      return;
    }
  }
}

void InteractorStyleNavigator::StartWindowLevel()
{
  if (this->State != VTKIS_NONE) 
  {
    return;
  }
  this->StartState(VTKIS_WINDOW_LEVEL);
  this->InvokeEvent(vtkCommand::StartWindowLevelEvent, this);
}

void InteractorStyleNavigator::EndWindowLevel()
{
  if (this->State != VTKIS_WINDOW_LEVEL) 
  {
    return;
  }
  this->InvokeEvent(vtkCommand::EndWindowLevelEvent, this);
  this->StopState();
}

void InteractorStyleNavigator::StartPick()
{
  if (this->State != VTKIS_NONE) 
  {
    return;
  }
  this->StartState(VTKIS_PICK);
  this->InvokeEvent(vtkCommand::StartPickEvent, this);
}

void InteractorStyleNavigator::EndPick()
{
  if (this->State != VTKIS_PICK) 
  {
    return;
  }
  this->InvokeEvent(vtkCommand::EndPickEvent, this);
  this->StopState();
}

void InteractorStyleNavigator::OnMouseMove()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  switch (this->State) 
  {
  case VTKIS_WINDOW_LEVEL:
    this->FindPokedRenderer(x, y);
    this->WindowLevel();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_PICK:
    this->FindPokedRenderer(x, y);
    this->Pick();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_PAN:
    this->FindPokedRenderer(x, y);
    this->Pan();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  case VTKIS_DOLLY:
    this->FindPokedRenderer(x, y);
    this->Dolly();
    this->InvokeEvent(vtkCommand::InteractionEvent, NULL);
    break;

  default:
    this->InvokeEvent(vtkCommand::UserEvent, NULL);
    break;
  }
}

void InteractorStyleNavigator::OnEnter()
{
  switch (this->State) 
  {
  case VTKIS_WINDOW_LEVEL:
    break;

  case VTKIS_PICK:
    break;

  case VTKIS_PAN:
    break;

  default:
    this->InvokeEvent(vtkCommand::EnterEvent, NULL);
    break;
  }
}

void InteractorStyleNavigator::OnLeave()
{
  switch (this->State) 
  {
  case VTKIS_WINDOW_LEVEL:
    break;

  case VTKIS_PICK:
    break;

  case VTKIS_PAN:
    break;

  default:
    this->InvokeEvent(vtkCommand::LeaveEvent, NULL);
    break;
  }
}

void InteractorStyleNavigator::OnRightButtonDown()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->FindPokedRenderer(x, y);
  if (this->CurrentRenderer == NULL) 
  {
  	return;
  }

  this->GrabFocus(this->EventCallbackCommand);
  if (!this->Interactor->GetShiftKey() && !this->Interactor->GetControlKey()) 
  {
  	this->WindowLevelStartPosition[0] = x;
  	this->WindowLevelStartPosition[1] = y;
  	this->StartWindowLevel();
  } 
  else 
  {
  	this->Superclass::OnRightButtonDown();
  }
  this->InvokeEvent(vtkCommand::RightButtonPressEvent, this);
}

void InteractorStyleNavigator::OnRightButtonUp()
{
  switch (this->State) 
  {
  case VTKIS_WINDOW_LEVEL:
  	this->EndWindowLevel();
  	if (this->Interactor) 
    {
  		this->ReleaseFocus();
  	}
  	break;
  }

  this->Superclass::OnRightButtonUp();
}

void InteractorStyleNavigator::OnLeftButtonDown()
{
  int x = this->Interactor->GetEventPosition()[0];
  int y = this->Interactor->GetEventPosition()[1];

  this->NumberOfClicks++;
  int pickPosition[2];
  this->GetInteractor()->GetEventPosition(pickPosition);

  int xdist = pickPosition[0] - this->PreviousPosition[0];
  int ydist = pickPosition[1] - this->PreviousPosition[1];

  this->PreviousPosition[0] = pickPosition[0];
  this->PreviousPosition[1] = pickPosition[1];

  int moveDistance = (int)sqrt((double)(xdist*xdist + ydist*ydist));
  if (moveDistance > this->ResetPixelDistance)
  {
    this->NumberOfClicks = 1;
  }

  this->FindPokedRenderer(x, y);
  if (this->CurrentRenderer == NULL) 
  {
    return;
  }

  this->GrabFocus(this->EventCallbackCommand);
  if (this->Interactor->GetControlKey()) 
  {
    this->CurrentRenderer->GetRenderWindow()->SetCurrentCursor(8);
    this->StartPan();
  }
  else if (this->Interactor->GetShiftKey()) 
  {
    this->Superclass::OnLeftButtonDown();
  }
  else 
  {
    this->StartPick();
  }
}

void InteractorStyleNavigator::OnLeftButtonUp()
{
  switch (this->State) 
  {
  case VTKIS_PICK:
    this->EndPick();
    if (this->Interactor) 
    {
      this->ReleaseFocus();
    }
    break;
  case VTKIS_PAN:
    this->EndPan();
    if (this->Interactor) 
    {
      this->Interactor->GetRenderWindow()->SetCurrentCursor(0);
      this->ReleaseFocus();
    }
    break;
  }

  this->Superclass::OnLeftButtonUp();
}

void InteractorStyleNavigator::OnMiddleButtonDown()
{
  this->Superclass::OnMiddleButtonDown();
}

void InteractorStyleNavigator::OnMiddleButtonUp()
{
  this->Superclass::OnMiddleButtonUp();
}

void InteractorStyleNavigator::OnChar()
{
}

void InteractorStyleNavigator::OnMouseWheelForward()
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }
  this->GrabFocus(this->EventCallbackCommand);
  if (this->Interactor->GetControlKey())
  {
    this->StartDolly();
    double factor = this->MotionFactor * 0.2 * this->MouseWheelMotionFactor;
    this->Dolly(pow((double)1.1, factor));
    this->EndDolly();
  }
  this->ReleaseFocus();
  this->InvokeEvent(vtkCommand::MouseWheelForwardEvent, this);
}

void InteractorStyleNavigator::OnMouseWheelBackward()
{
  this->FindPokedRenderer(this->Interactor->GetEventPosition()[0], this->Interactor->GetEventPosition()[1]);
  if (this->CurrentRenderer == NULL)
  {
    return;
  }

  this->GrabFocus(this->EventCallbackCommand);
  if (this->Interactor->GetControlKey()) 
  {
    this->StartDolly();
    double factor = this->MotionFactor * -0.2 * this->MouseWheelMotionFactor;
    this->Dolly(pow((double)1.1, factor));
    this->EndDolly();
  }
  this->ReleaseFocus();
  this->InvokeEvent(vtkCommand::MouseWheelBackwardEvent, this);
}

void InteractorStyleNavigator::WindowLevel()
{
  vtkRenderWindowInteractor *rwi = this->Interactor;

  this->WindowLevelCurrentPosition[0] = rwi->GetEventPosition()[0];
  this->WindowLevelCurrentPosition[1] = rwi->GetEventPosition()[1];

  this->InvokeEvent(vtkCommand::WindowLevelEvent, this);
}

void InteractorStyleNavigator::Pick()
{
  this->InvokeEvent(vtkCommand::PickEvent, this);
}

void InteractorStyleNavigator::Pan()
{
  if (this->CurrentRenderer == NULL) 
  {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;

  double viewFocus[4], focalDepth, viewPoint[3];
  double newPickPoint[4], oldPickPoint[4], motionVector[3];

  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();
  camera->GetFocalPoint(viewFocus);
  this->ComputeWorldToDisplay(viewFocus[0], viewFocus[1], viewFocus[2],
    viewFocus);
  focalDepth = viewFocus[2];

  this->ComputeDisplayToWorld((double)rwi->GetEventPosition()[0], (double)rwi->GetEventPosition()[1], focalDepth, newPickPoint);

  this->ComputeDisplayToWorld((double)rwi->GetLastEventPosition()[0], (double)rwi->GetLastEventPosition()[1], focalDepth, oldPickPoint);

  motionVector[0] = oldPickPoint[0] - newPickPoint[0];
  motionVector[1] = oldPickPoint[1] - newPickPoint[1];
  motionVector[2] = oldPickPoint[2] - newPickPoint[2];

  camera->GetFocalPoint(viewFocus);
  camera->GetPosition(viewPoint);
  camera->SetFocalPoint(motionVector[0] + viewFocus[0], motionVector[1] + viewFocus[1], motionVector[2] + viewFocus[2]);
  camera->SetPosition(motionVector[0] + viewPoint[0], motionVector[1] + viewPoint[1], motionVector[2] + viewPoint[2]);

  if (rwi->GetLightFollowCamera())
  {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
  }

  this->InvokeEvent(vtkCommand::EndInteractionEvent, this);
}

void InteractorStyleNavigator::Dolly()
{
  if (this->CurrentRenderer == NULL) {
    return;
  }

  vtkRenderWindowInteractor *rwi = this->Interactor;
  double *center = this->CurrentRenderer->GetCenter();
  int dy = rwi->GetEventPosition()[1] - rwi->GetLastEventPosition()[1];
  double dyf = this->MotionFactor * (double)(dy) / (double)(center[1]);
  this->Dolly(pow((double)1.1, dyf));
}

void InteractorStyleNavigator::Dolly(double factor)
{
  if (this->CurrentRenderer == NULL)
  {
    return;
  }

  vtkCamera *camera = this->CurrentRenderer->GetActiveCamera();

  if (camera->GetParallelProjection()) 
  {
    camera->SetParallelScale(camera->GetParallelScale() / factor);
  }
  else {
    camera->Dolly(factor);
    if (this->AutoAdjustCameraClippingRange) 
    {
      this->CurrentRenderer->ResetCameraClippingRange();
    }
  }

  if (this->Interactor->GetLightFollowCamera()) 
  {
    this->CurrentRenderer->UpdateLightsGeometryToFollowCamera();
  }
  this->CurrentRenderer->ResetCameraClippingRange();

  this->InvokeEvent(vtkCommand::EndInteractionEvent, this);
}
