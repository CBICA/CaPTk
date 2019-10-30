/**
\file  InteractorStyleNavigator.h

\brief Declaration of InteractorStyleNavigator class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _InteractorStyleNavigator_h_
#define _InteractorStyleNavigator_h_


//#include "CAPTk.h"
#include "vtkInteractorStyle.h"

#define VTKIS_WINDOW_LEVEL 1024
//#define VTKIS_PICK         1025


/**
\class InteractorStyleNavigator

\brief Interactor related stuff
*/
class InteractorStyleNavigator : public vtkInteractorStyle
{
public:
  static InteractorStyleNavigator *New();
#if VTK_MAJOR_VERSION >= 6
  vtkTypeMacro(InteractorStyleNavigator, vtkInteractorStyle);
#else
  vtkTypeRevisionMacro(InteractorStyleNavigator, vtkInteractorStyle);
#endif

  vtkGetVector2Macro(WindowLevelStartPosition, int);
  vtkGetVector2Macro(WindowLevelCurrentPosition, int);

  virtual void OnMouseMove();
  virtual void OnLeftButtonDown();
  virtual void OnLeftButtonUp();
  virtual void OnRightButtonDown();
  virtual void OnRightButtonUp();
  virtual void OnMiddleButtonDown();
  virtual void OnMiddleButtonUp();
  virtual void OnEnter();
  virtual void OnLeave();
  virtual void OnMouseWheelForward();
  virtual void OnMouseWheelBackward();
  virtual void OnChar();
  virtual void WindowLevel();
  virtual void Pick();
  virtual void StartWindowLevel();
  virtual void EndWindowLevel();
  virtual void StartPick();
  virtual void EndPick();
  virtual void Dolly();
  virtual void Pan();
  virtual void FindPokedRenderer(int, int);
  bool isDoubleClick()
  {
    if (NumberOfClicks > 1)
    {
      NumberOfClicks = 0;
      return true;
    }
    return false;
  }
protected:
  InteractorStyleNavigator();
  ~InteractorStyleNavigator();

  static void ProcessEvents(vtkObject* object,
    unsigned long event,
    void* clientdata,
    void* calldata);

  double MotionFactor;

  virtual void Dolly(double factor);

  int WindowLevelStartPosition[2];
  int WindowLevelCurrentPosition[2];

private:
  InteractorStyleNavigator(const InteractorStyleNavigator&);
  void operator=(const InteractorStyleNavigator&);
private:
  unsigned int NumberOfClicks;
  int PreviousPosition[2];
  int ResetPixelDistance;
};


#endif
