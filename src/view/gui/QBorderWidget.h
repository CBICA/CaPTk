#pragma once

#include <vtkSmartPointer.h>
#include <vtkBorderRepresentation.h>
#include <vtkBorderWidget.h>
#include <vtkCommand.h>
#include <vtkImageActor.h>
#include <vtkRenderer.h>

class vtkBorderCallback : public vtkCommand
{
public:
  vtkBorderCallback()
  {
  }

  //static vtkBorderCallback *New()
  //{
  //	return new vtkBorderCallback();
  //}

  virtual void Execute(vtkObject *caller, unsigned long, void*)
  {

    vtkBorderWidget *borderWidget = reinterpret_cast<vtkBorderWidget*>(caller);

    // Get the display coordinates of the two corners of the box
    vtkCoordinate* lowerLeftCoordinate = static_cast<vtkBorderRepresentation*>(borderWidget->GetRepresentation())->GetPositionCoordinate();
    double* lowerLeft;
    lowerLeft = lowerLeftCoordinate->GetComputedWorldValue(this->Renderer);
    std::cout << "Lower left coordinate: " << lowerLeft[0] << ", " << lowerLeft[1] << std::endl;
    lowerLeft[2] = 0;

    vtkCoordinate* upperRightCoordinate = static_cast<vtkBorderRepresentation*>(borderWidget->GetRepresentation())->GetPosition2Coordinate();
    double* upperRight;
    upperRight = upperRightCoordinate->GetComputedWorldValue(this->Renderer);
    std::cout << "Upper right coordinate: " << upperRight[0] << ", " << upperRight[1] << std::endl;
    upperRight[2] = 0;
  }

  void SetRenderer(vtkSmartPointer<vtkRenderer> ren) { this->Renderer = ren; }
  void SetImageActor(vtkSmartPointer<vtkImageActor> im) { this->ImageActor = im; }

private:
  vtkSmartPointer<vtkRenderer> Renderer;
  vtkSmartPointer<vtkImageActor> ImageActor;

};
