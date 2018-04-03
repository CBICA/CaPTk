#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkImageViewer2.h>
#include <vtkInteractorStyleImage.h>
#include <vtkActor2D.h>
#include <vtkTextProperty.h>
#include <vtkTextMapper.h>
#include <sstream>


class StatusMessage
{
public:
  static std::string Format(int slice, int maxSlice)
  {
    std::stringstream tmp;
    tmp << "Slice Number  " << slice + 1 << "/" << maxSlice + 1;
    return tmp.str();
  }
};

class myVtkInteractorStyleImage : public vtkInteractorStyleImage
{
public:
  static myVtkInteractorStyleImage* New();
  vtkTypeMacro(myVtkInteractorStyleImage, vtkInteractorStyleImage);

protected:
  vtkImageViewer2* _ImageViewer;
  vtkTextMapper* _StatusMapper;
  int _Slice;
  int _MinSlice;
  int _MaxSlice;

public:
  void SetImageViewer(vtkImageViewer2* imageViewer)
  {
    _ImageViewer = imageViewer;
    _MinSlice = imageViewer->GetSliceMin();
    _MaxSlice = imageViewer->GetSliceMax();
    _Slice = _MinSlice;
    cout << "Slicer: Min = " << _MinSlice << ", Max = " << _MaxSlice << std::endl;
  }

  void SetStatusMapper(vtkTextMapper* statusMapper)
  {
    _StatusMapper = statusMapper;
  }


protected:
  void MoveSliceForward()
  {
    if (_Slice < _MaxSlice)
    {
      _Slice += 1;
      cout << "MoveSliceForward::Slice = " << _Slice << std::endl;
      _ImageViewer->SetSlice(_Slice);
      std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
      _StatusMapper->SetInput(msg.c_str());
      _ImageViewer->Render();
    }
  }

  void MoveSliceBackward()
  {
    if (_Slice > _MinSlice)
    {
      _Slice -= 1;
      cout << "MoveSliceBackward::Slice = " << _Slice << std::endl;
      _ImageViewer->SetSlice(_Slice);
      std::string msg = StatusMessage::Format(_Slice, _MaxSlice);
      _StatusMapper->SetInput(msg.c_str());
      _ImageViewer->Render();
    }
  }


  virtual void OnKeyDown()
  {
    std::string key = this->GetInteractor()->GetKeySym();
    if (key.compare("Up") == 0)
      MoveSliceForward();
    else if (key.compare("Down") == 0)
      MoveSliceBackward();
    vtkInteractorStyleImage::OnKeyDown();
  }


  virtual void OnMouseWheelForward()
  {
    MoveSliceForward();
  }


  virtual void OnMouseWheelBackward()
  {
    if (_Slice > _MinSlice)
      MoveSliceBackward();
  }
};