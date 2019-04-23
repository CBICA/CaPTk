#ifndef _ComparisonViewerCommand_h_
#define _ComparisonViewerCommand_h_


//#include "CAPTk.h"
#include "SlicerManager.h"
#include "CaPTkUtils.h"
#include "vtkCommand.h"
#include "vtkSmartPointer.h"

class fMainWindow;
class vtkRenderWindowInteractor;

#if __GNUC__
#pragma GCC visibility push(hidden)
#endif
class ComparisonViewerCommand : public vtkCommand
{
public:
  static ComparisonViewerCommand *New()
  {
    return new ComparisonViewerCommand;
  }

  void Execute(vtkObject *caller, unsigned long event, void *vtkNotUsed(callData));

  SlicerManager *SM;
  void Dolly(double factor, vtkRenderWindowInteractor *interactor);
  void SetSlicerNumber(int slicer) { mSlicerNumber = slicer; }
  void AddActions();
  void moveCursor(/*int VisibleInWindow, */double x, double y);

  static std::pair<int, int > point3Dto2D(const PointVal& pt3D, const int orientation);

  static PointVal point2Dto3D(const std::pair<int, int >& pt, const int orientation, const int slice, const int value = 0);

  static std::vector< std::pair<int, int> > points3Dto2D(const std::vector<PointVal>& points3D, const int orientation);

  static std::vector<PointVal> points2Dto3D(const std::vector< std::pair<int, int> >& points2D, const int orientation, const int slice, const int value = 0);

  static std::vector< std::pair<int, int> > getCirclePoints(int x, int y, int r);

  static std::vector< std::pair<int, int> > getLinePoints(int x0, int y0, int x1, int y1, float wd);

  static std::vector<PointVal> drawLine(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int width);

  static std::vector<PointVal> fillShape(PointVal centerPt, int labelToFill, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice);

static std::vector<PointVal> drawRectangle(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int thickness);

static std::vector<PointVal> drawCircle(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int width);

static std::vector<PointVal> drawSphere(PointVal startPt, vtkSmartPointer<vtkImageData> image, int radius);

static PointVal drawPoint(PointVal pt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, const int i = 0, const int j = 0);

static PointVal drawPoint(PointVal pt, vtkSmartPointer<vtkImageData> image);

void makeStroke(double x, double y);

//! Set Comparison Viewers
void SetComparisonViewers(std::vector<vtkSmartPointer<Slicer>>);

//! Set current comparison viewer
void SetCurrentViewer(vtkSmartPointer<Slicer>);

protected:
  ComparisonViewerCommand();
  ~ComparisonViewerCommand() {}

private:
  int FindSlicerNumber(vtkRenderWindow* renwin);
  double InitialWindow;
  double InitialLevel;
  int mStartSlicer;
  int mSlicerNumber;
  fMainWindow* mw;
  std::vector< PointVal > m_undoBuffer;
  std::vector< PointVal > m_shapeBuffer;
  PointVal m_startPoint;

  //! Comparison Viewers
  std::vector<vtkSmartPointer<Slicer>> m_ComparisonViewers;

  //! Current Viewer
  vtkSmartPointer<Slicer> m_currentViewer;

};
#if __GNUC__
#pragma GCC visibility pop
#endif


#endif
