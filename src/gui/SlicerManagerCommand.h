#ifndef _SlicerManagerCommand_h_
#define _SlicerManagerCommand_h_


#include "CAPTk.h"
#include "SlicerManager.h"
#include "InteractorStyleNavigator.h"
#include "fMainWindow.h"

#include "itkConfidenceConnectedImageFilter.h"
#include "itkVTKImageToImageFilter.h"

#if __GNUC__
#pragma GCC visibility push(hidden)
#endif
class SlicerManagerCommand : public vtkCommand
{
public:
  static SlicerManagerCommand *New()
  {
    return new SlicerManagerCommand;
  }

  void Execute(vtkObject *caller, unsigned long event, void *vtkNotUsed(callData));

  SlicerManager *SM;
  void Dolly(double factor, vtkRenderWindowInteractor *interactor);
  void SetSlicerNumber(int slicer) { mSlicerNumber = slicer; }
  void AddActions();
  void moveCursor(int VisibleInWindow, double x, double y)
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
  //TBD move these functions to utility class Bresenham Line Drawing Algorithm and others
  static std::pair<int, int > point3Dto2D(const PointVal& pt3D, const int orientation)
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
  static PointVal point2Dto3D(const std::pair<int, int >& pt, const int orientation, const int slice, const int value = 0)
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
  static std::vector< std::pair<int, int> > points3Dto2D(const std::vector<PointVal>& points3D, const int orientation)
  {
    std::vector< std::pair<int, int> > points2D;
    for (size_t i = 0; i < points3D.size(); i++)
    {
      points2D.push_back(point3Dto2D(points3D[i], orientation));
    }
    return points2D;
  }
  static std::vector<PointVal> points2Dto3D(const std::vector< std::pair<int, int> >& points2D, const int orientation, const int slice, const int value = 0)
  {
    std::vector<PointVal> points3D;
    for (size_t i = 0; i < points2D.size(); i++)
    {
      points3D.push_back(point2Dto3D(points2D[i], orientation, slice, value));
    }
    return points3D;
  }
  static std::vector< std::pair<int, int> > getCirclePoints(int x, int y, int r)
  {
    if (r <= 0) r = 1;

    std::vector< std::pair<int, int> > points;
    double dtheta = 2 * M_PI / 8 / r;
    int n = 2 * M_PI / dtheta;
    points.push_back(std::make_pair(x + r, y));
    for (int i = 1; i <= n; i++)
    {
      double theta = i*dtheta;
      int x1 = int(x + r*cos(theta) + 0.5);
      int y1 = int(y + r*sin(theta) + 0.5);
      points.push_back(std::make_pair(x1, y1));
    }
    points.push_back(points[0]);
    return points;
  }
  static std::vector< std::pair<int, int> > getLinePoints(int x0, int y0, int x1, int y1, float wd)
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
  static std::vector<PointVal> drawLine(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int width)
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

  static std::vector<PointVal> fillShape(PointVal centerPt, int labelToFill, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice)
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

  static std::vector<PointVal> drawRectangle(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int thickness)
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
  static std::vector<PointVal> drawCircle(PointVal startPt, PointVal endPt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, int width)
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
  static PointVal drawPoint(PointVal pt, vtkSmartPointer<vtkImageData> image, const int orientation, const int slice, const int i = 0, const int j = 0)
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
  void makeStroke(int VisibleInWindow, double x, double y)
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
protected:
  SlicerManagerCommand();
  ~SlicerManagerCommand() {}

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


};
#if __GNUC__
#pragma GCC visibility pop
#endif


#endif
