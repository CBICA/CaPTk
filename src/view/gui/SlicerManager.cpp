/**
\file  SlicerManager.cpp

\brief Implementation of SlicerManager class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "SlicerManager.h"
#include "Slicer.h"
#include "Landmarks.h"
#include "SlicerManagerCommand.h"
#include "itkImageToVTKImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionConstIterator.h"

#include "vtkInteractorStyle.h"
#include "vtkMatrix4x4.h"
#include "vtkRenderWindow.h"
#include "vtkTransform.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty2D.h"
#include "vtkCamera.h"
#include "vtkLookupTable.h"
#include "vtkImageMapToWindowLevelColors.h"

#include <QDirIterator>
#include <QDir>
#include <QVariant>

#include "cbicaITKImageInfo.h"
#include "CaPTkEnums.h"

SlicerManager::SlicerManager(int numberOfSlicers, Landmarks* landmarks, Landmarks* seedPoints, Landmarks* tissuePoints)
{

  mPathFileName = "";
  mFileName = "";
  mBaseFileName = "";
  mId = "";
  mLastError = "";
  mPreset = 0;

  mLinkedId.resize(0);

  mPreviousSlice.resize(numberOfSlicers);
  mPreviousTSlice.resize(numberOfSlicers);

  mOrder = 0;
  mThresholdIndex = 20;
  mThreshold = 0;

  mImage = vtkSmartPointer<vtkImageData>::New();
  mTransform = vtkSmartPointer<vtkTransform>::New();
  for (int i = 0; i < numberOfSlicers; i++) {
    mSlicers.push_back(vtkSmartPointer<Slicer>::New());
  }

  mLandmarks = landmarks;
  mSeedPoints = seedPoints;
  mTissuePoints = tissuePoints;

  mCurrentLandmarksType = LANDMARK_TYPE::NONE;
  mCurrentLandmarksRow = -1;
  mCurrentLandmarksCol = -1;
  mTissueSelectedEntry = -1;

  mImageType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;
  mImageSubType = CAPTK::ImageModalityType::IMAGE_TYPE_UNDEFINED;;
}

SlicerManager::~SlicerManager()
{
}

int SlicerManager::GetDimension()
{
  if (mImage) {
    int dim = 0;
#if VTK_MAJOR_VERSION <= 5
    int* extent = mImage->GetWholeExtent();
#else
    int* extent = mImage->GetExtent();
#endif
    if (extent[4] != extent[5]) {
      dim = 3;
    }
    else if (extent[3] != extent[4]) {
      dim = 2;
    }
    else if (extent[0] != extent[1]) {
      dim = 1;
    }
    return dim;
  }
  else {
    return -1;
  }
}

void SlicerManager::AddLink(const std::string &newId)
{
  mLinkedId.push_back(newId);
  mLinkedId.unique();
}

void SlicerManager::RemoveLink(const std::string &oldId)
{
  mLinkedId.remove(oldId);
}

void SlicerManager::SetPerfImage(ImageTypeFloat4D::Pointer image)
{
  mPerfusionImagePointer = image;
  //crop one 3D image from the 4D series
  ImageTypeFloat4D::RegionType region = image->GetLargestPossibleRegion();
  ImageTypeFloat4D::IndexType regionIndex;
  ImageTypeFloat4D::SizeType regionSize;
  regionSize[0] = region.GetSize()[0];
  regionSize[1] = region.GetSize()[1];
  regionSize[2] = region.GetSize()[2];
  regionSize[3] = 0;
  regionIndex[0] = 0;
  regionIndex[1] = 0;
  regionIndex[2] = 0;
  regionIndex[3] = 0;
  ImageTypeFloat4D::RegionType desiredRegion(regionIndex, regionSize);

  typedef itk::ExtractImageFilter< ImageTypeFloat4D, ImageTypeFloat3D > FilterType;
  FilterType::Pointer filter = FilterType::New();
  filter->SetExtractionRegion(desiredRegion);
  filter->SetInput(mPerfusionImagePointer);

  filter->SetDirectionCollapseToIdentity(); // This is required.
  filter->Update();
  mITKImage = filter->GetOutput();
  UpdateVtkImage();
}

void SlicerManager::SetOriginalOrigin(ImageTypeFloat3D::PointType origin)
{
  std::vector< double > actualOriginVector;
  actualOriginVector.resize(3);
  for (size_t i = 0; i < 3; i++)
  {
    actualOriginVector[i] = origin[i];
  }
  SetOriginalOrigin(actualOriginVector);
}
void SlicerManager::SetOriginalOrigin(std::vector< double > origin)
{
  mOrigin[0] = origin[0];
  mOrigin[1] = origin[1];
  if (origin.size() == 2)
  {
    mOrigin[2] = 0;
  }
  else
  {
    mOrigin[2] = origin[2];
  }
}
void SlicerManager::SetImage(ImageTypeFloat3D::Pointer t1ceimage)
{
  mITKImage = t1ceimage;
  UpdateVtkImage();

}

ImageTypeFloat3D::Pointer SlicerManager::GetITKImage()
{
  return this->mITKImage;
}

void SlicerManager::SetComparisonMode(bool mode)
{
  this->m_ComparisonMode = mode;
  for (int i = 0; i < (int)mSlicers.size(); i++)
  {
    mSlicers[i]->SetComparisonMode(mode);
  }
}

bool SlicerManager::GetComparisonMode()
{
  return this->m_ComparisonMode;
}

void SlicerManager::SetImageOnSlicer(vtkSmartPointer<vtkImageData> image, vtkTransform* transform, int slicer)
{
  mSlicers[slicer]->SetImage(image, transform);
}

void SlicerManager::UpdateVtkImage()
{
  typedef itk::Image<float, 3> OutputImageType3D;
  // Convert from ITK object to VTK object
  typedef itk::ImageToVTKImageFilter<OutputImageType3D> ConverterType;
  ConverterType::Pointer converter = ConverterType::New();
  mItkToVtkConverters = dynamic_cast<itk::ProcessObject *>(converter.GetPointer());
  converter->SetInput(mITKImage);
  converter->Update();
  mImage = converter->GetOutput();

  vtkSmartPointer<vtkMatrix4x4> matrix = vtkSmartPointer<vtkMatrix4x4>::New();
  matrix->Identity();
  // for (unsigned int i = 0; i < mITKImage->GetImageDimension(); i++)
  // {
  //   for (unsigned int j = 0; j < mITKImage->GetImageDimension(); j++)
  //   {
  //     (*matrix)[i][j] = mITKImage->GetDirection()[i][j];
  //     (*matrix)[i][3] -= (*matrix)[i][j] * mITKImage->GetOrigin()[j];
  //   }
  //   (*matrix)[i][3] += mITKImage->GetOrigin()[i];
  // }
  for (unsigned int i = 0; i < mITKImage->GetImageDimension(); i++)
  {
    float temp = 1;
    for (unsigned int j = 0; j < mITKImage->GetImageDimension(); j++)
    {
      //(*matrix)[i][j] = mITKImage->GetDirection()[i][j];
      //(*matrix)[i][3] -= (*matrix)[i][j] * mITKImage->GetOrigin()[j];
      matrix->SetElement(i, j, mITKImage->GetDirection()[i][j]);
      temp -= matrix->GetElement(i, j) * mITKImage->GetOrigin()[j];
    }
    //(*matrix)[i][3] += mITKImage->GetOrigin()[i];
    temp += mITKImage->GetOrigin()[i];
    matrix->SetElement(i, 3, temp);
  }

  matrix->Invert();
  mTransform->SetMatrix(matrix);
  for (int i = 0; i < (int)mSlicers.size(); i++)
  {
    mSlicers[i]->SetImage(mImage, mTransform);
    mSlicers[i]->GetImageActor()->SetInterpolate(0);
  }
}
void SlicerManager::SetFilename(const std::string &filename)
{
  mPathFileName = filename;

  std::string path_out, base_out, ext_out;
  cbica::splitFileName(filename, path_out, base_out, ext_out);

  mFileName = base_out + ext_out;
  mBaseFileName = base_out;

}

void SlicerManager::SetMask( vtkSmartPointer< vtkImageData > mask)
{
  mMask = mask;
}

Slicer* SlicerManager::GetSlicer(int i)
{
  return mSlicers[i];
}

void SlicerManager::UpdateSlicer(int num, bool state)
{
  if (mSlicers[num]->GetImage()) {
    mSlicers[num]->SetDisplayMode(state);
  }
}

void SlicerManager::SetSlicerWindow(int i, vtkRenderWindow* RW)
{
  mSlicers[i]->SetRenderWindow(i, RW);
  mSlicers[i]->SetImageSeriesDescription(this->mBaseFileName);
}

void SlicerManager::SetInteractorStyleNavigator(int i, vtkInteractorStyle* style)
{
  SlicerManagerCommand *smc = SlicerManagerCommand::New();
  smc->SM = this;
  smc->SetSlicerNumber(i);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->SetInteractorStyle(style);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::KeyPressEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::WindowLevelEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::EndWindowLevelEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::StartWindowLevelEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::PickEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::StartPickEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::LeaveEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::UserEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::MouseWheelForwardEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::MouseWheelBackwardEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::LeftButtonReleaseEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::EndPickEvent, smc);
  mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::EndInteractionEvent, smc);
  //mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->AddObserver(vtkCommand::RightButtonPressEvent, smc);

  if (this->mImageSubType == CAPTK::ImageModalityType::IMAGE_TYPE_T1CE )
  {
    mSlicers[i]->borderWidget->SetInteractor(mSlicers[i]->GetRenderWindow()->GetInteractor()->GetInteractorStyle()->GetInteractor());
    static_cast<vtkBorderRepresentation*>(mSlicers[i]->borderWidget->GetRepresentation())->GetBorderProperty()->SetColor(1, 0, 0);
    mSlicers[i]->borderWidget->SelectableOff();
    mSlicers[i]->borderCallback->SetRenderer(mSlicers[i]->GetRenderer());
    mSlicers[i]->borderCallback->SetImageActor(mSlicers[i]->GetImageActor());
    mSlicers[i]->borderWidget->AddObserver(vtkCommand::InteractionEvent, smc/*mSlicers[i]->borderCallback*/);
    mSlicers[i]->borderWidget->On();
  }
  smc->Delete();
}

void SlicerManager::LeftButtonReleaseEvent(int slicer)
{
  emit LeftButtonReleaseSignal(slicer);
}

void SlicerManager::SetSliceOrientation(int slicer, int orientation)
{
  mSlicers[slicer]->SetSliceOrientation(orientation);
  emit UpdateOrientation(slicer, orientation);
}

void SlicerManager::ToggleInterpolation()
{
  bool interpolate = !(mSlicers[0]->GetImageActor()->GetInterpolate());
  for (int i = 0; i < (int)mSlicers.size(); i++) 
  {
    mSlicers[i]->GetImageActor()->SetInterpolate(interpolate);
  }
}

void SlicerManager::SetColorWindow(double s)
{
  for (int i = 0; i < (int)mSlicers.size(); i++) 
  {
    mSlicers[i]->SetColorWindow(s);
  }
}

void SlicerManager::SetColorLevel(double s)
{
  for (int i = 0; i < (int)mSlicers.size(); i++)
  {
    mSlicers[i]->SetColorLevel(s);
  }
}

void SlicerManager::SetOpacity(int i, double factor)
{
  mSlicers[i]->SetOpacity(1 / factor);
}

void SlicerManager::UpdateViews(int slicer)
{
  double x = (mSlicers[slicer]->GetCurrentPosition()[0] - mSlicers[slicer]->GetInput()->GetOrigin()[0]) / mSlicers[slicer]->GetInput()->GetSpacing()[0];
  double y = (mSlicers[slicer]->GetCurrentPosition()[1] - mSlicers[slicer]->GetInput()->GetOrigin()[1]) / mSlicers[slicer]->GetInput()->GetSpacing()[1];
  double z = (mSlicers[slicer]->GetCurrentPosition()[2] - mSlicers[slicer]->GetInput()->GetOrigin()[2]) / mSlicers[slicer]->GetInput()->GetSpacing()[2];

#if VTK_MAJOR_VERSION <= 5
  if (x >= mSlicers[slicer]->GetInput()->GetWholeExtent()[0] && x <= mSlicers[slicer]->GetInput()->GetWholeExtent()[1] &&
    y >= mSlicers[slicer]->GetInput()->GetWholeExtent()[2] && y <= mSlicers[slicer]->GetInput()->GetWholeExtent()[3] &&
    z >= mSlicers[slicer]->GetInput()->GetWholeExtent()[4] && z <= mSlicers[slicer]->GetInput()->GetWholeExtent()[5])
#else
  if (x >= mSlicers[slicer]->GetInput()->GetExtent()[0] && x <= mSlicers[slicer]->GetInput()->GetExtent()[1] &&
    y >= mSlicers[slicer]->GetInput()->GetExtent()[2] && y <= mSlicers[slicer]->GetInput()->GetExtent()[3] &&
    z >= mSlicers[slicer]->GetInput()->GetExtent()[4] && z <= mSlicers[slicer]->GetInput()->GetExtent()[5])
#endif
  {
    mSlicers[slicer]->SetActive(true);
    //
    mSlicers[slicer]->UpdateCursorPosition();

    switch (mSlicers[slicer]->GetSliceOrientation()) 
    {
    case vtkImageViewer2::SLICE_ORIENTATION_XY:
      if (mSlicers[slicer]->GetSlice() != (int)ROUND(z))
      {
        mSlicers[slicer]->SetSlice((int)ROUND(z));
      }
      break;
    case vtkImageViewer2::SLICE_ORIENTATION_XZ:
      if (mSlicers[slicer]->GetSlice() != (int)ROUND(y)) 
      {
        mSlicers[slicer]->SetSlice((int)ROUND(y));
      }
      break;
    case vtkImageViewer2::SLICE_ORIENTATION_YZ:
      if (mSlicers[slicer]->GetSlice() != (int)ROUND(x)) 
      {
        mSlicers[slicer]->SetSlice((int)ROUND(x));
      }
      break;
    }
    mSlicers[slicer]->Render();

    for (int i = 0; i < (int)mSlicers.size(); i++) 
    {
      if (i != slicer) {
        mSlicers[i]->SetActive(false);
      }
      //
      bool bDraw = mSlicers[i]->GetRenderer()->GetDraw();
      int RW0 = mSlicers[i]->GetRenderWindow()->GetSize()[0];
      int RW1 = mSlicers[i]->GetRenderWindow()->GetSize()[1];
      //
      if (i != slicer && bDraw && RW0 > 2 && RW1 > 2)
      {
        mSlicers[i]->SetCurrentPosition(mSlicers[slicer]->GetCurrentPosition()[0], mSlicers[slicer]->GetCurrentPosition()[1], mSlicers[slicer]->GetCurrentPosition()[2]);
        mSlicers[i]->UpdateCursorPosition();

        switch (mSlicers[i]->GetSliceOrientation()) 
        {
        case vtkImageViewer2::SLICE_ORIENTATION_XY:
          if (mSlicers[i]->GetSlice() != (int)ROUND(z)) 
          {
            mSlicers[i]->SetSlice((int)ROUND(z));
          }
          break;
        case vtkImageViewer2::SLICE_ORIENTATION_XZ:
          if (mSlicers[i]->GetSlice() != (int)ROUND(y)) 
          {
            mSlicers[i]->SetSlice((int)ROUND(y));
          }
          break;
        case vtkImageViewer2::SLICE_ORIENTATION_YZ:
          if (mSlicers[i]->GetSlice() != (int)ROUND(x)) 
          {
            mSlicers[i]->SetSlice((int)ROUND(x));
          }
          break;
        }

        mSlicers[i]->Render();
      }

      UpdateSlice(i);
    }
  }
}

void SlicerManager::UpdateLinked(int slicer)
{
  double x = (mSlicers[slicer]->GetCurrentPosition()[0] - mSlicers[slicer]->GetInput()->GetOrigin()[0]) / mSlicers[slicer]->GetInput()->GetSpacing()[0];
  double y = (mSlicers[slicer]->GetCurrentPosition()[1] - mSlicers[slicer]->GetInput()->GetOrigin()[1]) / mSlicers[slicer]->GetInput()->GetSpacing()[1];
  double z = (mSlicers[slicer]->GetCurrentPosition()[2] - mSlicers[slicer]->GetInput()->GetOrigin()[2]) / mSlicers[slicer]->GetInput()->GetSpacing()[2];

#if VTK_MAJOR_VERSION <= 5
  if (x >= mSlicers[slicer]->GetInput()->GetWholeExtent()[0] && x <= mSlicers[slicer]->GetInput()->GetWholeExtent()[1] &&
    y >= mSlicers[slicer]->GetInput()->GetWholeExtent()[2] && y <= mSlicers[slicer]->GetInput()->GetWholeExtent()[3] &&
    z >= mSlicers[slicer]->GetInput()->GetWholeExtent()[4] && z <= mSlicers[slicer]->GetInput()->GetWholeExtent()[5])
#else
  if (x >= mSlicers[slicer]->GetInput()->GetExtent()[0] && x <= mSlicers[slicer]->GetInput()->GetExtent()[1] &&
    y >= mSlicers[slicer]->GetInput()->GetExtent()[2] && y <= mSlicers[slicer]->GetInput()->GetExtent()[3] &&
    z >= mSlicers[slicer]->GetInput()->GetExtent()[4] && z <= mSlicers[slicer]->GetInput()->GetExtent()[5])
#endif
  {
    for (std::list<std::string>::const_iterator i = mLinkedId.begin(); i != mLinkedId.end(); i++) 
    {
      emit UpdateLinkManager(*i, slicer, mSlicers[slicer]->GetCurrentPosition()[0], mSlicers[slicer]->GetCurrentPosition()[1], mSlicers[slicer]->GetCurrentPosition()[2]);
    }
  }
}

void SlicerManager::updateToRefCam(Slicer *refSlicer)
{
  vtkCamera *refCam = refSlicer->GetRenderer()->GetActiveCamera();
  double refPosition[3], refFocal[3];
  refCam->GetPosition(refPosition);
  refCam->GetFocalPoint(refFocal);

  for (int i = 0; i < (int)mSlicers.size(); i++) {
    vtkCamera *camera = mSlicers[i]->GetRenderer()->GetActiveCamera();
    camera->SetParallelScale(refCam->GetParallelScale());

    double position[3], focal[3];
    camera->GetPosition(position);
    camera->GetFocalPoint(focal);

    if (refSlicer->GetSliceOrientation() == mSlicers[i]->GetSliceOrientation()) 
    {
      for (int i = 0; i<3; i++) {
        position[i] = refPosition[i];
        focal[i] = refFocal[i];
      }
    }

    if (refSlicer->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_XY) 
    {
      if (mSlicers[i]->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_XZ) 
      {
        position[0] = refPosition[0];
        focal[0] = refFocal[0];
      }
      if (mSlicers[i]->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_YZ) 
      {
        position[1] = refPosition[1];
        focal[1] = refFocal[1];
      }
    }

    if (refSlicer->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_XZ) 
    {
      if (mSlicers[i]->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_YZ) 
      {
        position[2] = refPosition[2];
        focal[2] = refFocal[2];
      }
      if (mSlicers[i]->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_XY) 
      {
        position[0] = refPosition[0];
        focal[0] = refFocal[0];
      }
    }

    if (refSlicer->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_YZ) 
    {
      if (mSlicers[i]->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_XY) 
      {
        position[1] = refPosition[1];
        focal[1] = refFocal[1];
      }
      if (mSlicers[i]->GetSliceOrientation() == vtkImageViewer2::SLICE_ORIENTATION_XZ) 
      {
        position[2] = refPosition[2];
        focal[2] = refFocal[2];
      }
    }

    camera->SetFocalPoint(focal);
    camera->SetPosition(position);
    mSlicers[i]->ForceUpdateDisplayExtent();
  }

}

double SlicerManager::GetColorWindow()
{
  if (mSlicers.size()) {
    return mSlicers[0]->GetColorWindow();
  }
  return -1;
}

double SlicerManager::GetColorLevel()
{
  if (mSlicers.size()) {
    return mSlicers[0]->GetColorLevel();
  }
  return -1;
}

void SlicerManager::ResetTransformationToIdentity()
{
  this->mTransform->Identity();
  for (int i = 0; i < this->GetNumberOfSlicers(); i++)
  {
    this->GetSlicer(i)->ForceUpdateDisplayExtent();
    this->GetSlicer(i)->ResetCamera();
    this->GetSlicer(i)->Render();
  }
}

void SlicerManager::Render()
{
  for (int i = 0; i < (int)mSlicers.size(); i++) 
  {
    mSlicers[i]->Render();
  }
}

void SlicerManager::GenerateDefaultLookupTable()
{
  SetPreset(mPreset);
}

void SlicerManager::RemoveActors()
{
  for (int i = 0; i < (int)mSlicers.size(); i++) 
  {
    mSlicers[i]->SetDisplayMode(0);
    mSlicers[i]->GetRenderer()->RemoveActor(mSlicers[i]->GetImageActor());
  }
}

void SlicerManager::UpdateInfoOnCursorPosition(int slicer)
{
  double x = mSlicers[slicer]->GetCursorPosition()[0];
  double y = mSlicers[slicer]->GetCursorPosition()[1];
  double z = mSlicers[slicer]->GetCursorPosition()[2];
  double X = (x - mSlicers[slicer]->GetInput()->GetOrigin()[0]) / mSlicers[slicer]->GetInput()->GetSpacing()[0];
  double Y = (y - mSlicers[slicer]->GetInput()->GetOrigin()[1]) / mSlicers[slicer]->GetInput()->GetSpacing()[1];
  double Z = (z - mSlicers[slicer]->GetInput()->GetOrigin()[2]) / mSlicers[slicer]->GetInput()->GetSpacing()[2];
  //
  // round up pixel values
  X = ROUND(X);
  Y = ROUND(Y);
  Z = ROUND(Z);
  x = X * mSlicers[slicer]->GetInput()->GetSpacing()[0] + mSlicers[slicer]->GetInput()->GetOrigin()[0];
  y = Y * mSlicers[slicer]->GetInput()->GetSpacing()[1] + mSlicers[slicer]->GetInput()->GetOrigin()[1];
  z = Z * mSlicers[slicer]->GetInput()->GetSpacing()[2] + mSlicers[slicer]->GetInput()->GetOrigin()[2];
  //
  double value = -VTK_DOUBLE_MAX;
#if VTK_MAJOR_VERSION <= 5
  if (X >= mSlicers[slicer]->GetInput()->GetWholeExtent()[0] && X <= mSlicers[slicer]->GetInput()->GetWholeExtent()[1] &&
    Y >= mSlicers[slicer]->GetInput()->GetWholeExtent()[2] && Y <= mSlicers[slicer]->GetInput()->GetWholeExtent()[3] &&
    Z >= mSlicers[slicer]->GetInput()->GetWholeExtent()[4] && Z <= mSlicers[slicer]->GetInput()->GetWholeExtent()[5])
#else
  if (X >= mSlicers[slicer]->GetInput()->GetExtent()[0] && X <= mSlicers[slicer]->GetInput()->GetExtent()[1] &&
    Y >= mSlicers[slicer]->GetInput()->GetExtent()[2] && Y <= mSlicers[slicer]->GetInput()->GetExtent()[3] &&
    Z >= mSlicers[slicer]->GetInput()->GetExtent()[4] && Z <= mSlicers[slicer]->GetInput()->GetExtent()[5])
#endif
  {
    //	mSlicers[slicer]->GetInput()
    value = this->GetScalarComponentAsDouble(mSlicers[slicer]->GetInput(), X, Y, Z);

    emit UpdatePosition(mSlicers[slicer]->GetCursorVisibility(), x, y, z, X, Y, Z, value);
  }
}

void SlicerManager::Activated()
{
  emit currentImageChanged(mId);
}

void SlicerManager::Picked()
{
  emit currentPickedImageChanged(mId);
}

void SlicerManager::UpdateWindowLevel()
{
  emit WindowLevelChanged();
}

void SlicerManager::UpdateSlice(int slicer)
{
  if (mPreviousSlice[slicer] == mSlicers[slicer]->GetSlice()) 
  {
    return;
  }
  emit UpdateSlice(slicer, mSlicers[slicer]->GetSlice());
  mSlicers[slicer]->Render();
  mPreviousSlice[slicer] = mSlicers[slicer]->GetSlice();
}

void SlicerManager::UpdateSliceRange(int slicer)
{
  emit UpdateSliceRange(slicer, mSlicers[slicer]->GetSliceRange()[0], mSlicers[slicer]->GetSliceRange()[1]);
}

void SlicerManager::SetPreset(int preset)
{
  double window = mSlicers[0]->GetColorWindow();
  double level = mSlicers[0]->GetColorLevel();
  double range[2];
  vtkLookupTable* LUT = NULL;
  int i;

  switch (preset) {
  case PRESET_AUTO:
    mImage->GetScalarRange(range);
    window = range[1] - range[0];
    level = (range[1] + range[0]) * 0.5;
    //
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->GetWindowLevel()->SetLookupTable(NULL);
    }
    this->SetColorWindow(window);
    this->SetColorLevel(level);
    break;
  case PRESET_USER:
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->GetWindowLevel()->SetLookupTable(NULL);
    }
    this->SetColorWindow(window);
    this->SetColorLevel(level);
    break;
  case PRESET_LABEL:
    window = 1;
    level = -1;
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetColorLevel(level);
      mSlicers[i]->SetColorWindow(window);
    }
    //
    LUT = vtkLookupTable::New();
    LUT->SetTableRange(0, 250);
    LUT->Build();
    //
    for (i = 0; i < LUT->GetNumberOfTableValues(); i++) 
    {
      LUT->SetTableValue(i, 0, 0, 0);
    }
    // CSF
    LUT->SetTableValue(LUT->GetIndex(10), 1, 1, 0);
    // VS or VEIN or VEINS
    LUT->SetTableValue(LUT->GetIndex(25), 1, 1, 0);
    // VT
    LUT->SetTableValue(LUT->GetIndex(50), 1, 1, 0);
    // GM
    LUT->SetTableValue(LUT->GetIndex(150), 0.5, 0.5, 0.5);
    // WM
    LUT->SetTableValue(LUT->GetIndex(250), 1, 1, 1);
    // ED
    LUT->SetTableValue(LUT->GetIndex(100), 0, 1, 1);
    // NCR or CAN
    LUT->SetTableValue(LUT->GetIndex(175), 1, 0, 0);
    // TU or CAE
    LUT->SetTableValue(LUT->GetIndex(200), 0, 0.5, 0.75);
    // NE
    LUT->SetTableValue(LUT->GetIndex(185), 0.75, 0.75, 0.75);
    //LUT->SetTableValue(LUT->GetIndex(185), 0.75, 0.5, 0.75);
    //LUT->SetTableValue(LUT->GetIndex(185), 1, 0.5, 1);
    // CB
    LUT->SetTableValue(LUT->GetIndex(5), 1, 1, 1);
    // RTN
    LUT->SetTableValue(LUT->GetIndex(210), 0.80, 0.80, 0.80);
    // RTE
    LUT->SetTableValue(LUT->GetIndex(220), 0.85, 0.85, 0.85);
    //
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->GetWindowLevel()->SetLookupTable(LUT);
    }
    LUT->Delete();
    break;
  case PRESET_LABEL2:
    window = 1;
    level = -1;
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetColorLevel(level);
      mSlicers[i]->SetColorWindow(window);
    }
    //
    LUT = vtkLookupTable::New();
    LUT->SetTableRange(0, 4);
    LUT->Build();
    //
    for (i = 0; i < LUT->GetNumberOfTableValues(); i++)
    {
      LUT->SetTableValue(i, 0, 0, 0);
    }
    // ED
    LUT->SetTableValue(LUT->GetIndex(2), 0, 1, 1);
    // TU
    LUT->SetTableValue(LUT->GetIndex(4), 0, 0.5, 0.75);
    // NCR
    LUT->SetTableValue(LUT->GetIndex(1), 1, 0, 0);
    // NE
    LUT->SetTableValue(LUT->GetIndex(3), 0.75, 0.75, 0.75);
    //
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->GetWindowLevel()->SetLookupTable(LUT);
    }
    LUT->Delete();
    break;
  case PRESET_THRESHOLD:
    window = 1;
    level = -1;
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetColorLevel(level);
      mSlicers[i]->SetColorWindow(window);
    }
    //
    mImage->GetScalarRange(range);
    range[1] *= 0.4;
    //
    LUT = vtkLookupTable::New();
    LUT->SetTableRange(range[0], range[1]);
    LUT->SetNumberOfTableValues(256);
    LUT->Build();
    //
    for (i = 0; i < LUT->GetNumberOfTableValues(); i++)
    {
      LUT->SetTableValue(i, 0, 0, 0);
    }
    //
    mThreshold = range[0] + (range[1] - range[0]) * mThresholdIndex / (LUT->GetNumberOfTableValues() - 1);
    //
    for (i = 0; i <= mThresholdIndex; i++)
    {
      LUT->SetTableValue(i, 0, 0.5, 0.75);
    }
    //
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->GetWindowLevel()->SetLookupTable(LUT);
    }
    LUT->Delete();
    break;
  case PRESET_PROB:
    window = 1;
    level = -1;
    for (i = 0; i < (int)mSlicers.size(); i++) 
    {
      mSlicers[i]->SetColorLevel(level);
      mSlicers[i]->SetColorWindow(window);
    }
    //
    mImage->GetScalarRange(range);
    //
    LUT = vtkLookupTable::New();
    LUT->SetTableRange(range[0], range[1]);
    LUT->SetNumberOfTableValues(256);
    LUT->Build();
    //
    //   0:   0,   0, 128
    //  32:   0,   0, 254
    //  96:   0, 254, 254
    // 160: 254, 254,   0
    // 224: 254,   0,   0
    // 255: 128,   0,   0
    LUT->SetTableValue(0, 0, 0, 0);
    for (i = 1; i < 32; i++)
    {
      LUT->SetTableValue(i, 0, 0, (128.0 + 126.0*i / 32.0) / 255.0);
    }
    for (i = 0; i < 64; i++)
    {
      LUT->SetTableValue(32 + i, 0, (254.0*i / 64.0) / 255.0, 254.0 / 255.0);
    }
    for (i = 0; i < 64; i++)
    {
      LUT->SetTableValue(96 + i, (254.0*i / 64.0) / 255.0, 254.0 / 255.0, (254.0 - 254.0*i / 64.0) / 255.0);
    }
    for (i = 0; i < 64; i++)
    {
      LUT->SetTableValue(160 + i, 254.0 / 255.0, (254.0 - 254.0*i / 64.0) / 255.0, 0);
    }
    for (i = 0; i < 32; i++)
    {
      LUT->SetTableValue(224 + i, (254.0 - 126.0*i / 31.0) / 255.0, 0, 0);
    }
    //
    for (i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->GetWindowLevel()->SetLookupTable(LUT);
    }
    LUT->Delete();
    break;
  }

  mPreset = preset;
}

void SlicerManager::AddLandmark(float x, float y, float z)
{
  double x_index = (x - mOrigin[0] /*mSlicers[0]->GetInput()->GetOrigin()[0]*/) / mSlicers[0]->GetInput()->GetSpacing()[0];
  double y_index = (y - mOrigin[1] /*mSlicers[0]->GetInput()->GetOrigin()[1]*/) / mSlicers[0]->GetInput()->GetSpacing()[1];
  double z_index = (z - mOrigin[2] /*mSlicers[0]->GetInput()->GetOrigin()[2]*/) / mSlicers[0]->GetInput()->GetSpacing()[2];
  //
  // round up pixel values
  x_index = std::abs(std::round(x_index));
  y_index = std::abs(std::round(y_index));
  z_index = std::abs(std::round(z_index));
  x = x_index * mSlicers[0]->GetInput()->GetSpacing()[0] + mOrigin[0] /*mSlicers[0]->GetInput()->GetOrigin()[0]*/;
  y = y_index * mSlicers[0]->GetInput()->GetSpacing()[1] + mOrigin[1] /*mSlicers[0]->GetInput()->GetOrigin()[1]*/;
  z = z_index * mSlicers[0]->GetInput()->GetSpacing()[2] + mOrigin[2] /*mSlicers[0]->GetInput()->GetOrigin()[2]*/;
  //
#if VTK_MAJOR_VERSION <= 5
  if (x_index >= mSlicers[0]->GetInput()->GetWholeExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetWholeExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetWholeExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetWholeExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetWholeExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetWholeExtent()[5])
#else
  if (x_index >= mSlicers[0]->GetInput()->GetExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetExtent()[5])
#endif
  {
    double value = 0 /*this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index)*/;

    if (mCurrentLandmarksType == LANDMARK_TYPE::DEFAULT)
    {
      int id = mLandmarks->GetNumberOfPoints();
      if (mLandmarks->AddLandmark(x, y, z, 0, value, id))
      {
        emit LandmarkAdded();
      }
    }
    else if ((mCurrentLandmarksType == LANDMARK_TYPE::TUMOR_POINTS) && (mSeedPoints->GetNumberOfPoints() > 0))
    {
      if ((mCurrentLandmarksRow < 0) || mCurrentLandmarksRow >= (int)mSeedPoints->GetNumberOfPoints())
      {
        int id = mSeedPoints->GetNumberOfPoints();
        if (mSeedPoints->AddLandmark(x, y, z, 0, value, id))
        {
          emit SeedPointsAdded();
        }
      }
      else
      {
        double r = mSeedPoints->mLandmarks[mCurrentLandmarksRow].radius;
        if (mSeedPoints->SetLandmark(mCurrentLandmarksRow, x, y, z, r, value, mCurrentLandmarksRow))
        {
          emit SeedPointsAdded(mCurrentLandmarksRow, true);
        }
      }
    }
    else if ((mCurrentLandmarksType == LANDMARK_TYPE::TISSUE_POINTS) && (mTissueSelectedEntry != -1))
    {
      emit TissuePointsRemoved(mTissueSelectedEntry);
      // Functionality taken from Shift+Space
      for (unsigned int j = 0; j < mTissuePoints->GetNumberOfPoints(); j++)
      {
        if (mTissuePoints->mLandmarks[j].coordinates[0] == x &&
          mTissuePoints->mLandmarks[j].coordinates[1] == y &&
          mTissuePoints->mLandmarks[j].coordinates[2] == z)
          return;
      }

      double value = 0 /*this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index)*/;
      if (mTissuePoints->AddLandmark(x, y, z, 0, value, mCurrentLandmarksRow + 1))
      {
        emit TissuePointsAdded(mCurrentLandmarksRow);

        VectorDouble tisseTypeCounter;
        tisseTypeCounter.resize(15);
        for (unsigned int i = 0; i < mTissuePoints->GetNumberOfPoints(); i++)
          tisseTypeCounter[mTissuePoints->mLandmarks[i].id - 1] = tisseTypeCounter[mTissuePoints->mLandmarks[i].id - 1] + 1;
        emit UpdateNumberOfPoints(mTissuePoints->GetNumberOfPoints(), tisseTypeCounter);
      }
      emit HighlightTableIndexOnDeletion(x, y, z, x, y, z, value);
    }
  }
}

//void SlicerManager::AddLoadedLandmark(double x, double y, double z, int landmarktype)
//{
//	double x_index = x;
//	double y_index = y;
//	double z_index = z;
//	//
//	// round up pixel values
//	x_index = ROUND(x_index);
//	y_index = ROUND(y_index);
//	z_index = ROUND(z_index);
//	x = x_index * mSlicers[0]->GetInput()->GetSpacing()[0] + mSlicers[0]->GetInput()->GetOrigin()[0];
//	y = y_index * mSlicers[0]->GetInput()->GetSpacing()[1] + mSlicers[0]->GetInput()->GetOrigin()[1];
//	z = z_index * mSlicers[0]->GetInput()->GetSpacing()[2] + mSlicers[0]->GetInput()->GetOrigin()[2];
//	//
//#if VTK_MAJOR_VERSION <= 5
//	if (x_index >= mSlicers[0]->GetInput()->GetWholeExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetWholeExtent()[1] &&
//		y_index >= mSlicers[0]->GetInput()->GetWholeExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetWholeExtent()[3] &&
//		z_index >= mSlicers[0]->GetInput()->GetWholeExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetWholeExtent()[5])
//#else
//	if (x_index >= mSlicers[0]->GetInput()->GetExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetExtent()[1] &&
//		y_index >= mSlicers[0]->GetInput()->GetExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetExtent()[3] &&
//		z_index >= mSlicers[0]->GetInput()->GetExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetExtent()[5])
//#endif
//	{
//		double value = this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index);
//		int id = landmarktype;
//		if (mDrawingPoints->AddLandmark(x, y, z, 0, value, id))
//		{
//			int nearpoints = 0;
//			int farpoints = 0;
//			for (int i = 0; i < mDrawingPoints->GetNumberOfPoints(); i++)
//			{
//				if (mDrawingPoints->mLandmarks[i].id == NEAR_POINT_LABEL)
//					nearpoints++;
//				else if (mDrawingPoints->mLandmarks[i].id == FAR_POINT_LABEL)
//					farpoints++;
//			}
//			emit UpdateNumberOfPoints(nearpoints, farpoints);
//		}
//	}
//}


void SlicerManager::AddLandmarkRadius(float x, float y, float z)
{
  double x_index = (x - mOrigin[0] /*mSlicers[0]->GetInput()->GetOrigin()[0]*/) / mSlicers[0]->GetInput()->GetSpacing()[0];
  double y_index = (y - mOrigin[1] /*mSlicers[0]->GetInput()->GetOrigin()[1]*/) / mSlicers[0]->GetInput()->GetSpacing()[1];
  double z_index = (z - mOrigin[2] /*mSlicers[0]->GetInput()->GetOrigin()[2]*/) / mSlicers[0]->GetInput()->GetSpacing()[2];
  //
  // round up pixel values
  x_index = std::abs(std::round(x_index));
  y_index = std::abs(std::round(y_index));
  z_index = std::abs(std::round(z_index));
  // TBD: the next 3 lines are redundant
  x = x_index * mSlicers[0]->GetInput()->GetSpacing()[0] + mOrigin[0] /*mSlicers[0]->GetInput()->GetOrigin()[0]*/;
  y = y_index * mSlicers[0]->GetInput()->GetSpacing()[1] + mOrigin[1] /*mSlicers[0]->GetInput()->GetOrigin()[1]*/;
  z = z_index * mSlicers[0]->GetInput()->GetSpacing()[2] + mOrigin[2] /*mSlicers[0]->GetInput()->GetOrigin()[2]*/;

  auto currentExtent = mSlicers[0]->GetInput()->GetExtent();
  //
#if VTK_MAJOR_VERSION <= 5
  if (x_index >= mSlicers[0]->GetInput()->GetWholeExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetWholeExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetWholeExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetWholeExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetWholeExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetWholeExtent()[5])
#else
  if (x_index >= currentExtent[0] && x_index <= currentExtent[1] &&
    y_index >= currentExtent[2] && y_index <= currentExtent[3] &&
    z_index >= currentExtent[4] && z_index <= currentExtent[5])
#endif
  {
    if (mCurrentLandmarksType == LANDMARK_TYPE::DEFAULT)
    {
    }
    else if ((mCurrentLandmarksType == LANDMARK_TYPE::TUMOR_POINTS) && (mSeedPoints->GetNumberOfPoints() > 0))
    {

      double r, x0, y0, z0, value = 0;
      if (mSeedPoints->mLandmarks[mCurrentLandmarksRow].bValid)
      {
        x0 = mSeedPoints->mLandmarks[mCurrentLandmarksRow].coordinates[0];
        y0 = mSeedPoints->mLandmarks[mCurrentLandmarksRow].coordinates[1];
        z0 = mSeedPoints->mLandmarks[mCurrentLandmarksRow].coordinates[2];
        //value = mSeedPoints->mLandmarks[mCurrentLandmarksRow].pixel_value;
        //
        r = sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0));
        //
        if (mSeedPoints->SetLandmark(mCurrentLandmarksRow, x0, y0, z0, r, value, mCurrentLandmarksRow))
        {
          emit SeedPointsAdded(mCurrentLandmarksRow, true);
        }
      }
    }
    else if (mCurrentLandmarksType == LANDMARK_TYPE::TISSUE_POINTS)
    {
      for (unsigned int j = 0; j < mTissuePoints->GetNumberOfPoints(); j++)
      {
        if (mTissuePoints->mLandmarks[j].coordinates[0] == x &&
          mTissuePoints->mLandmarks[j].coordinates[1] == y &&
          mTissuePoints->mLandmarks[j].coordinates[2] == z)
          return;
      }

      double value = 0 /*this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index)*/;
      if (mTissuePoints->AddLandmark(x, y, z, 0, value, mCurrentLandmarksRow/* + 1*/))
      {
        emit TissuePointsAdded(mCurrentLandmarksRow);

        VectorDouble tisseTypeCounter;
        tisseTypeCounter.resize(NumberOfTissueTypes);
        for (unsigned int i = 0; i < mTissuePoints->GetNumberOfPoints(); i++)
          tisseTypeCounter[mTissuePoints->mLandmarks[i].id - 1] = tisseTypeCounter[mTissuePoints->mLandmarks[i].id - 1] + 1;
        emit UpdateNumberOfPoints(mTissuePoints->GetNumberOfPoints(), tisseTypeCounter);
      }
    }
  }
}

void SlicerManager::NextImageWithOrder(int order)
{
  emit ChangeImageWithOrder(this, order);
}

void SlicerManager::VerticalSliderHasChanged(int slicer, int slice)
{
  emit AVerticalSliderHasChanged(slicer, slice);
}

double SlicerManager::GetScalarComponentAsDouble(vtkSmartPointer< vtkImageData > image, double X, double Y, double Z, int component)
{
  int ix, iy, iz;
  return mSlicers[0]->GetScalarComponentAsDouble(image, X, Y, Z, ix, iy, iz, component);
}

void SlicerManager::SetCurrentLandmarksType(int type, int row, int col)
{

  if (type == LANDMARK_TYPE::DEFAULT)
  {
    for (int i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetLandmarks(mLandmarks, LANDMARK_TYPE::DEFAULT);
    }
  }
  else if (type == LANDMARK_TYPE::TUMOR_POINTS)
  {
    for (int i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetLandmarks(mSeedPoints, LANDMARK_TYPE::TUMOR_POINTS);
    }
  }
  else if (type == LANDMARK_TYPE::TISSUE_POINTS)
  {
    for (int i = 0; i < (int)mSlicers.size(); i++)
    {
      mSlicers[i]->SetLandmarks(mTissuePoints, LANDMARK_TYPE::TISSUE_POINTS);
    }
  }
  mCurrentLandmarksType = type;
  mCurrentLandmarksRow = row;
  mCurrentLandmarksCol = col;
}

void SlicerManager::RemoveLandmark(float x, float y, float z)
{
  double x_index = ROUND((x - mSlicers[0]->GetInput()->GetOrigin()[0]) / mSlicers[0]->GetInput()->GetSpacing()[0]);
  double y_index = ROUND((y - mSlicers[0]->GetInput()->GetOrigin()[1]) / mSlicers[0]->GetInput()->GetSpacing()[1]);
  double z_index = ROUND((z - mSlicers[0]->GetInput()->GetOrigin()[2]) / mSlicers[0]->GetInput()->GetSpacing()[2]);


  x = x_index * mSlicers[0]->GetInput()->GetSpacing()[0] + mSlicers[0]->GetInput()->GetOrigin()[0];
  y = y_index * mSlicers[0]->GetInput()->GetSpacing()[1] + mSlicers[0]->GetInput()->GetOrigin()[1];
  z = z_index * mSlicers[0]->GetInput()->GetSpacing()[2] + mSlicers[0]->GetInput()->GetOrigin()[2];
  //
#if VTK_MAJOR_VERSION <= 5
  if (x_index >= mSlicers[0]->GetInput()->GetWholeExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetWholeExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetWholeExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetWholeExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetWholeExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetWholeExtent()[5])
#else
  if (x_index >= mSlicers[0]->GetInput()->GetExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetExtent()[5])
#endif
  {
    // double value = this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index);

    /*if (mCurrentLandmarksType == LANDMARK_TYPE::DEFAULT)
    {
    if (mLandmarks->RemoveLandmark(x, y, z))
    {
    emit LandmarkAdded();
    }
    }
    else*/ if (mCurrentLandmarksType == LANDMARK_TYPE::TUMOR_POINTS)
    {
      //if (mNearPoints->RemoveLandmark(x, y, z))
      //{
      //	//emit UpdateNumberOfNearPoints(mNearPoints->GetNumberOfPoints()-1);
      //	//emit SeedPointsAdded();
      //}
    }
    else if (mCurrentLandmarksType == LANDMARK_TYPE::TISSUE_POINTS)
    {
      //if (mDrawingPoints->RemoveLandmark(x, y, z))
      //{
      //	int nearpoints = 0;
      //	int farpoints = 0;
      //	for (int i = 0; i < mDrawingPoints->GetNumberOfPoints(); i++)
      //	{
      //		if (mDrawingPoints->mLandmarks[i].id == NEAR_POINT_LABEL)
      //			nearpoints++;
      //		else if (mDrawingPoints->mLandmarks[i].id == FAR_POINT_LABEL)
      //			farpoints++;
      //	}
      //	emit UpdateNumberOfPoints(nearpoints, farpoints);
      //}

    }
  }
}
void SlicerManager::ActionAdded(std::vector<PointVal>& points)
{
  QVariantList list;
  for (size_t i = 0; i < points.size(); i++)
  {
    list.append(qVariantFromValue(points[i]));
  }
  emit UpdateActionInMain(list);
}
void SlicerManager::UpdateBorderCoordinates(double startZ, double endZ)
{
  emit UpdateBorderWidgetInMain(startZ, endZ);
}
void SlicerManager::UpdateBorderCoordinates(double startX, double startY, double endX, double endY)
{
  emit UpdateBorderWidgetInMain(startX, startY, endX, endY);
}
void SlicerManager::AddLandmarkShift(float x, float y, float z)
{
  double x_index = (x - mOrigin[0] /*mSlicers[0]->GetInput()->GetOrigin()[0]*/) / mSlicers[0]->GetInput()->GetSpacing()[0];
  double y_index = (y - mOrigin[1] /*mSlicers[0]->GetInput()->GetOrigin()[1]*/) / mSlicers[0]->GetInput()->GetSpacing()[1];
  double z_index = (z - mOrigin[2] /*mSlicers[0]->GetInput()->GetOrigin()[2]*/) / mSlicers[0]->GetInput()->GetSpacing()[2];
  //
  // round up pixel values
  x_index = std::abs(std::round(x_index));
  y_index = std::abs(std::round(y_index));
  z_index = std::abs(std::round(z_index));
  x = x_index * mSlicers[0]->GetInput()->GetSpacing()[0] + mOrigin[0] /*mSlicers[0]->GetInput()->GetOrigin()[0]*/;
  y = y_index * mSlicers[0]->GetInput()->GetSpacing()[1] + mOrigin[1] /*mSlicers[0]->GetInput()->GetOrigin()[1]*/;
  z = z_index * mSlicers[0]->GetInput()->GetSpacing()[2] + mOrigin[2] /*mSlicers[0]->GetInput()->GetOrigin()[2]*/;
  //
#if VTK_MAJOR_VERSION <= 5
  if (x_index >= mSlicers[0]->GetInput()->GetWholeExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetWholeExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetWholeExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetWholeExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetWholeExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetWholeExtent()[5])
#else
  if (x_index >= mSlicers[0]->GetInput()->GetExtent()[0] && x_index <= mSlicers[0]->GetInput()->GetExtent()[1] &&
    y_index >= mSlicers[0]->GetInput()->GetExtent()[2] && y_index <= mSlicers[0]->GetInput()->GetExtent()[3] &&
    z_index >= mSlicers[0]->GetInput()->GetExtent()[4] && z_index <= mSlicers[0]->GetInput()->GetExtent()[5])
#endif
  {
    double value = 0 /*this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index)*/;

    if (mCurrentLandmarksType == LANDMARK_TYPE::DEFAULT)
    {
      int id = mLandmarks->GetNumberOfPoints();
      if (mLandmarks->AddLandmark(x, y, z, 0, value, id))
      {
        emit LandmarkAdded();
      }
    }
    else if (mCurrentLandmarksType == LANDMARK_TYPE::TUMOR_POINTS)
    {
      for (unsigned int j = 0; j < mSeedPoints->GetNumberOfPoints(); j++)
      {
        if (mSeedPoints->mLandmarks[j].coordinates[0] == x &&
          mSeedPoints->mLandmarks[j].coordinates[1] == y &&
          mSeedPoints->mLandmarks[j].coordinates[2] == z)
          return;
      }
      int id = mSeedPoints->GetNumberOfPoints();
      if (mSeedPoints->AddLandmark(x, y, z, 0, value, id))
      {
        emit SeedPointsAdded();
      }
      mCurrentLandmarksRow = id;
    }
    else if (mCurrentLandmarksType == LANDMARK_TYPE::TISSUE_POINTS)
    {
      for (unsigned int j = 0; j < mTissuePoints->GetNumberOfPoints(); j++)
      {
        if (mTissuePoints->mLandmarks[j].coordinates[0] == x &&
          mTissuePoints->mLandmarks[j].coordinates[1] == y &&
          mTissuePoints->mLandmarks[j].coordinates[2] == z)
          return;
      }

      double value = 0 /*this->GetScalarComponentAsDouble(mSlicers[0]->GetInput(), x_index, y_index, z_index)*/;
      if (mTissuePoints->AddLandmark(x, y, z, 0, value, mCurrentLandmarksRow/* + 1*/))
      {
        emit TissuePointsAdded(mCurrentLandmarksRow);

        VectorDouble tisseTypeCounter;
        tisseTypeCounter.resize(15);
        for (unsigned int i = 0; i < mTissuePoints->GetNumberOfPoints(); i++)
          tisseTypeCounter[mTissuePoints->mLandmarks[i].id - 1] = tisseTypeCounter[mTissuePoints->mLandmarks[i].id - 1] + 1;
        emit UpdateNumberOfPoints(mTissuePoints->GetNumberOfPoints(), tisseTypeCounter);
      }
    }
  }
}
void SlicerManager::Get3DImageAtCurrentPerfusionIndex(int sliderindex)
{
  typedef ImageTypeFloat3D ImageType3D;
  typedef ImageTypeFloat4D ImageType4D;
  ImageType4D::SizeType cropSize;
  cropSize[0] = mPerfusionImagePointer->GetLargestPossibleRegion().GetSize()[0];
  cropSize[1] = mPerfusionImagePointer->GetLargestPossibleRegion().GetSize()[1];
  cropSize[2] = mPerfusionImagePointer->GetLargestPossibleRegion().GetSize()[2];
  cropSize[3] = 0;

  ImageType4D::IndexType cropIndex;
  cropIndex[0] = 0;
  cropIndex[1] = 0;
  cropIndex[2] = 0;
  cropIndex[3] = sliderindex;

  ImageType4D::RegionType regionToCrop;
  regionToCrop.SetSize(cropSize);
  regionToCrop.SetIndex(cropIndex);

  typedef itk::ExtractImageFilter< ImageType4D, ImageType3D > ExtractionFilterType;
  ExtractionFilterType::Pointer croppingFilter = ExtractionFilterType::New();
  croppingFilter->SetInput(mPerfusionImagePointer);
  croppingFilter->SetExtractionRegion(regionToCrop);
  croppingFilter->SetDirectionCollapseToSubmatrix(); // This is required.
  croppingFilter->Update();
  mITKImage->Graft(croppingFilter->GetOutput());//TBD this is crashing 
 
  UpdateVtkImage();
  this->Render();
}

