/**
\file Landmarks.cpp

\brief Implementation of Landmarks class
*/

#include "Landmarks.h"
//#include "CAPTk.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"


Landmarks::Landmarks(int type)
{
  mType = type;
  mLandmarks.resize(0);
  //
  mPoint = vtkPoints::New();
  mRadius = vtkFloatArray::New();
  mLandID = vtkFloatArray::New();
  mLandDataWithID = vtkPolyData::New();
  mLandDataWithRadius = vtkPolyData::New();
}

Landmarks::~Landmarks()
{
  if (mPoint) {
    mPoint->Delete();
  //}
  //if (mRadius) {
  //  mRadius->Delete();
  //}
  //if (mLandID) {
  //  mLandID->Delete();
  //}
  //if (mLandDataWithID) {
  //  mLandDataWithID->Delete();
  //}
  //if (mLandDataWithRadius) {
  //  mLandDataWithRadius->Delete();
  //}
  }
  if (mRadius) {
    mRadius->Delete();
  }
  if (mLandID) {
    mLandID->Delete();
  }
  if (mLandDataWithID) {
    mLandDataWithID->Delete();
  }
  if (mLandDataWithRadius) {
    mLandDataWithRadius->Delete();
  }
}

void Landmarks::Clear()
{
  mLandmarks.resize(0);
  //
  mPoint->Reset();
  mRadius->Reset();
  mLandID->Reset();
  mLandDataWithID->Reset();
  mLandDataWithRadius->Reset();
}

bool Landmarks::AddLandmark(float x, float y, float z, float r, double value, int id)
{
  switch (mType)
  {
  case LANDMARK_TYPE::DEFAULT:
    if (mLandmarks.size() >= LANDMARK_MAX_TUMOR_POINTS)
    {
      return false;
    }
    break;
  case LANDMARK_TYPE::TUMOR_POINTS:
    if (mLandmarks.size() >= LANDMARK_MAX_TUMOR_POINTS)
    {
      return false;
    }
    break;
  case LANDMARK_TYPE::TISSUE_POINTS:
    if (mLandmarks.size() >= LANDMARK_MAX_TISSUE_POINTS)
    {
      return false;
    }
    break;
  default:
    break;
  }

  sLandmark point;
  point.bValid = true;
  point.coordinates[0] = x;
  point.coordinates[1] = y;
  point.coordinates[2] = z;
  point.pixel_value = value;
  point.radius = r;
  point.id = id;
  mLandmarks.push_back(point);
  mPoint->InsertNextPoint(x, y, z);
  mRadius->InsertNextTuple1(r);
  mLandID->InsertNextTuple1(id);

  UpdateData();

  return true;
}

bool Landmarks::SetLandmark(int index, float x, float y, float z, float r, double value, int id)
{
  if (mType == LANDMARK_TYPE::TUMOR_POINTS)
  {
    if (index < 0 || index >= (int)mLandmarks.size())
    {
      return false;
    }

    mLandmarks[index].bValid = true;
    mLandmarks[index].coordinates[0] = x;
    mLandmarks[index].coordinates[1] = y;
    mLandmarks[index].coordinates[2] = z;
    mLandmarks[index].pixel_value = value;
    mLandmarks[index].radius = r;
    mLandmarks[index].id = id;

    mPoint->Reset();
    mRadius->Reset();
    mLandID->Reset();
    for (int i = 0; i < (int)mLandmarks.size(); i++)
    {
      if (mLandmarks[i].bValid)
      {
        mPoint->InsertNextPoint(mLandmarks[i].coordinates[0], mLandmarks[i].coordinates[1], mLandmarks[i].coordinates[2]);
        mRadius->InsertNextTuple1(mLandmarks[i].radius);
        mLandID->InsertNextTuple1(mLandmarks[i].id);
      }
    }

    UpdateData();

    return true;
  }
  else if (mType == LANDMARK_TYPE::TISSUE_POINTS)
  {
    if (index < 0 || index >= (int)mLandmarks.size())
    {
      return false;
    }

    mLandmarks[index].bValid = true;
    mLandmarks[index].coordinates[0] = x;
    mLandmarks[index].coordinates[1] = y;
    mLandmarks[index].coordinates[2] = z;
    mLandmarks[index].pixel_value = value;
    mLandmarks[index].radius = r;
    mLandmarks[index].id = id;

    mPoint->Reset();
    mRadius->Reset();
    mLandID->Reset();
    for (int i = 0; i < (int)mLandmarks.size(); i++)
    {
      if (mLandmarks[i].bValid) {
        mPoint->InsertNextPoint(mLandmarks[i].coordinates[0], mLandmarks[i].coordinates[1], mLandmarks[i].coordinates[2]);
        mRadius->InsertNextTuple1(mLandmarks[i].radius);
        mLandID->InsertNextTuple1(mLandmarks[i].id);
      }
    }

    UpdateData();

    return true;
  }

  return false;
}

bool Landmarks::RemoveLastLandmark()
{
  if (mType == LANDMARK_TYPE::DEFAULT || mType == LANDMARK_TYPE::TUMOR_POINTS)
  {
    mPoint->SetNumberOfPoints(mPoint->GetNumberOfPoints() - 1);
    mRadius->RemoveLastTuple();
    mLandID->RemoveLastTuple();
    mLandmarks.pop_back();

    UpdateData();

    return true;
  }

  return false;
}

bool Landmarks::RemoveLandmark(int index)
{
  if (mType == LANDMARK_TYPE::DEFAULT || mType == LANDMARK_TYPE::TUMOR_POINTS) 
  {
    int npoints = mPoint->GetNumberOfPoints();
    for (int i = index; i < npoints - 1; i++) {
      mPoint->InsertPoint(i, mPoint->GetPoint(i + 1));
      mRadius->InsertTuple1(i, mRadius->GetTuple1(i + 1));
    }
    mPoint->SetNumberOfPoints(npoints - 1);
    mRadius->SetNumberOfTuples(npoints - 1);
    mLandID->RemoveLastTuple();
    mLandmarks.erase(mLandmarks.begin() + index);

    UpdateData();

    return true;
  }
  else if (mType == LANDMARK_TYPE::TISSUE_POINTS) 
  {
    int npoints = mPoint->GetNumberOfPoints();
    for (int i = index; i < npoints - 1; i++) {
      mPoint->InsertPoint(i, mPoint->GetPoint(i + 1));
      mRadius->InsertTuple1(i, mRadius->GetTuple1(i + 1));
      mLandID->InsertTuple1(i, mLandID->GetTuple1(i + 1));
    }
    mPoint->SetNumberOfPoints(npoints - 1);
    mRadius->SetNumberOfTuples(npoints - 1);
    mLandID->SetNumberOfTuples(npoints - 1);
    mLandmarks.erase(mLandmarks.begin() + index);

    UpdateData();

    return true;
  }

  return false;
}

void Landmarks::UpdateData()
{
  mLandDataWithID->SetPoints(mPoint);
  mLandDataWithID->GetPointData()->SetScalars(mLandID);
  mLandDataWithID->Modified();
#if VTK_MAJOR_VERSION <= 5
  mLandDataWithID->Update();
#endif
  //
  mLandDataWithRadius->SetPoints(mPoint);
  mLandDataWithRadius->GetPointData()->SetScalars(mRadius);
  mLandDataWithRadius->Modified();
#if VTK_MAJOR_VERSION <= 5
  mLandDataWithRadius->Update();
#endif
}
