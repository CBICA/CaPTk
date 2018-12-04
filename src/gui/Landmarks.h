/**
\file Landmarks.h

\brief Declaration of Landmarks class
*/


#ifndef _Landmarks_h_
#define _Landmarks_h_


//#include "CAPTk.h"
#include <vector>

class vtkPoints;
class vtkFloatArray;
class vtkPolyData;

enum LANDMARK_TYPE
{
  NONE = -1, DEFAULT, TUMOR_POINTS, TISSUE_POINTS
};

#define LANDMARK_MAX 350
#define LANDMARK_MAX_TUMOR_POINTS 20
#define LANDMARK_MAX_TISSUE_POINTS 300

#define NumberOfTissueTypes 15

static const char labels_laconic[NumberOfTissueTypes][32] = { "BG", "CSF", "GM", "WM", "VS", "ED", "NCR", "TU", "NE", "CB", "VT", "CAN", "CAE", "RTN", "RTE" };
static const char labels_verbose[NumberOfTissueTypes][32] = { "Background", "Cerebrospinal Fluid", "Gray Matter", "White Matter", "Vessels", "Edema", "Necrosis", "Tumor: Enhancing", "Tumor: Non Enhacing", "Cerebellum", "Ventricular CSF", "Cavity: Non Enhancing", "Cavity: Enhancing", "Recurrence Tumor: Non Enhancing", "Recurrence Tumor: Enhancing" };
static const char genericLabels_laconic[NumberOfTissueTypes][32] = { "L01", "L02", "L03", "L04", "L05", "L06", "L07", "L08", "L09", "L10", "L11", "L12", "L13", "L14", "L15" };
static const char genericLabels_verbose[NumberOfTissueTypes][32] = { "Label-01", "Label-02", "Label-03", "Label-04", "Label-05", "Label-06", "Label-07", "Label-08", "Label-09", "Label-10", "Label-11", "Label-12", "Label-13", "Label-14", "Label-15" };

//! Enum for different Tissue Types
enum TissueTypes
{
  BG, CSF, GM, WM, VS, ED, NCR, TU, NE, CB, VT, CAN, CAE, RTN, RTE
};


struct sLandmark
{
  bool bValid;
  float coordinates[3];
  float image_coordinates[3];
  double pixel_value;
  float radius;
  int id;
};


/**
\class Landmarks

\brief Tumor landmarks
*/
class Landmarks
{
public:
  Landmarks(int type);
  ~Landmarks();
  bool AddLandmark(float x, float y, float z, float r, double value, int id);
  bool SetLandmark(int index, float x, float y, float z, float r, double value, int id);
  bool RemoveLastLandmark();
  bool RemoveLandmark(int index);
  unsigned int GetNumberOfPoints() { return static_cast<unsigned int>(mLandmarks.size()); }
  void UpdateData();
  void Clear();

  int mType;
  std::vector<sLandmark> mLandmarks;
  vtkPoints* mPoint;
  vtkFloatArray* mRadius;
  vtkFloatArray* mLandID;
  vtkPolyData* mLandDataWithID;
  vtkPolyData* mLandDataWithRadius;
};


#endif