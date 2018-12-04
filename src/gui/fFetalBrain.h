///////////////////////////////////////////////////////////////////////////////////////
// fRecurrenceDialog.h
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
// License Agreement: https://www.cbica.upenn.edu/sbia/software/license.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fFetalBrain_h_
#define _fFetalBrain_h_


//#include "CAPTk.h"
#include "ui_fFetalBrain.h"
//#include "FetalBrain.h"

#define SUBJECT_CLASSIFICATION 0
#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2


class fFetalBrain : public QDialog, private Ui::fFetalBrain
{
  Q_OBJECT

public:
  fFetalBrain();
  ~fFetalBrain();
  int mode;
  std::map<std::string, float> featurevec;
  void SetCurrentImagePath(const std::string &inputPath)
  {
    mInputPathName = inputPath;
  }

  std::string mInputPathName;

 
  public slots:
  void SkullSegmentation();
  void Prediction();

  void SingleSubject();
  void TrainingRadioButtonChecked();
  
  void OpenSVMModelFile();
  void OpenTestSubjectsDirectory();
  std::string  OpenTrainSubjectDirectory();

  void SelectOutputDirectory();

  void CancelButtonPressed();
  void ConfirmButtonPressed();

  //void Prediction();


signals:
  void SubjectBasedFetalVent();
  void skullstripfun();
  void drawlinear();
  void TrainNewFetalModel(std::string, std::string);
};


#endif





