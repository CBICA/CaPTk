///////////////////////////////////////////////////////////////////////////////////////
// fPerfusionEstimator.h
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
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fPerfusionMeasuresDialog_h_
#define _fPerfusionMeasuresDialog_h_


//#include "CAPTk.h"
#include "ui_fPerfusionMeasuresDialog.h"

/**
\class fPerfusionEstimator

\brief This class controls the elements in the recurrence dialog
*/
class fPerfusionEstimator : public QDialog, private Ui::fPerfusionEstimator
{
  Q_OBJECT

public:
  fPerfusionEstimator();
  ~fPerfusionEstimator();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName);
  }

  QString mInputPathName;
  QString mOutputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void SelectOutputImage();
  void SelectInputImage();

signals:
  void RunPerfusionMeasuresCalculation(const bool rcbv, const bool psr,const bool ph, const std::string inputfile,  std::string outputFolder);
};


#endif
