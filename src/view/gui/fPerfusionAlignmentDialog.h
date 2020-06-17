///////////////////////////////////////////////////////////////////////////////////////
// fPerfusionAligner.h
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

#ifndef _fPerfusionAlignmentDialog_h_
#define _fPerfusionAlignmentDialog_h_


//#include "CAPTk.h"
#include "ui_fPerfusionAlignmentDialog.h"

/**
\class fPerfusionAligner

\brief This class controls the elements in the recurrence dialog
*/
class fPerfusionAligner : public QDialog, private Ui::fPerfusionAligner
{
  Q_OBJECT

public:
  fPerfusionAligner();
  ~fPerfusionAligner();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName);
  }

  QString mInputPathName;
  QString mInputT1cePathName;
  QString mOutputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void SelectOutputImage();
  void SelectInputImage();
  void SelectT1ceInputImage();

signals:
  void RunPerfusionAlignmentCalculation(double echotime, int before,int after, const std::string inputfile, const std::string inputt1cefile, std::string outputFolder);
};


#endif
