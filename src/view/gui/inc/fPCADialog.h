///////////////////////////////////////////////////////////////////////////////////////
// fPCAEstimator.h
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

#ifndef _fPCADialog_h_
#define _fPCADialog_h_


//#include "CAPTk.h"
#include "ui_fPCADialog.h"

/**
\class fPCAEstimator

\brief This class controls the elements in the recurrence dialog
*/
class fPCAEstimator : public QDialog, private Ui::fPCAEstimator
{
  Q_OBJECT

public:
  fPCAEstimator();
  ~fPCAEstimator();

  QString mInputPathName;
  QString mOutputPathName;
public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenInputImage();
  void SelectOutputImage();

signals:
  void RunPCAEstimation(const int numberOfPCAs, const std::string outputFolder);
};


#endif
