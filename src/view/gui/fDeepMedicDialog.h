///////////////////////////////////////////////////////////////////////////////////////
// fDeepMedicDialog.h
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

#ifndef _fDeepMedicDialog_h_
#define _fDeepMedicDialog_h_


//#include "CAPTk.h"
#include "ui_fDeepMedicDialog.h"

/**
\class fDeepMedicDialog

\brief This class controls the elements in the DICOM converter
*/
class fDeepMedicDialog : public QDialog, private Ui::fDeepMedicDialog
{
  Q_OBJECT

public:

  //! Default models available
  enum ModelTypes
  {
    Tumor,
    SkullStripping,
    Custom,
    Max
  };
  fDeepMedicDialog();
  ~fDeepMedicDialog();

  QString m_exe, m_dataDir, m_modelDir; // contains full path and exe name of dcm2nii
  std::string m_baseModelDir;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
  }


  QString mInputPathName;


public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void SelectModelDirectory();
  void SelectOutputDirectory();
  void SetDefaultModel();
  void SetDefaultModel(int modelType);

signals:
  void RunDeepMedic(const std::string modelDirectory, const std::string outputDirectory);
};


#endif
