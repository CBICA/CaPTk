///////////////////////////////////////////////////////////////////////////////////////
// fSkullStripper.h
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
// License Agreement: https://www.med.upenn.edu/cbica/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fSkullStripDialog_h_
#define _fSkullStripDialog_h_


//#include "CAPTk.h"
#include "ui_fSkullStripDialog.h"

/**
\class fSkullStripper

\brief This class controls the elements in the recurrence dialog
*/
class fSkullStripper : public QDialog, private Ui::fSkullStripper
{
  Q_OBJECT

public:
  fSkullStripper();
  ~fSkullStripper();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName + "skullStrip_Output.nii.gz");
  }

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void OpenReferenceImage();
  void OpenReferenceMaskImage();
  void OpenInputImage();
  void SelectOutputImage();

signals:
  void RunSkullStripping(const std::string referenceAtlasFile, const std::string referenceAtlasMaskFile, const std::string inputImageFile, const std::string outputImageFile);
};


#endif
