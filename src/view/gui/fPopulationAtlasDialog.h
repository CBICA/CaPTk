///////////////////////////////////////////////////////////////////////////////////////
// fPopulationAtlasDialog.h
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

#ifndef _fPopulationAtlasDialog_h_
#define _fPopulationAtlasDialog_h_


//#include "CAPTk.h"
#include "ui_fPopulationAtlasDialog.h"

/**
\class fPopulationAtlasDialog

\brief This class controls the elements in the recurrence dialog
*/
class fPopulationAtlasDialog : public QDialog, private Ui::fPopulationAtlasDialog
{
  Q_OBJECT

public:
  fPopulationAtlasDialog();
  ~fPopulationAtlasDialog();
  int mode;

  void SetCurrentImagePath(const QString &inputPath)
  {
	  mInputPathName = inputPath;
  }

  QString mInputPathName;


public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();
  void OpenOutputDirectory();
  void OpenInputFile();
  void OpenInputAtlasFile();

signals:
  void GeneratePopualtionAtlas(std::string inputdirectory, std::string atlasfile, std::string outputdirectory);
};


#endif





