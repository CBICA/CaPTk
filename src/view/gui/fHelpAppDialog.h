///////////////////////////////////////////////////////////////////////////////////////
// fHelpAppDialog.h
//
// Copyright (c) 2016. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/cbica/captk/license.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fHelpAppDialog_h_
#define _fHelpAppDialog_h_


#include "CAPTk.h"
#include "ui_fHelpAppDialog.h"


/**
\class fHelpAppDialog

\brief This class controls the elements in the help dialog
*/
class fHelpAppDialog : public QDialog, private Ui::fHelpAppDialog
{
  Q_OBJECT

public:
  fHelpAppDialog();
  ~fHelpAppDialog();

  void SetLoadingTab();
  void SetTumorTab();
  void SetShortcutsTab();
  void SetDrawingTab();


  public slots:

public:

};
#endif
