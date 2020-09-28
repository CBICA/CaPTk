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

#ifndef _fSBRTNoduleDialog_h_
#define _fSBRTNoduleDialog_h_


//#include "CAPTk.h"
#include "ui_fSBRTNoduleDialog.h"

/**
\class fSBRTNoduleDialog

\brief This class displays the popup dialog for LungField part of SBRT applications
*/
class fSBRTNoduleDialog : public QDialog, private Ui::fSBRTNoduleDialog
{
  Q_OBJECT

public:

 /**
 \brief Constructor
 */
  fSBRTNoduleDialog();

 /**
 \brief Destructor
 */
  ~fSBRTNoduleDialog();

 /**
 \brief Set the path of the currently loaded image
 \param inputPath path of currently loaded image as QString
*/
  void SetCurrentImagePath(const QString &inputPath)
  {
	  mInputPathName = inputPath;
  }

 /**
 \brief ivars
 */
  QString mInputPathName;	//! path of currently loaded data

 /*
 \brief slots
 */
public slots:
 /**
 \brief slot callback when 'OK' button on UI is clicked
 */
  void OnOKButtonClicked();

 /**
 \brief slot callback when 'Cancel' button on UI is clicked
 */
  void OnCancelButtonClicked();

 /**
 \brief slot callback when 'Browse' button for Seed Image on UI is clicked
 */
  void OnSeedImageBrowseButtonClicked();

 /**
 \brief slot callback when 'Browse' button for output directory is clicked
 */
  //void OnOutputDirectoryBrowseButtonClicked();

 /*
 \brief signals
 */
signals:
/**
 \brief signal emitted when 'OK' button is clicked and when required parameters
 are for lung field computation are ready
 */
	void SBRTNoduleParamReady(std::string seedImage, int labelValue);
};

#endif





