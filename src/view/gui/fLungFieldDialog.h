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
// License Agreement: https://www.med.upenn.edu/sbia/software/license.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fLungFieldDialog_h_
#define _fLungFieldDialog_h_


#include "CAPTk.h"
#include "ui_fLungFieldDialog.h"

#define SUBJECT_CLASSIFICATION 0
#define EXISTING_CLASSIFICATION 1
#define TRAIN_MODEL 2

/**
\class fLungFieldDialog

\brief This class displays the popup dialog for LungField part of SBRT applications
*/
class fLungFieldDialog : public QDialog, private Ui::fLungFieldDialog
{
  Q_OBJECT

public:

 /**
 \brief Constructor
 */
  fLungFieldDialog();

 /**
 \brief Destructor
 */
  ~fLungFieldDialog();

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
 \brief slot callback when 'Browse' button for Mask on UI is clicked
 */
  void OnMaskFileBrowseButtonClicked();

 /**
 \brief slot callback when 'Browse' button for output directory is clicked
 */
  void OnOutputDirectoryBrowseButtonClicked();

 /*
 \brief signals
 */
signals:
/**
 \brief signal emitted when 'OK' button is clicked and when required parameters
 are for lung field computation are ready
 */
	void LungFieldParamReady(std::string maskFileName, std::string outputdirectory);
};

#endif





