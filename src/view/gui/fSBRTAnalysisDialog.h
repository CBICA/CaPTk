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
// License Agreement: https://www.med.upenn.edu/sbia/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fSBRTAnalysisDialog_h_
#define _fSBRTAnalysisDialog_h_


//#include "CAPTk.h"
#include "ui_fSBRTAnalysisDialog.h"

/**
\class fSBRTAnalysisDialog

\brief This class displays the popup dialog for LungField part of SBRT applications
*/
class fSBRTAnalysisDialog : public QDialog, private Ui::fSBRTAnalysisDialog
{
  Q_OBJECT

public:

 /**
 \brief Constructor
 */
  fSBRTAnalysisDialog();

 /**
 \brief Destructor
 */
  ~fSBRTAnalysisDialog();

 /**
 \brief Set the path of the currently loaded image
 \param inputPath path of currently loaded image as QString
*/
  void SetCurrentImagePath(const QString &inputPath)
  {
	  mInputPathName = inputPath;
  }

 /**
 \brief Set the predicted survival value
 \param val as double
*/
  void SetSurvivalValue(float val);
  
 /**
 \brief Set the predicted nodal failure value
 \param val as double
*/
  void SetNodalFailureValue(float val);
  
 /**
 \brief ivars
 */
  QString mInputPathName;	//! path of currently loaded data

  void SetTrainedModelLink(std::string link)
  {
    m_trainedModelLink = link;
  }
  std::string m_trainedModelLink;

 /*
 \brief slots
 */
public slots:
 /**
 \brief slot callback when 'OK' button on UI is clicked
 */
  void OnOKButtonClicked();

  /**
\brief slot callback when 'download' url on UI is clicked
*/
  void OnURLClicked(const QString&);

  /**
\brief slot callback when 'Browse' button on UI is clicked
*/
  void OnSelectModelBtnClicked();

  /**
  \brief slot callback when 'Cancel' button on UI is clicked
  */
  void OnCancelButtonClicked();

};

#endif





