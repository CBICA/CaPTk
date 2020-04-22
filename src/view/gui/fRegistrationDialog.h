///////////////////////////////////////////////////////////////////////////////////////
// fRegistrationDialog.h
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

#ifndef _fRegistrationDialog_h_
#define _fRegistrationDialog_h_


//#include "CAPTk.h"
#include "ui_fRegistrationDialog.h"

/**
\class fRegistrationDialog

\brief This class controls the elements in the registration dialog
*/
class fRegistrationDialog : public QDialog, private Ui::fRegistrationDialog
{
    Q_OBJECT

public:
    fRegistrationDialog();
    ~fRegistrationDialog();
    bool affineMode = true, rigidMode = true, deformMode = false;
    std::string metric = "NMI";
    std::string radius = "5x5x5";
    std::string m_iterations = "100x50x5";
    bool radii = false;
    int clicked = 0;
    QString matExtn = ".mat";

    QString mInputPathName;

    public slots:
    void CancelButtonPressed();
    void ConfirmButtonPressed();
    void SelectFixedFile();
    void populateFields(std::string & configfilePath, std::string & matfilePath);    
    void SelectMovingFile1();
    void SelectMovingFile2();
    void SelectMovingFile3();
    void SelectMovingFile4();
    void SelectMovingFile5();    
    void SelectMovingOutputFile1();
    void SelectMovingOutputFile2();
    void SelectMovingOutputFile3();
    void SelectMovingOutputFile4();    
    void SelectMovingOutputFile5();
    void SelectMatrixFile1();
    void SelectMatrixFile2();
    void SelectMatrixFile3();
    void SelectMatrixFile4();
    void SelectMatrixFile5();
    void SelectedAffineMode();
    void SelectedRigidMode();
    void SelectedDeformMode();
    void SelectedMetric(int index);
    void addMoreImages();
    void SelectGenerateMatrix(bool checked);
    void ResetButtonPressed();
    void setRadii(const QString nccRadii);
    void getIterations(const QString iterations);
 
signals:
    void RegistrationSignal(std::string fixedfilename, 
      std::vector<std::string> inputfilenames, 
      std::vector<std::string> outputfilenames, 
      std::vector<std::string> matrixfilenames, 
      std::string metrics, 
      bool rigidMode, bool affineMode, bool deformMode, 
      std::string radii, std::string iterations);
};

#endif