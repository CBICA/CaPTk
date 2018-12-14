/********************************************************************************
** Form generated from reading UI file 'fMainWindow.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FMAINWINDOW_H
#define UI_FMAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include <QtWidgets/QDockWidget>
#include <QSpacerItem>
#include <QGroupBox>
#include <QRadioButton>
#include <QtWidgets/QSpacerItem>
#include "fTumorPanel.h"
#include "fImagesPanel.h"
#include "fDrawingPanel.h"
#include "fTrainingDialog.h"
#include "fFeaturePanel.h"
#include "fRecurrenceDialog.h"
#include "fPseudoProgressionDialog.h"
#include "fRegistrationDialog.h"
#include "fPreprocessingDialog.h"
#include "fSurvivalDialog.h"
#include "fEGFRvIIIDialog.h"
#include "fSkullStripDialog.h"
#include "fPerfusionMeasuresDialog.h"
#include "fDiffusionMeasuresDialog.h"
#include "fPCADialog.h"
#include "fHistoMatchDialog.h"
#include "fWhiteStripeDialog.h"
#include "fDirectionalityDialog.h"
#include "fPopulationAtlasDialog.h"
#include "fImagingSubtypeDialog.h"
#include "fMolecularSubtypeDialog.h"
#include "fDCM2NIfTI.h"
#include "fDeepMedicDialog.h"
#include "fDeepMedicNormDialog.h"
#include "fFetalBrain.h"
#include "fSBRTNoduleDialog.h"
#include "fSBRTAnalysisDialog.h"

#include "QVTKOpenGLWidget.h"
#include "vtkGenericOpenGLRenderWindow.h"
#include "fBottomImageInfoTip.h"

QT_BEGIN_NAMESPACE

class Ui_fMainWindow
{
private:

public:
  fRecurrenceDialog			recurrencePanel;
  fPseudoProgressionDialog pseudoPanel;
  fPopulationAtlasDialog	atlasPanel;
  fRegistrationDialog		registrationPanel;
  fPreprocessingDialog	preprocessingPanel;
  fSurvivalPredictor survivalPanel;
  fEGFRvIIIPredictor egfrv3Panel;
  fMolecularSubtypePredictor msubtypePanel;
  fImagingSubtypePredictor isubtypePanel;
  fFetalBrain fetalbrainpanel;
  fSBRTNoduleDialog nodulePanel;
  fSBRTAnalysisDialog analysisPanel;

  fSkullStripper skullStrippingPanel;
  fPCAEstimator pcaPanel;
  fTrainingSimulator trainingPanel;
  fPerfusionEstimator perfmeasuresPanel;
  fDiffusionEstimator diffmeasuresPanel;
  fDCM2NIfTIConverter dcmConverter;
  fDeepMedicDialog deepMedicDialog;
  fHistoMatcher histoMatchPanel;
  fDeepMedicNormalizer deepMedicNormPanel;
  fWhiteStripeObj whiteStripeNormalizer;
  fDirectionalityDialog directionalityEstimator;

  fDrawingPanel *drawingPanel;
  fFeaturePanel *featurePanel;
  fImagesPanel *imagesPanel;
  fBottomImageInfoTip *infoPanel;
  fTumorPanel *tumorPanel;

  QSlider *image4DSlider;
  QGroupBox *preferencesGroupBox;
  QWidget *centralwidget;

  QWidget *AxialWidget;
  QWidget *SaggitalWidget;
  QWidget *CoronalWidget;

  QVTKWidget *SaggitalViewWidget;
  QVTKWidget *AxialViewWidget;
  QVTKWidget *CoronalViewWidget;

  QSlider *SaggitalViewSlider;
  QSlider *AxialViewSlider;
  QSlider *CoronalViewSlider;



  QWidget *infoWidget;
  QGridLayout *gridLayout_5;
  QTabWidget *m_tabWidget;
  QDockWidget *m_toolTabdock;

  //preference group box members
  QGridLayout *PrefGridLayout;
  QLabel *thresholdLabel;
  QDoubleSpinBox *levelSpinBox;
  QLabel *presetLabel;
  QDoubleSpinBox *windowSpinBox;
  QLabel *windowLabel;
  QDoubleSpinBox *thresholdSpinBox;
  QComboBox *presetComboBox;
  QLabel *levelLabel;


  QStatusBar *statusbar;
  QGridLayout *overallGridLayout;




  //-------------menu-----------
  QMenuBar *menubar;
  QMenu* menuFile;
  QMenu* menuLoadFile;
  QMenu* menuSaveFile;
  QMenu* menuExit;
  QMenu* menuLoadFileDicom;
  QMenu* menuLoadFileNifti;
  QMenu* menuDownload;

  QMenu* menuApp;
  QMenu* menuPreprocessing;
  QMenu* menuHelp;

  QAction *help_discussion;
  QAction *help_download;
  QAction *help_forum;
  QAction *help_bugs;
  QAction *help_features;
  //-------------actions-------------

  QAction *actionLoad_Recurrence_Images;
  QAction *actionLoad_Nifti_Images;
  QAction *actionLoad_Nifti_ROI;


  QAction *actionSave_Nifti_Images;
  QAction *actionSave_Dicom_Images;
  QAction *actionSave_ROI_Images;
  QAction *actionSave_ROI_Dicom_Images;

  QAction *actionHelp_Interactions;
  QAction *actionSave_Images;
  QAction *actionAbout;
  QAction *actionExit;

  QAction *actionAppEGFR;
  QAction *actionAppRecurrence;
  QAction *actionAppGeodesic;

  // obtain list from CMake variables using populateStringListInMenu() function
  std::vector< std::string >
    m_nativeApps, // native CPP applications
    m_preprocessApps, // native pre-processing routines
    m_pyCLIApps, // python command line applications
    m_pyGUIApps; // python graphical applications

  std::map< std::string, std::string > m_allNonNativeApps;

  void setupUi(QMainWindow *fMainWindow)
  {
    actionLoad_Recurrence_Images = new QAction(fMainWindow);
    actionLoad_Nifti_Images = new QAction(fMainWindow);
    actionLoad_Nifti_ROI = new QAction(fMainWindow);
    actionSave_Nifti_Images = new QAction(fMainWindow);
    actionSave_Dicom_Images = new QAction(fMainWindow);
    actionSave_ROI_Images = new QAction(fMainWindow);
    actionSave_ROI_Dicom_Images = new QAction(fMainWindow);

    actionExit = new QAction(fMainWindow);

    actionAppEGFR = new QAction(fMainWindow);
    actionAppRecurrence = new QAction(fMainWindow);
    actionAppGeodesic = new QAction(fMainWindow);

    centralwidget = new QWidget(fMainWindow);
    //--------------------main vertical layout--------------------------
    overallGridLayout = new QGridLayout(centralwidget);
    overallGridLayout->setSpacing(2);
    overallGridLayout->setContentsMargins(2, 2, 2, 2);
    overallGridLayout->setSizeConstraint(QLayout::SetMinimumSize);

    //TopLeveloverallGridLayout: item 1: view widget
    SaggitalWidget = new QWidget(centralwidget);

    QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy2.setHorizontalStretch(20);
    sizePolicy2.setVerticalStretch(20);
    sizePolicy2.setHeightForWidth(SaggitalWidget->sizePolicy().hasHeightForWidth());
    SaggitalWidget->setSizePolicy(sizePolicy2);

    QGridLayout * SaggitalWidgetGridLayout = new QGridLayout(SaggitalWidget);
    SaggitalWidgetGridLayout->setContentsMargins(0, 0, 0, 0);
    SaggitalWidgetGridLayout->setSizeConstraint(QLayout::SetMinimumSize);
    SaggitalWidgetGridLayout->setHorizontalSpacing(0);
    SaggitalViewWidget = new QVTKWidget(SaggitalWidget);
    SaggitalViewWidget->setMouseTracking(true);
    SaggitalWidgetGridLayout->addWidget(SaggitalViewWidget, 0, 0, 1, 1);
    SaggitalViewSlider = new QSlider(SaggitalWidget);
    SaggitalViewSlider->setOrientation(Qt::Vertical);
    SaggitalWidgetGridLayout->addWidget(SaggitalViewSlider, 0, 1, 1, 1);
    overallGridLayout->addWidget(SaggitalWidget, 2, 2, 1, 1);

    AxialWidget = new QWidget(centralwidget);
    QSizePolicy sizePolicy6(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy6.setHorizontalStretch(20);
    sizePolicy6.setVerticalStretch(20);
    sizePolicy6.setHeightForWidth(AxialWidget->sizePolicy().hasHeightForWidth());
    AxialWidget->setSizePolicy(sizePolicy6);
    AxialWidget->setMinimumSize(QSize(256, 256));//TBD not working

    QGridLayout * AxialWidgetGridLayout = new QGridLayout(AxialWidget);
    AxialWidgetGridLayout->setContentsMargins(0, 0, 0, 0);
    AxialWidgetGridLayout->setSizeConstraint(QLayout::SetMinimumSize);
    AxialWidgetGridLayout->setHorizontalSpacing(0);

    AxialViewWidget = new QVTKWidget(AxialWidget);
    AxialViewWidget->setMouseTracking(true);
    AxialWidgetGridLayout->addWidget(AxialViewWidget, 0, 0, 1, 1);
    AxialViewSlider = new QSlider(AxialWidget);
    AxialViewSlider->setOrientation(Qt::Vertical);
    AxialWidgetGridLayout->addWidget(AxialViewSlider, 0, 1, 1, 1);
    overallGridLayout->addWidget(AxialWidget, 2, 1, 1, 1);

    CoronalWidget = new QWidget(centralwidget);
    QSizePolicy sizePolicy7(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy7.setHorizontalStretch(20);
    sizePolicy7.setVerticalStretch(20);
    sizePolicy7.setHeightForWidth(CoronalWidget->sizePolicy().hasHeightForWidth());
    CoronalWidget->setSizePolicy(sizePolicy7);
    QGridLayout * CoronalWidgetGridLayout = new QGridLayout(CoronalWidget);
    CoronalWidgetGridLayout->setContentsMargins(0, 0, 0, 0);
    CoronalWidgetGridLayout->setSizeConstraint(QLayout::SetMinimumSize);
    CoronalWidgetGridLayout->setHorizontalSpacing(0);
    CoronalViewWidget = new QVTKWidget(CoronalWidget);
    CoronalViewWidget->setMouseTracking(true);
    CoronalWidgetGridLayout->addWidget(CoronalViewWidget, 0, 0, 1, 1);
    CoronalViewSlider = new QSlider(CoronalWidget);
    CoronalViewSlider->setOrientation(Qt::Vertical);
    CoronalWidgetGridLayout->addWidget(CoronalViewSlider, 0, 1, 1, 1);

    overallGridLayout->addWidget(CoronalWidget, 2, 0, 1, 1);
    //-------------------------------second line of fMainWindow---------------------------
    image4DSlider = new QSlider(centralwidget);
    image4DSlider->setMaximum(45);
    image4DSlider->setPageStep(1);
    image4DSlider->setSliderPosition(5);
    image4DSlider->setOrientation(Qt::Horizontal);
    overallGridLayout->addWidget(image4DSlider, 1, 0, 1, 3);
    //---------------------------------------------------------------------------------------
    preferencesGroupBox = new QGroupBox(centralwidget);
    QHBoxLayout* bottomLayout = new QHBoxLayout();


    infoPanel = new fBottomImageInfoTip(centralwidget);


    bottomLayout->addWidget(infoPanel);
    bottomLayout->addStretch();
    bottomLayout->addWidget(preferencesGroupBox);
    overallGridLayout->addLayout(bottomLayout, 4, 0, 2, 3);
    //-------------------preferences related objects----------------------

    PrefGridLayout = new QGridLayout(preferencesGroupBox);

    thresholdLabel = new QLabel(preferencesGroupBox);
    levelSpinBox = new QDoubleSpinBox(preferencesGroupBox);
    levelSpinBox->setDecimals(3);
    levelSpinBox->setMinimum(-66000);
    levelSpinBox->setMaximum(66000);
    levelSpinBox->setSingleStep(10);



    presetLabel = new QLabel(preferencesGroupBox);
    windowSpinBox = new QDoubleSpinBox(preferencesGroupBox);
    windowSpinBox->setDecimals(3);
    windowSpinBox->setMinimum(-66000);
    windowSpinBox->setMaximum(66000);
    windowSpinBox->setSingleStep(10);


    windowLabel = new QLabel(preferencesGroupBox);
    thresholdSpinBox = new QDoubleSpinBox(preferencesGroupBox);
    thresholdSpinBox->setDecimals(3);
    thresholdSpinBox->setMinimum(-99999);
    thresholdSpinBox->setMaximum(99999);
    thresholdSpinBox->setSingleStep(1);


    presetComboBox = new QComboBox(preferencesGroupBox);
    QSizePolicy sizePolicy8(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy8.setHorizontalStretch(0);
    sizePolicy8.setVerticalStretch(0);
    sizePolicy8.setHeightForWidth(presetComboBox->sizePolicy().hasHeightForWidth());
    presetComboBox->setSizePolicy(sizePolicy8);
    presetComboBox->setMaximumSize(QSize(100, 16777215));
    levelLabel = new QLabel(preferencesGroupBox);
    levelLabel->setAlignment(Qt::AlignCenter);
    thresholdLabel->setAlignment(Qt::AlignCenter);


    PrefGridLayout->addWidget(windowLabel, 0, 0, 1, 1);
    PrefGridLayout->addWidget(windowSpinBox, 0, 1, 1, 1);

    PrefGridLayout->addWidget(levelLabel, 0, 2, 1, 1);
    PrefGridLayout->addWidget(levelSpinBox, 0, 3, 1, 1);

    PrefGridLayout->addWidget(presetLabel, 1, 0, 1, 1);
    PrefGridLayout->addWidget(presetComboBox, 1, 1, 1, 1);


    PrefGridLayout->addWidget(thresholdLabel, 1, 2, 1, 1);
    PrefGridLayout->addWidget(thresholdSpinBox, 1, 3, 1, 1);

    m_tabWidget = new QTabWidget();


    //------------------------------final recurrence related material------------------------

    QSizePolicy sizePolicy5(QSizePolicy::Preferred, QSizePolicy::Expanding);
    sizePolicy5.setHorizontalStretch(0);
    sizePolicy5.setVerticalStretch(0);

    imagesPanel = new fImagesPanel(); // New Images Panel
    m_tabWidget->addTab(imagesPanel, QString());
    tumorPanel = new fTumorPanel();
    m_tabWidget->addTab(tumorPanel, QString());
    drawingPanel = new fDrawingPanel();
    featurePanel = new fFeaturePanel();
    m_tabWidget->addTab(drawingPanel, QString());
    m_tabWidget->addTab(featurePanel, "Feature Extraction");
    int minheight = /*std::max(drawingPanel->sizeHint().height(), featurePanel->sizeHint().height())*/featurePanel->sizeHint().height() + 25;
    m_tabWidget->setMinimumHeight(minheight);
    m_tabWidget->setMaximumHeight(m_tabWidget->minimumHeight());
    m_toolTabdock = new QDockWidget(fMainWindow);
    m_toolTabdock->setWindowFlags(Qt::Window);

#ifdef Q_OS_WIN
    m_toolTabdock->setFeatures(QDockWidget::DockWidgetFloatable);
#else
    //TBD fix this - work around untill solved
    m_toolTabdock->setFeatures(QDockWidget::NoDockWidgetFeatures);
#endif
    m_toolTabdock->setWidget(m_tabWidget);
    overallGridLayout->addWidget(m_toolTabdock, 0, 0, 1, 3);

    QFrame * frame = new QFrame(fMainWindow);
    sizePolicy5.setHeightForWidth(frame->sizePolicy().hasHeightForWidth());
    frame->setSizePolicy(sizePolicy5);
    frame->setFrameShape(QFrame::HLine);
    frame->setFrameShadow(QFrame::Sunken);

    overallGridLayout->addWidget(frame, 3, 0, 1, 3);


    fMainWindow->setCentralWidget(centralwidget);
    AxialViewWidget->raise();
    CoronalViewWidget->raise();
    SaggitalViewWidget->raise();
    infoPanel->raise();
    m_tabWidget->raise();


    //---------------setting menu and status bar for the main window---------------
    statusbar = new QStatusBar(fMainWindow);
    fMainWindow->setStatusBar(statusbar);

    menubar = new QMenuBar(fMainWindow);
    menuFile = new QMenu("File");
    menuLoadFile = new QMenu("Load");
    menuSaveFile = new QMenu("Save");
    menuExit = new QMenu("Exit");

    menuLoadFileDicom = new QMenu("Dicom");
    menuLoadFileNifti = new QMenu("Nifti");
    menuFile->addMenu(menuLoadFile);
    menuFile->addMenu(menuSaveFile);

    menuApp = new QMenu("Applications");

    menuPreprocessing = new QMenu("Preprocessing");

    actionHelp_Interactions = new QAction(fMainWindow);

    actionAbout = new QAction(fMainWindow);

    menuHelp = new QMenu("Help");
    menuHelp->addAction(actionHelp_Interactions);
    menuDownload = menuHelp->addMenu("Sample Data");
    auto supportMenu = menuHelp->addMenu("Support Links");
    menuHelp->addAction(actionAbout);

    help_discussion = new QAction(fMainWindow);
    help_forum = new QAction(fMainWindow);
    help_bugs = new QAction(fMainWindow);
    help_features = new QAction(fMainWindow);
    help_download = new QAction(fMainWindow);
    supportMenu->addAction(help_bugs);
    supportMenu->addAction(help_download);
    
    menubar->addMenu(menuFile);
    menubar->addMenu(menuPreprocessing);
#ifndef PACKAGE_VIEWER
    menubar->addMenu(menuApp);
#endif
    menubar->addMenu(menuHelp);
    fMainWindow->setMenuBar(menubar);

    menubar->addAction(menuFile->menuAction());
    menubar->addAction(menuPreprocessing->menuAction());
#ifndef PACKAGE_VIEWER
    menubar->addAction(menuApp->menuAction());
#endif
    menubar->addAction(menuHelp->menuAction());

    menuLoadFile->addAction(actionLoad_Nifti_Images);
    menuLoadFile->addAction(actionLoad_Nifti_ROI);

    menuSaveFile->addAction(actionSave_Nifti_Images);
    menuSaveFile->addAction(actionSave_ROI_Images);

    menuFile->addAction(actionExit);

    menuDownload->addAction("GreedyRegistration");

    retranslateUi(fMainWindow);

    m_tabWidget->setCurrentIndex(0);


    QMetaObject::connectSlotsByName(fMainWindow);
  }

  void retranslateUi(QMainWindow *fMainWindow)
  {
 
    fMainWindow->setWindowTitle(QApplication::translate("fMainWindow", "fMainWindow", 0));
    actionLoad_Nifti_Images->setText(QApplication::translate("fMainWindow", "Image(s)", 0));
    actionLoad_Nifti_ROI->setText(QApplication::translate("fMainWindow", "ROI", 0));

    actionSave_Nifti_Images->setText(QApplication::translate("fMainWindow", "Image (NIfTI)", 0));
    actionSave_Dicom_Images->setText(QApplication::translate("fMainWindow", "Image (DICOM)", 0));
    actionSave_ROI_Images->setText(QApplication::translate("fMainWindow", "ROI (NIfTI)", 0));
    actionSave_ROI_Dicom_Images->setText(QApplication::translate("fMainWindow", "ROI (DICOM)", 0));
    actionHelp_Interactions->setText(QApplication::translate("fMainWindow", "Usage", 0));
    help_discussion->setText(QApplication::translate("fMainWindow", "Discussion Forum", 0));
    help_forum->setText(QApplication::translate("fMainWindow", "Help Forum", 0));
    help_bugs->setText(QApplication::translate("fMainWindow", "Bugs and Feature", 0));
    help_features->setText(QApplication::translate("fMainWindow", "Feature Requests", 0));
    help_download->setText(QApplication::translate("fMainWindow", "Latest Downloads", 0));
    actionAbout->setText(QApplication::translate("fMainWindow", "About", 0));
    actionExit->setText(QApplication::translate("fMainWindow", "Exit", 0));
    actionAppGeodesic->setText(QApplication::translate("fMainWindow", "Geodesic segmentation", 0));

    thresholdLabel->setText(QApplication::translate("fMainWindow", "Threshold", 0));
    presetLabel->setText(QApplication::translate("fMainWindow", "Preset", 0));
    windowLabel->setText(QApplication::translate("fMainWindow", "Window", 0));
    levelLabel->setText(QApplication::translate("fMainWindow", "Level", 0));
    m_tabWidget->setTabText(m_tabWidget->indexOf(tumorPanel), QApplication::translate("fMainWindow", "Seed Points", 0));
    m_tabWidget->setTabText(m_tabWidget->indexOf(drawingPanel), QApplication::translate("fMainWindow", "Drawing", 0));
    m_tabWidget->setTabText(m_tabWidget->indexOf(imagesPanel), QApplication::translate("fMainWindow", "Images", 0));
  }

};

namespace Ui {
  class fMainWindow : public Ui_fMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FMAINWINDOW_H
