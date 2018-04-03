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
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QSlider>
#include <QtGui/QStatusBar>
#include <QtGui/QTabWidget>
#include <QtGui/QTableWidget>
#include <QtGui/QToolBar>
#include <QtGui/QWidget>
#include <QtGui/QDockWidget>
#include <QSpacerItem>
#include <QGroupBox>
#include <QRadioButton>
#include <QtGui/QSpacerItem>
#include "fTumorPanel.h"
#include "fImagesPanel.h"
#include "fDrawingPanel.h"
#include "fFeaturePanel.h"
#include "fRecurrenceDialog.h"
#include "fRegistrationDialog.h"
#include "fPreprocessingDialog.h"
#include "fSurvivalDialog.h"
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
//#include "fDeepMedicDialog.h"
#include "fFetalBrain.h"

#include "QVTKWidget.h"
#include "fBottomImageInfoTip.h"

QT_BEGIN_NAMESPACE

class Ui_fMainWindow
{
private:

  /**
  \struct ActionAndName

  \brief This is a helper struct to tie an action with its name as a std::string
  */
  struct ActionAndName
  {
    QAction* action;
    std::string name;
  };

  //! Wrap to ensure previous functionality doesn't break
  inline std::vector< ActionAndName > populateStringListInMenu(const std::string &inputList, QMainWindow* inputFMainWindow, QMenu* menuToPopulate, std::string menuAppSubGroup, bool ExcludeGeodesic)
  {
    if (!inputList.empty())
    {
      std::string inputList_wrap = inputList;
      if (inputList_wrap[0] == ' ')
      {
        inputList_wrap.erase(0, 1);
      }
      std::vector< std::string > vectorOfInputs = cbica::stringSplit(inputList_wrap, " ");
      return populateStringListInMenu(vectorOfInputs, inputFMainWindow, menuToPopulate, menuAppSubGroup, ExcludeGeodesic);
    }
    else
    {
      return std::vector< ActionAndName >{};
    }
  }
  /**
  \brief Takes a list of application variables from CMake defines and put it in specified window and menu

  \param inputList The list obtained from CMake variable which is added to cache
  \param inputFMainWindow The current fMainWindow from which the QActions need to inherit
  \param menuToPopulate The QMenu in which the QActions need to be populated *visualizationInputImagesLabel
  \return A vector of ActionAndName structs which ties a QAction to the corresponding name from inputList
  */
  inline std::vector< ActionAndName > populateStringListInMenu(const std::vector< std::string > &vectorOfInputs, QMainWindow* inputFMainWindow, QMenu* menuToPopulate, std::string menuAppSubGroup, bool ExcludeGeodesic)
  {
    std::vector< ActionAndName > returnVector;
    //if (!inputList.empty())
    //{
    //  std::string inputList_wrap = inputList;
    //  if (inputList_wrap[0] == ' ')
    //  {
    //    inputList_wrap.erase(0, 1);
    //  }
    //  std::vector< std::string > vectorOfInputs = cbica::stringSplit(inputList_wrap, " ");
    if (ExcludeGeodesic)
    {
      returnVector.resize(vectorOfInputs.size());
    }
    else
    {
      returnVector.resize(vectorOfInputs.size() + 1);
    }
    size_t returnVecCounter = 0;

    if (!menuAppSubGroup.empty())
    {
      returnVector[returnVecCounter].action = new QAction(inputFMainWindow);
      returnVector[returnVecCounter].action->setObjectName(QString::fromUtf8(std::string("action" + menuAppSubGroup).c_str()));
      returnVector[returnVecCounter].action->setIconText(QString(menuAppSubGroup.c_str()));
      returnVector[returnVecCounter].action->setText(QString(menuAppSubGroup.c_str()));
      returnVector[returnVecCounter].action->setEnabled(false);
      returnVector[returnVecCounter].name = menuAppSubGroup;
      returnVector[returnVecCounter].action->setEnabled(false);
      menuToPopulate->addAction(returnVector[returnVecCounter].action);
      returnVecCounter++;
    }

    for (size_t i = 0; i < vectorOfInputs.size(); i++)
    {
      if ((vectorOfInputs[i] != "FeatureExtraction") && (vectorOfInputs[i] != "Register"))
      {
        returnVector[returnVecCounter].action = new QAction(inputFMainWindow);
        returnVector[returnVecCounter].action->setObjectName(QString::fromUtf8(std::string("action" + vectorOfInputs[i]).c_str()));
        returnVector[returnVecCounter].action->setIconText(QString(vectorOfInputs[i].c_str()));
        returnVector[returnVecCounter].name = vectorOfInputs[i];
#ifdef CAPTK_BUILD_CONSOLE_ONLY
        returnVector[returnVecCounter].action->setEnabled(false);
#endif
        menuToPopulate->addAction(returnVector[returnVecCounter].action);

        //if (vectorOfInputs[i] == "GeodesicSegmentation")
        //{
        //  menuToPopulate->addSeparator();
        //}
        returnVecCounter++;
      }
    }
    //}

    return returnVector;
  }

public:
  fRecurrenceDialog			recurrencePanel;
  fPopulationAtlasDialog	atlasPanel;
  fRegistrationDialog		registrationPanel;
  fPreprocessingDialog	preprocessingPanel;
  fSurvivalPredictor survivalPanel;
  fMolecularSubtypePredictor msubtypePanel;
  fImagingSubtypePredictor isubtypePanel;
  fFetalBrain fetalbrainpanel;

  fSkullStripper skullStrippingPanel;
  fPCAEstimator pcaPanel;
  fPerfusionEstimator perfmeasuresPanel;
  fDiffusionEstimator diffmeasuresPanel;
  fDCM2NIfTIConverter dcmConverter;
  //fDeepMedicDialog deepMedicDialog;
  fHistoMatcher histoMatchPanel;
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
  //QMenu* menuAbout;
  //QMenu* menuShortcuts;
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

  //QAction *actionDownload_Full;
  //QAction *actionDownload_WhiteStripe;
  //QAction *actionDownload_Survival;
  //QAction *actionDownload_SBRT;
  //QAction *actionDownload_LIBRA;
  //QAction *actionDownload_Recurrence;
  //QAction *actionDownload_PopulationAtlas;
  //QAction *actionDownload_MolecularSubTypes;
  //QAction *actionDownload_Geodesic;
  //QAction *actionDownload_EGFRvIII_PHI;
  //QAction *actionDownload_Confetti;

  QAction *actionAppEGFR;
  QAction *actionAppRecurrence;
  QAction *actionAppGeodesic;

  // initialize vectors of Actions and Names so that the process can be automated and the QAction is tied to its corresponding Name
  std::vector< ActionAndName >
    vectorOfGBMApps, // GBM-specific applications
    vectorOfBreastApps, // breast-specific applications
    vectorOfLungApps, // lung-specific applications
    vectorOfMiscApps, // the rest
    vectorOfPreprocessingActionsAndNames; // for preprocessing algorithms

  //ActionAndName app_directionalityEstimation;

  // obtain list from CMake variables using populateStringListInMenu() function
  std::vector< std::string >
    m_nativeApps, // native CPP applications
    m_preprocessApps, // native pre-processing routines
    m_pyCLIApps, // python command line applications
    m_pyGUIApps; // python graphical applications

  //std::vector< NonNativeApp > m_allNonNativeApps;
  std::map< std::string, std::string > m_allNonNativeApps;


  //std::vector< ActionAndName > vectorOfDownloadLinks;


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
    //overallGridLayout->addWidget(InfoPanelWidget, 4, 0, 1, 2);
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

    //overallGridLayout->addWidget(preferencesGroupBox, 4, 2, 1, 1);


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

    //actionDownload_Full = new QAction(fMainWindow);
    //actionDownload_Full->setText("Full");
    //actionDownload_Full->setToolTip("Download Full Data");
    //actionDownload_WhiteStripe = new QAction(fMainWindow);
    //actionDownload_WhiteStripe->setText("WhiteStripe");
    //actionDownload_WhiteStripe->setToolTip("Download WhiteStripe");
    //actionDownload_Survival = new QAction(fMainWindow);
    //actionDownload_Survival->setText("Survival");
    //actionDownload_Survival->setToolTip("Download Survival");
    //actionDownload_SBRT = new QAction(fMainWindow);
    //actionDownload_SBRT->setText("SBRT");
    //actionDownload_SBRT->setToolTip("Download SBRT");
    //actionDownload_LIBRA = new QAction(fMainWindow);
    //actionDownload_LIBRA->setText("LIBRA");
    //actionDownload_LIBRA->setToolTip("Download LIBRA");
    //actionDownload_Recurrence = new QAction(fMainWindow);
    //actionDownload_Recurrence->setText("Recurrence");
    //actionDownload_Recurrence->setToolTip("Download Recurrence");
    //actionDownload_PopulationAtlas = new QAction(fMainWindow);
    //actionDownload_PopulationAtlas->setText("PopulationAtlas");
    //actionDownload_PopulationAtlas->setToolTip("Download PopulationAtlas");
    //actionDownload_MolecularSubTypes = new QAction(fMainWindow);
    //actionDownload_MolecularSubTypes->setText("MolecularSubTypes");
    //actionDownload_MolecularSubTypes->setToolTip("Download MolecularSubTypes");
    //actionDownload_Geodesic = new QAction(fMainWindow);
    //actionDownload_Geodesic->setText("Geodesic");
    //actionDownload_Geodesic->setToolTip("Download Geodesic");
    //actionDownload_EGFRvIII_PHI = new QAction(fMainWindow);
    //actionDownload_EGFRvIII_PHI->setText("EGFRvIII_PHI");
    //actionDownload_EGFRvIII_PHI->setToolTip("Download EGFRvIII_PHI");
    //actionDownload_Confetti = new QAction(fMainWindow);
    //actionDownload_Confetti->setText("Confetti");
    //actionDownload_Confetti->setToolTip("Download Confetti");

    actionAbout = new QAction(fMainWindow);

    menuHelp = new QMenu("Help");
    //menuAbout = new QMenu("About");
    //menuShortcuts = new QMenu("Shortcuts");
    menuHelp->addAction(actionHelp_Interactions);
    menuDownload = menuHelp->addMenu("Sample Data");
    auto supportMenu = menuHelp->addMenu("Support Links");
    menuHelp->addAction(actionAbout);

    help_discussion = new QAction(fMainWindow);
    help_forum = new QAction(fMainWindow);
    help_bugs = new QAction(fMainWindow);
    help_features = new QAction(fMainWindow);
    help_download = new QAction(fMainWindow);
    supportMenu->addAction(help_discussion);
    supportMenu->addAction(help_forum);
    supportMenu->addAction(help_bugs);
    supportMenu->addAction(help_features);
    supportMenu->addAction(help_download);

    //downloadMenu->addAction(actionDownload_Full);
    //downloadMenu->addAction(actionDownload_EGFRvIII_PHI);
    //downloadMenu->addAction(actionDownload_Recurrence);
    //downloadMenu->addAction(actionDownload_Survival);
    //downloadMenu->addAction(actionDownload_LIBRA);
    //downloadMenu->addAction(actionDownload_WhiteStripe);
    //downloadMenu->addAction(actionDownload_PopulationAtlas);
    //downloadMenu->addAction(actionDownload_MolecularSubTypes);
    //downloadMenu->addAction(actionDownload_Geodesic);
    //downloadMenu->addAction(actionDownload_SBRT);
    //downloadMenu->addAction(actionDownload_Confetti);

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
    //menuSaveFile->addAction(actionSave_Dicom_Images);
    menuSaveFile->addAction(actionSave_ROI_Images);
    //menuSaveFile->addAction(actionSave_ROI_Dicom_Images);

    menuFile->addAction(actionExit);
    //menuAbout->addAction(actionAbout);
    //menuShortcuts->addAction(actionHelp);

    std::string nonNativeAppPaths_wrap = std::string(CAPTK_APP_LIST_PY_GUI);
    if (nonNativeAppPaths_wrap[0] == ' ')
    {
      nonNativeAppPaths_wrap.erase(0, 1);
    }
    nonNativeAppPaths_wrap = nonNativeAppPaths_wrap + " itksnap";
    m_pyGUIApps = cbica::stringSplit(nonNativeAppPaths_wrap, " ");
    nonNativeAppPaths_wrap = std::string(CAPTK_APP_LIST_PY_CLI);
    if (nonNativeAppPaths_wrap[0] == ' ')
    {
      nonNativeAppPaths_wrap.erase(0, 1);
    }

    m_pyCLIApps = cbica::stringSplit(nonNativeAppPaths_wrap, " ");
    //m_allNonNativeApps.resize(m_pyGUIApps.size() + m_pyCLIApps.size());
    size_t allAppCounter = 0;
    for (size_t i = 0; i < m_pyGUIApps.size(); i++)
    {
      if (m_pyGUIApps[i] == "confetti")
      {
        m_pyGUIApps[i] = "ConfettiGUI";
      }
      if ((m_pyGUIApps[i] == "librabatch") || (m_pyGUIApps[i] == "librasingle"))
      {
        m_pyGUIApps[i] = "libra";
      }
      if (m_pyGUIApps[i] == "SBRT_Segment")
      {
        m_pyGUIApps[i] = "SBRT_Lung_Segment";
      }
      if (m_pyGUIApps[i] == "SBRT_Analyze")
      {
        m_pyGUIApps[i] = "SBRT_Lung_Analyze";
      }
      //m_allNonNativeApps[allAppCounter].name = m_pyGUIApps[i];
#ifndef CAPTK_PACKAGE_PROJECT
      if (m_pyGUIApps[i] == "itksnap")
      {
        m_allNonNativeApps[m_pyGUIApps[i]] = std::string(PROJECT_SOURCE_DIR) + "src/applications/individualApps/itksnap/bin/";
        //m_allNonNativeApps[allAppCounter].path = std::string(PROJECT_SOURCE_DIR) + "src/applications/individualApps/" + m_allNonNativeApps[allAppCounter].name + "/bin/";
      }
      else
      {
        m_allNonNativeApps[m_pyGUIApps[i]] = std::string(PROJECT_SOURCE_DIR) + "src/applications/individualApps/";
        if ((m_pyGUIApps[i] == "SBRT_Segment") || (m_pyGUIApps[i] == "SBRT_Analyze"))
        {
          m_allNonNativeApps[m_pyGUIApps[i]] = QApplication::applicationDirPath().toStdString() + "/";
        }
        //m_allNonNativeApps[allAppCounter].path = std::string(PROJECT_SOURCE_DIR) + "src/applications/individualApps/";
      }
#else
      //m_allNonNativeApps[allAppCounter].path = QApplication::applicationDirPath().toStdString() + "/";
      m_allNonNativeApps[m_pyGUIApps[i]] = QApplication::applicationDirPath().toStdString() + "/";
#endif

      m_allNonNativeApps[m_pyGUIApps[i]] += m_pyGUIApps[i]
#ifdef _WIN32
        + ".exe"
#endif
        ;
      allAppCounter++;
    }

    for (size_t i = 0; i < m_pyCLIApps.size(); i++)
    {
      //m_allNonNativeApps[allAppCounter].name = m_pyCLIApps[i];
      m_allNonNativeApps[m_pyCLIApps[i]] =
#ifndef CAPTK_PACKAGE_PROJECT
        std::string(PROJECT_SOURCE_DIR) + "src/applications/individualApps/"
#else
        QApplication::applicationDirPath().toStdString() + "/"
#endif
        ;

      m_allNonNativeApps[m_pyCLIApps[i]] += m_pyCLIApps[i]
#ifdef _WIN32
        + ".exe"
#endif
        ;
    }

    // TBD: this needs to be controlled from CMake and not hard-coded here
    auto brainAppList = " WhiteStripe PopulationAtlases confetti EGFRvIIISurrogateIndex RecurrenceEstimator SurvivalPredictor MolecularSubtypePredictor";
    auto breastAppList = " librasingle librabatch";
    auto lungAppList = /*" SBRT_Segment SBRT_Analyze"*/" SBRT_Segment";
    auto miscAppList = " itksnap GeodesicSegmentation DirectionalityEstimate PerfusionDerivatives PerfusionPCA DiffusionDerivatives";

    vectorOfGBMApps = populateStringListInMenu(brainAppList, fMainWindow, menuApp, "Brain Cancer", false);
    //app_directionalityEstimation.action = new QAction(fMainWindow);
    //app_directionalityEstimation.action->setObjectName(QString::fromUtf8(std::string("actionDirectionalityEstimate").c_str()));
    //app_directionalityEstimation.action->setIconText(QString("Directionality Estimation"));
    //app_directionalityEstimation.action->setText("Directionality Estimation");
    //app_directionalityEstimation.name = "DirectionalityEstimate";
    //menuApp->addAction(app_directionalityEstimation.action);
    menuApp->addSeparator();
    vectorOfBreastApps = populateStringListInMenu(breastAppList, fMainWindow, menuApp, "Breast Cancer", false);
    menuApp->addSeparator();
    vectorOfLungApps = populateStringListInMenu(lungAppList, fMainWindow, menuApp, "Lung Cancer", false);
    menuApp->addSeparator();
    vectorOfMiscApps = populateStringListInMenu(miscAppList, fMainWindow, menuApp, "Miscellaneous", false);
    vectorOfPreprocessingActionsAndNames = populateStringListInMenu(std::string(PREPROCESS_ALGOS), fMainWindow, menuPreprocessing, "", false);

    for (const auto &currentActionAndName : vectorOfGBMApps)
    {
      if (currentActionAndName.name != "Brain Cancer")
      {
        if (currentActionAndName.name == "confetti")
        {
          menuDownload->addAction("Confetti");
        }
        else
        {
          menuDownload->addAction(currentActionAndName.name.c_str());
        }
      }
    }

    for (const auto &currentActionAndName : vectorOfMiscApps)
    {
      if (currentActionAndName.name != "Miscellaneous")
      {
        if (currentActionAndName.name != "itksnap")
        {
          if (!currentActionAndName.name.empty())
          {
            menuDownload->addAction(currentActionAndName.name.c_str());
          }
        }
      }
    }

    bool libraCheck = false;
    for (const auto &currentActionAndName : vectorOfBreastApps)
    {
      if (currentActionAndName.name != "Breast Cancer")
      {
        if (!libraCheck)
        {
          if (currentActionAndName.name.find("libra") != std::string::npos)
          {
            libraCheck = true;
            menuDownload->addAction("LIBRA");
          }
        }
      }
    }

    bool sbrtCheck = false;
    for (const auto &currentActionAndName : vectorOfLungApps)
    {
      if (currentActionAndName.name != "Lung Cancer")
      {
        if (!sbrtCheck)
        {
          if ((currentActionAndName.name.find("SBRT") != std::string::npos))
          {
            sbrtCheck = true;
            menuDownload->addAction("SBRT");
          }
        }
      }
    }

    retranslateUi(fMainWindow);

    m_tabWidget->setCurrentIndex(0);


    QMetaObject::connectSlotsByName(fMainWindow);
  }

  void retranslateUi(QMainWindow *fMainWindow)
  {
    fMainWindow->setWindowTitle(QApplication::translate("fMainWindow", "fMainWindow", 0, QApplication::UnicodeUTF8));
    actionLoad_Nifti_Images->setText(QApplication::translate("fMainWindow", "Image(s)", 0, QApplication::UnicodeUTF8));
    actionLoad_Nifti_ROI->setText(QApplication::translate("fMainWindow", "ROI", 0, QApplication::UnicodeUTF8));

    actionSave_Nifti_Images->setText(QApplication::translate("fMainWindow", "Image (NIfTI)", 0, QApplication::UnicodeUTF8));
    actionSave_Dicom_Images->setText(QApplication::translate("fMainWindow", "Image (DICOM)", 0, QApplication::UnicodeUTF8));
    actionSave_ROI_Images->setText(QApplication::translate("fMainWindow", "ROI (NIfTI)", 0, QApplication::UnicodeUTF8));
    actionSave_ROI_Dicom_Images->setText(QApplication::translate("fMainWindow", "ROI (DICOM)", 0, QApplication::UnicodeUTF8));
    actionHelp_Interactions->setText(QApplication::translate("fMainWindow", "Usage", 0, QApplication::UnicodeUTF8));
    help_discussion->setText(QApplication::translate("fMainWindow", "Discussion Forum", 0, QApplication::UnicodeUTF8));
    help_forum->setText(QApplication::translate("fMainWindow", "Help Forum", 0, QApplication::UnicodeUTF8));
    help_bugs->setText(QApplication::translate("fMainWindow", "Bug Tracker", 0, QApplication::UnicodeUTF8));
    help_features->setText(QApplication::translate("fMainWindow", "Feature Requests", 0, QApplication::UnicodeUTF8));
    help_download->setText(QApplication::translate("fMainWindow", "Latest Downloads", 0, QApplication::UnicodeUTF8));
    actionAbout->setText(QApplication::translate("fMainWindow", "About", 0, QApplication::UnicodeUTF8));
    actionExit->setText(QApplication::translate("fMainWindow", "Exit", 0, QApplication::UnicodeUTF8));
    actionAppGeodesic->setText(QApplication::translate("fMainWindow", "Geodesic segmentation", 0, QApplication::UnicodeUTF8));

    thresholdLabel->setText(QApplication::translate("fMainWindow", "Threshold", 0, QApplication::UnicodeUTF8));
    presetLabel->setText(QApplication::translate("fMainWindow", "Preset", 0, QApplication::UnicodeUTF8));
    windowLabel->setText(QApplication::translate("fMainWindow", "Window", 0, QApplication::UnicodeUTF8));
    levelLabel->setText(QApplication::translate("fMainWindow", "Level", 0, QApplication::UnicodeUTF8));
    m_tabWidget->setTabText(m_tabWidget->indexOf(tumorPanel), QApplication::translate("fMainWindow", "Seed Points", 0, QApplication::UnicodeUTF8));
    m_tabWidget->setTabText(m_tabWidget->indexOf(drawingPanel), QApplication::translate("fMainWindow", "Drawing", 0, QApplication::UnicodeUTF8));
    m_tabWidget->setTabText(m_tabWidget->indexOf(imagesPanel), QApplication::translate("fMainWindow", "Images", 0, QApplication::UnicodeUTF8));
  }

};

namespace Ui {
  class fMainWindow : public Ui_fMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FMAINWINDOW_H
