
#ifndef UI_fFeaturePanel_H
#define UI_fFeaturePanel_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QFrame>
// #include <QtGui/QGridLayout>
// #include <QtGui/QHBoxLayout>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QPushButton>
// #include <QtGui/QRadioButton>
// #include <QtGui/QSpacerItem>
// #include <QtGui/QTableWidget>
// #include <QtGui/QVBoxLayout>
// #include <QtGui/QWidget>
// NEW CHANGES
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <qcombobox.h>
#include <QToolButton>
#include <qgroupbox.h>
#include <qsize.h>
#include <map>
//#include "CAPTk.h"
#include <QCheckBox>
#include "CaPTkGUIUtils.h"
#include <QLineEdit>

QT_BEGIN_NAMESPACE

class Ui_fFeaturePanel
{
public:

  QPushButton *HelpButton;
  QPushButton * m_btnCompute;
  QPushButton * m_btnAdvanced;
  QComboBox * m_cmbFeatureType;
  QLineEdit* m_txtSaveFileName;


  QCheckBox* m_generic;
  QCheckBox* m_firstOrderStatistics;
  QCheckBox* m_Volumetric;
  QCheckBox* m_Morphologic;
  QCheckBox* m_Histogram;
  QCheckBox* m_GLRLM;
  QCheckBox* m_GLCM;
  QCheckBox* m_LBP;
  QCheckBox* m_NGTDM;
  QCheckBox* m_NGLDM;
  QCheckBox* m_GLSZM;
  QCheckBox* m_Lattice;
  QCheckBox* m_Laws;
  QCheckBox* m_PowerSpectrum;
  QCheckBox* m_Gabor;
  QCheckBox* m_Edges;

  QLineEdit* m_patientID_label;

  QRadioButton *radio1;
  QRadioButton *radio2;

  QLineEdit* m_roi;
  QLineEdit* m_roi_label;

  QCheckBox* m_verticalConcat;
  QPushButton* m_btnBrowseSaveFile;
  QCheckBox* csv_format;
  QCheckBox* xml_format;

  //std::vector< std::string > m_featureFiles;

  std::string m_tempFolderLocation;
  std::map<std::string,QCheckBox*> m_featureCheckBoxMap;
  void setupFeatureCheckBoxMap()
  {
    m_featureCheckBoxMap.clear();
    m_featureCheckBoxMap["Generic"] = m_generic;
    m_featureCheckBoxMap["Intensity"] = m_firstOrderStatistics;
    m_featureCheckBoxMap["Histogram"] = m_Histogram;
    m_featureCheckBoxMap["GLRLM"] = m_GLRLM;
    m_featureCheckBoxMap["Morphologic"] = m_Morphologic;//TBD
    m_featureCheckBoxMap["Volumetric"] = m_Volumetric;
    m_featureCheckBoxMap["GLCM"] = m_GLCM;
    m_featureCheckBoxMap["LBP"] = m_LBP;
    m_featureCheckBoxMap["NGTDM"] = m_NGTDM;
    m_featureCheckBoxMap["NGLDM"] = m_NGLDM;
    m_featureCheckBoxMap["GLSZM"] = m_GLSZM;
    m_featureCheckBoxMap["Lattice"] = m_Lattice;
    m_featureCheckBoxMap["Laws"] = m_Laws;
    m_featureCheckBoxMap["PowerSpectrum"] = m_PowerSpectrum;
    m_featureCheckBoxMap["GaborWavelets"] = m_Gabor;
    m_featureCheckBoxMap["EdgeEnhancement"] = m_Edges;
  }
  void setupUi(QWidget *parent)
  {
    int buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");
    m_cmbFeatureType = new QComboBox();
    m_cmbFeatureType->setToolTip(QString("Set the type of Feature to be extracted/saved"));
    m_cmbFeatureType->addItem("Custom");
    m_cmbFeatureType->addItem("Custom_Lattice_2D");
    m_cmbFeatureType->addItem("SBRT_Lung");

    //m_featureFiles = cbica::stringSplit(std::string(FeatureDefinitions), " ");
    //m_featureFiles.erase(m_featureFiles.begin()); // because an empty string always gets appended first

    //for (size_t i = 0; i < m_featureFiles.size(); i++)
    //{
    //  if (m_featureFiles[i].find("1_params") != std::string::npos)
    //  {
    //    m_cmbFeatureType->addItem("Custom");
    //  }
    //  else if (m_featureFiles[i].find("2_params") != std::string::npos)
    //  {
    //    m_cmbFeatureType->addItem("Custom_Lattice");
    //  }
    //  else if (m_featureFiles[i].find("SBRT_Lung") != std::string::npos)
    //  {
    //    m_cmbFeatureType->addItem("SBRT_Lung");
    //  }
    //}
    m_cmbFeatureType->setFixedWidth(buttonWidth + buttonWidth / 5);
    fixComboBox(m_cmbFeatureType);
	  m_btnCompute = new QPushButton(parent);
    m_btnCompute->setText(QString("Compute + Save"));
    m_btnAdvanced = new QPushButton(parent);
    m_btnAdvanced->setText(QString("Advanced ..."));
    m_btnCompute->setToolTip(QString("Compute and Save features"));
    m_btnCompute->setFixedWidth(buttonWidth + buttonWidth / 5);

    QGroupBox* featureGroup = new QGroupBox("Features");
    QHBoxLayout* featureLayout = new QHBoxLayout();
    QVBoxLayout* featureLayout1 = new QVBoxLayout();
    QVBoxLayout* featureLayout2 = new QVBoxLayout();
    QVBoxLayout* featureLayout3 = new QVBoxLayout();
    QVBoxLayout* featureLayout4 = new QVBoxLayout();
    QVBoxLayout* featureLayout5 = new QVBoxLayout();

    featureLayout1->addWidget(m_cmbFeatureType);
    featureLayout1->addWidget(m_btnAdvanced);
    //featureLayout1->addStretch();

    m_generic = new QCheckBox("Setup");
    m_generic->setToolTip(QString("Resampling, Quantization, Interpolation"));
    m_generic->setChecked(true);
    m_generic->setEnabled(false);
    m_firstOrderStatistics = new QCheckBox("1st Order Statistics");
    m_firstOrderStatistics->setToolTip(QString("Minimum, Maximum, Mean, Variance, Standard Deviation, Skewness, Kurtosis"));
    m_firstOrderStatistics->setChecked(true);
    m_firstOrderStatistics->setEnabled(false);
    m_Morphologic = new QCheckBox("Morphologic");
    m_Morphologic->setToolTip(QString("Elongation, Flatness, Perimeter and First 2 Principal Components."));
    m_Volumetric = new QCheckBox("Volumetric");
    m_Volumetric->setToolTip(QString("Number_of_pixels/voxels, Volume"));
    m_Histogram = new    QCheckBox("Histogram-based");
    m_Histogram->setToolTip(QString("Binning of intensities and peak estimation (click on '?' for details)."));
    m_GLRLM = new QCheckBox("GLRLM (Run Length)");
    m_GLRLM->setToolTip(QString("Creates a Gray-Level Run Length Matrix (GLRLM) for an offset of 2 voxels."));
    m_GLCM = new QCheckBox("GLCM (Co-occurrence)");
    m_GLCM->setToolTip(QString("Creates a Gray-Level Co-occurrence Matrix (GLCM) for an offset of 1 voxel."));
    m_LBP = new QCheckBox("LBP (Local Binary Patterns)");
    m_LBP->setToolTip(QString("Calculates LBPs for image and given mask."));
    m_LBP->setCheckState(Qt::CheckState::Unchecked);
    m_LBP->setEnabled(false);
    m_NGTDM = new QCheckBox("NGTDM (Tone Difference)");
    m_NGTDM->setToolTip(QString("Calculates Neighborhood Gray-Tone Difference Matrix."));
    m_NGLDM = new QCheckBox("NGLDM (Level Dependence)");
    m_NGLDM->setToolTip(QString("Calculates Neighborhood Gray-Level Dependence Matrix."));
    m_GLSZM = new QCheckBox("GLSZM (Size Zone)");
    m_GLSZM->setToolTip(QString("Calculates Gray-Level Size Zone Matrix."));
    m_Lattice = new QCheckBox("Lattice Computation");
    m_Lattice->setToolTip(QString("Calculates selected features based on a lattice framework."));
    m_Lattice->setCheckState(Qt::CheckState::Unchecked);
    m_Lattice->setEnabled(false);
    m_Laws = new QCheckBox("Laws");
    m_Laws->setToolTip(QString("Calculates Laws feature(s) for image on ROI."));
    m_Laws->setCheckState(Qt::CheckState::Unchecked);
    m_Laws->setEnabled(false);
    m_PowerSpectrum = new QCheckBox("Power Spectrum");
    m_PowerSpectrum->setToolTip(QString("Calculates Power Spectrum feature(s) for image on ROI."));
    m_PowerSpectrum->setCheckState(Qt::CheckState::Unchecked);
    m_PowerSpectrum->setEnabled(false);
    m_Gabor = new QCheckBox("Gabor Wavelets");
    m_Gabor->setToolTip(QString("Calculates Gabro Wavelet Feature(s) for image on ROI."));
    m_Gabor->setCheckState(Qt::CheckState::Unchecked);
    m_Gabor->setEnabled(false);
    m_Edges = new QCheckBox("Edge Enhancement");
    m_Edges->setToolTip(QString("Calculates Edge Enhancement Feature(s) for image on ROI."));
    m_Edges->setCheckState(Qt::CheckState::Unchecked);
    m_Edges->setEnabled(false);

    QGroupBox *customization = new QGroupBox(("Customization"));
    QGroupBox *texturegroup = new QGroupBox(("Textural"));

    QGroupBox *intensitygroup = new QGroupBox(("Intensity-Based"));
    QGroupBox *Volumetric = new QGroupBox(("Miscellaneous"));
    QGroupBox *extraGroup = new QGroupBox(("Extras"));

    QVBoxLayout* featureLayouthbox = new QVBoxLayout();
    featureLayouthbox->addWidget(customization);
    featureLayouthbox->addWidget(intensitygroup);
    featureLayouthbox->addWidget(Volumetric);

    featureLayout2->addWidget(m_generic);
    featureLayout2->addWidget(m_firstOrderStatistics);
    featureLayout2->addWidget(m_Histogram);

    featureLayout3->addWidget(m_Volumetric);
    featureLayout3->addWidget(m_Morphologic);
    //featureLayout3->addStretch();

    featureLayout4->addWidget(m_GLCM);
    featureLayout4->addWidget(m_GLRLM);
    featureLayout4->addWidget(m_GLSZM);
    featureLayout4->addWidget(m_NGTDM);
    featureLayout4->addWidget(m_NGLDM);
    featureLayout4->addWidget(m_LBP);
    featureLayout4->addStretch();

    featureLayout5->addWidget(m_Laws);
    featureLayout5->addWidget(m_PowerSpectrum);
    featureLayout5->addWidget(m_Gabor);
    featureLayout5->addWidget(m_Edges);
    featureLayout5->addWidget(m_Lattice);
    featureLayout5->addStretch();

    customization->setLayout(featureLayout1);
    texturegroup->setLayout(featureLayout4);
    intensitygroup->setLayout(featureLayout2);
    Volumetric->setLayout(featureLayout3);
    extraGroup->setLayout(featureLayout5);

    //featureLayout->addWidget(customization);
    featureLayout->addLayout(featureLayouthbox);
    featureLayout->addWidget(texturegroup);
    featureLayout->addWidget(extraGroup);

    featureGroup->setLayout(featureLayout);

    QGroupBox *selectionGroup = new QGroupBox(("Input Selection"));
    QGroupBox *subIDselect = new QGroupBox(("Subject ID"));
    QGroupBox *imageselect = new QGroupBox(("Image Selection"));
    QGroupBox *maskselect = new QGroupBox(("Mask Selection"));

    QVBoxLayout* subIDlayout = new QVBoxLayout();
    m_patientID_label = new QLineEdit(std::string("SubjectID").c_str());
    m_patientID_label->setToolTip(QString("Subject ID of the image(s) being processed"));
    m_patientID_label->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    m_patientID_label->setFixedWidth(buttonWidth + buttonWidth / 5);
    subIDlayout->addWidget(m_patientID_label);
    subIDselect->setLayout(subIDlayout);

    QVBoxLayout* featureLayoutin = new QVBoxLayout();
    radio1 = new QRadioButton("Selected Image");
    radio2 = new QRadioButton("All Images");
    featureLayoutin->addWidget(radio1);
    featureLayoutin->addWidget(radio2);
    imageselect->setLayout(featureLayoutin);

  //  QLabel* label1 = new QLabel("Mask Selector:");
    maskselect->setToolTip(QString("Allows user to select a particular mask value for which the features has to be exracted.ex:If mask has 3 labels and user wants label 2"));

    m_roi= new QLineEdit(std::string("1").c_str());
    m_roi->setToolTip(QString("Value(s) of Label to extract features from"));
    m_roi->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    m_roi->setFixedWidth(buttonWidth + 25);
    m_roi_label = new QLineEdit(std::string("LabelText").c_str());
    m_roi_label->setToolTip(QString("Name(s) of corresponding Labels"));
    m_roi_label->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    m_roi_label->setFixedWidth(buttonWidth + buttonWidth / 5);

    QVBoxLayout* featureLayoutmask = new QVBoxLayout();
   // featureLayoutmask->addWidget(label1);
    featureLayoutmask->addWidget(m_roi);
    featureLayoutmask->addWidget(m_roi_label);
    featureLayoutmask->addStretch();
    maskselect->setLayout(featureLayoutmask);

    QVBoxLayout* selectionLayout = new QVBoxLayout();
    selectionLayout->addWidget(subIDselect);
    selectionLayout->addWidget(imageselect);
    selectionLayout->addWidget(maskselect);
    selectionGroup->setMaximumWidth(buttonWidth + buttonWidth / 1.5);
    selectionGroup->setLayout(selectionLayout);


    // Browse output fileName
    QGroupBox* saveGroup = new QGroupBox("Output");
    QGroupBox* saveGroup_file = new QGroupBox("File Selector");
    QHBoxLayout* flLayout =   new QHBoxLayout();
  //  QLabel* label = new QLabel("FileName:");
    m_btnBrowseSaveFile = new QPushButton("...");
    m_btnBrowseSaveFile->setToolTip("Browse to select file");

    m_verticalConcat = new QCheckBox("Vertically Concatenated");
    m_verticalConcat->setToolTip("If not selected, this presents the output in the format for Training Module to use as input");
    //csv_format = new QCheckBox(".csv");
    //csv_format->setToolTip(QString("Saves file in csv format"));
    //xml_format = new QCheckBox(".xml");
    //xml_format->setToolTip(QString("Saves file in xml format"));

    // this is done solely for the reason for saving everything in the tempDir
    // this folder is deleted after this command to ensure no conflict with m_tempFolderLocation from fMainWindow
    m_txtSaveFileName = new QLineEdit(std::string(cbica::getUserHomeDirectory() + "/captk_features.csv").c_str());
    m_txtSaveFileName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

    //flLayout->addWidget(m_verticalConcat);
    flLayout->addWidget(m_txtSaveFileName);
    flLayout->addWidget(m_btnBrowseSaveFile);
    //flLayout->addWidget(csv_format);
    //flLayout->addWidget(xml_format);
    saveGroup_file->setLayout(flLayout);

    QVBoxLayout* saveLayout = new QVBoxLayout();
    //saveLayout->addLayout(flLayout);
    saveLayout->addWidget(m_verticalConcat);
    saveLayout->addWidget(saveGroup_file);
    saveLayout->addWidget(m_btnCompute);
    saveLayout->addStretch();
    saveGroup->setLayout(saveLayout);
    saveGroup->setMaximumWidth(buttonWidth * 3);


    HelpButton = new QPushButton();
    std::string iconDir = getCaPTkDataDir() + "/icons/";
//#ifndef CAPTK_PACKAGE_PROJECT
//    iconDir = captk_currentApplicationPath + "/../../data/icons/";
//#else
//    iconDir = getCaPTkDataDir() + "/icons/";
//#endif
    HelpButton->setIcon(QIcon((iconDir + "help.png").c_str()));
    HelpButton->setToolTip("Get Help");
    //HelpButton->setIconSize(QSize(30, 30));
    QVBoxLayout *helpLayout = new QVBoxLayout();
    helpLayout->addWidget(HelpButton);
    helpLayout->addStretch();

    QHBoxLayout* subLayout = new QHBoxLayout();
    subLayout->addWidget(featureGroup);
    subLayout->addWidget(selectionGroup);
    subLayout->addWidget(saveGroup);
    subLayout->addStretch();
    subLayout->addLayout(helpLayout);

    QVBoxLayout* mainLayout = new QVBoxLayout(parent);
    mainLayout->addLayout(subLayout);
    setupFeatureCheckBoxMap();
    //mainLayout->addStretch();
  }


};

namespace Ui {
  class fFeaturePanel : public Ui_fFeaturePanel {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_fFeaturePanel_H
