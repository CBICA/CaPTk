
#ifndef UI_fFeaturePanel_H
#define UI_fFeaturePanel_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTableWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include <qcombobox.h>
#include <QToolButton>
#include <qgroupbox.h>
#include <qsize.h>
#include <map>
#include "CAPTk.h"

QT_BEGIN_NAMESPACE

class Ui_fFeaturePanel
{
public:

  QPushButton *HelpButton;
  QPushButton * m_btnCompute;
  QPushButton * m_btnAdvanced;
  QComboBox * m_cmbFeatureType;
  QLineEdit* m_txtSaveFileName;


  QCheckBox* m_1st_oderstatistics;
  QCheckBox* m_Volumetricfeature;
  QCheckBox* m_Morphologicfeature;
  QCheckBox* m_Histogramfeature;
  QCheckBox* m_Runlenghtfeature;
  QCheckBox* m_Co_occurancefeature;
  QCheckBox* m_LBPfeature;
  QCheckBox* m_NeighbGrayTone;
  QCheckBox* m_GraylevelSizeZone;

  QRadioButton *radio1;
  QRadioButton *radio2;

  QLineEdit* m_roi;
  QLineEdit* m_roi_label;

  QPushButton* m_btnBrowseSaveFile;
  QCheckBox* csv_format;
  QCheckBox* xml_format;

  std::string tempFolderLocation;
  std::map<std::string,QCheckBox*> m_featureCheckBoxMap;
  void setupFeatureCheckBoxMap()
  {
    m_featureCheckBoxMap.clear();// TBD patma double check the names again 
    m_featureCheckBoxMap["Intensity"] = m_1st_oderstatistics;
    m_featureCheckBoxMap["Histogram"] = m_Histogramfeature;
    m_featureCheckBoxMap["GLRLM"] = m_Runlenghtfeature;
    m_featureCheckBoxMap["Morphologic"] = m_Morphologicfeature;//TBD 
    m_featureCheckBoxMap["Volumetric"] = m_Volumetricfeature;
    m_featureCheckBoxMap["GLCM"] = m_Co_occurancefeature;
    m_featureCheckBoxMap["LBP2D"] = m_LBPfeature;
    m_featureCheckBoxMap["NGTDM"] = m_NeighbGrayTone;
    m_featureCheckBoxMap["GLSZM"] = m_GraylevelSizeZone;
  }
  void setupUi(QWidget *parent)
  {
    int buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");
    m_cmbFeatureType = new QComboBox();
    m_cmbFeatureType->setToolTip(QString("Set the type of Feature to be extracted/saved"));
    m_cmbFeatureType->addItem("Custom");
    //m_cmbFeatureType->addItem("Lung_SBRT");
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
    QVBoxLayout* featureLayout2= new QVBoxLayout();
    QVBoxLayout* featureLayout3 = new QVBoxLayout();
    QVBoxLayout* featureLayout4= new QVBoxLayout();
    
    featureLayout1->addWidget(m_cmbFeatureType);
    featureLayout1->addWidget(m_btnAdvanced);
    featureLayout1->addStretch();
  
    m_1st_oderstatistics = new QCheckBox("1st Order Statistics");
    m_1st_oderstatistics->setToolTip(QString("Minimum, Maximum, Mean, Variance, Standard Deviation, Skewness, Kurtosis"));
    m_1st_oderstatistics->setEnabled(true);
    
    m_Morphologicfeature = new QCheckBox("Morphologic");
    m_Morphologicfeature->setToolTip(QString("Elongation, Flatness, Perimeter and First 2 Principal Components."));

    m_Volumetricfeature = new QCheckBox("Volumetric");
    m_Volumetricfeature->setToolTip(QString("Number_of_pixels/voxels, Volume"));
    m_Histogramfeature = new    QCheckBox("Histogram-related");
    m_Histogramfeature->setToolTip(QString("Binning of intensities and peak estimation (click on '?' for details)."));
    m_Runlenghtfeature = new QCheckBox("GLRLM (Run Length)");
    m_Runlenghtfeature->setToolTip(QString("Creates a Gray-Level Run Length Matrix (GLRLM) for an offset of 2 voxels in 3D."));
    m_Co_occurancefeature = new QCheckBox("GLCM (Co-occurence)");
    m_Co_occurancefeature->setToolTip(QString("Creates a Gray-Level Co-occurence Matrix (GLCM) for an offset of 1 voxel in 3D."));
    m_LBPfeature = new QCheckBox("LBP (Local Binary Patterns)");
    m_LBPfeature->setToolTip(QString("Calculates LBPs for 3D image and given mask."));
    m_LBPfeature->setEnabled(false);
    m_NeighbGrayTone = new QCheckBox("NGTDM (Tone Difference)");
    m_NeighbGrayTone->setToolTip(QString("Calculates Neighborhood Gray-Tone Difference Matrix."));
    m_GraylevelSizeZone = new QCheckBox("GLSZM (Size Zone)");
    m_GraylevelSizeZone->setToolTip(QString("Calculates Gray-Level Size Zone Matrix."));

    QGroupBox *customization = new QGroupBox(("Customization"));
    QGroupBox *texturegroup = new QGroupBox(("Textural"));

    QGroupBox *intensitygroup = new QGroupBox(("Intensity-Based"));
    QGroupBox *Volumetric = new QGroupBox(("Miscellaneous"));
    
    QVBoxLayout* featureLayouthbox = new QVBoxLayout();
    featureLayouthbox->addWidget(intensitygroup);
    featureLayouthbox->addWidget(Volumetric);

    featureLayout2->addWidget(m_1st_oderstatistics);
    featureLayout2->addWidget(m_Histogramfeature);

    featureLayout3->addWidget(m_Volumetricfeature);
    featureLayout3->addWidget(m_Morphologicfeature);
    featureLayout3->addStretch();

    featureLayout4->addWidget(m_Co_occurancefeature);
    featureLayout4->addWidget(m_Runlenghtfeature);
    featureLayout4->addWidget(m_GraylevelSizeZone);
    featureLayout4->addWidget(m_NeighbGrayTone);
    featureLayout4->addWidget(m_LBPfeature);
    featureLayout4->addStretch();
    
    customization->setLayout(featureLayout1);
    texturegroup->setLayout(featureLayout4);  
    intensitygroup->setLayout(featureLayout2);  
    Volumetric->setLayout(featureLayout3);
  
    featureLayout->addWidget(customization);
    featureLayout->addLayout(featureLayouthbox);
   // featureLayout->addWidget(intensitygroup);
    featureLayout->addWidget(texturegroup);
  //  featureLayout->addWidget(Volumetric);

    featureGroup->setLayout(featureLayout);
    
    QGroupBox *selectionGroup = new QGroupBox(("Input Selection"));
    QGroupBox *imageselect = new QGroupBox(("Image Selection"));
    QGroupBox *maskselect = new QGroupBox(("Mask Selection"));

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
    m_roi_label = new QLineEdit(std::string("ED").c_str());
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

    csv_format = new QCheckBox(".csv");
    //csv_format->setToolTip(QString("Saves file in csv format"));
    xml_format = new QCheckBox(".xml");
    //xml_format->setToolTip(QString("Saves file in xml format"));
    
    // this is done solely for the reason for saving everything in the tempDir
    // this folder is deleted after this command to ensure no conflict with tempFolderLocation from fMainWindow
    m_txtSaveFileName = new QLineEdit(std::string(loggerFolder + "features.csv").c_str());    
    m_txtSaveFileName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

   // flLayout->addWidget(label);
    flLayout->addWidget(m_txtSaveFileName);
    flLayout->addWidget(m_btnBrowseSaveFile);
    //flLayout->addWidget(csv_format);
    //flLayout->addWidget(xml_format);
    saveGroup_file->setLayout(flLayout);

    QVBoxLayout* saveLayout = new QVBoxLayout();
    //saveLayout->addLayout(flLayout);
    saveLayout->addWidget(saveGroup_file);
    saveLayout->addWidget(m_btnCompute);
    saveLayout->addStretch();
    saveGroup->setLayout(saveLayout);


    HelpButton = new QPushButton();
    std::string iconDir;
#ifdef DEVELOPER_MODE
    iconDir = QApplication::applicationDirPath().toStdString() + "/../../data/icons/";
#else
    iconDir = QApplication::applicationDirPath().toStdString() + "/../data/icons/";
#endif
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
