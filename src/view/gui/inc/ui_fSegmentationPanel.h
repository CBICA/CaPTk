
#ifndef UI_fSegmentationPanel_H
#define UI_fSegmentationPanel_H

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
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QLineEdit>
#include <qcombobox.h>
#include <QToolButton>
#include <qgroupbox.h>
#include <qsize.h>
#include <map>
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"

QT_BEGIN_NAMESPACE

class Ui_fSegmentationPanel
{
public:

  QPushButton *HelpButton;
  QPushButton * m_btnCompute;
  QLineEdit* m_txtSaveFileName;
  
  QCheckBox* m_seg1;
  QCheckBox* m_seg2;
  QCheckBox* m_seg3;
  QCheckBox* m_seg4;
  QCheckBox* m_seg5;

  std::vector< std::string > m_featureFiles;

  std::string m_tempFolderLocation;
  std::map< std::string,QCheckBox*> m_algoCheckBoxMap;
  void setupFeatureCheckBoxMap()
  {
    m_algoCheckBoxMap.clear();
    m_algoCheckBoxMap["Algorithm_1"] = m_seg1;
    m_algoCheckBoxMap["Algorithm_2"] = m_seg2;
    m_algoCheckBoxMap["Algorithm_3"] = m_seg3;
    m_algoCheckBoxMap["Algorithm_4"] = m_seg4;
    m_algoCheckBoxMap["Algorithm_5"] = m_seg5;
  }
  void setupUi(QWidget *parent)
  {
    int buttonWidth = QLabel().fontMetrics().width("ButtonSize------------------");

	  m_btnCompute = new QPushButton(parent);
    m_btnCompute->setText(QString("Fusion + Save"));
    m_btnCompute->setToolTip(QString("Fusion and Save features"));
    m_btnCompute->setFixedWidth(buttonWidth + buttonWidth / 5);

    QGroupBox* algoGroup = new QGroupBox("Features");
    QHBoxLayout* algoLayout = new QHBoxLayout();
    QVBoxLayout* algoLayout2 = new QVBoxLayout();
    QVBoxLayout* algoLayout3 = new QVBoxLayout();
    QVBoxLayout* algoLayout4 = new QVBoxLayout();
    QVBoxLayout* algoLayout5 = new QVBoxLayout();

    m_seg1 = new QCheckBox("BraTS Algo 1");
    m_seg2 = new QCheckBox("BraTS Algo 2");
    m_seg3 = new QCheckBox("BraTS Algo 3");
    m_seg4 = new QCheckBox("BraTS Algo 4");
    m_seg5 = new QCheckBox("BraTS Algo 5");

    m_seg1->setChecked(true);
    m_seg2->setChecked(true);
    m_seg3->setChecked(true);
    m_seg4->setChecked(false);
    m_seg5->setChecked(true);

    QGroupBox *firstGroup = new QGroupBox(("Algorithms"));

    QVBoxLayout* algoLayouthbox = new QVBoxLayout();
    algoLayouthbox->addWidget(firstGroup);

    algoLayout2->addWidget(m_seg1);
    algoLayout2->addWidget(m_seg2);
    algoLayout2->addWidget(m_seg3);
    algoLayout2->addWidget(m_seg4);
    algoLayout2->addWidget(m_seg5);

    //customization->setLayout(algoLayout1);
    //texturegroup->setLayout(algoLayout4);
    firstGroup->setLayout(algoLayout2);
    //Volumetric->setLayout(algoLayout3);
    //extraGroup->setLayout(algoLayout5);

    algoLayout->addLayout(algoLayouthbox);

    algoGroup->setLayout(algoLayout);

    QGroupBox *selectionGroup = new QGroupBox(("Input Selection"));
    QGroupBox *imageselect = new QGroupBox(("Image Selection"));
    QGroupBox *maskselect = new QGroupBox(("Mask Selection"));


    // Browse output fileName
    QGroupBox* saveGroup = new QGroupBox("Output");
    QGroupBox* saveGroup_file = new QGroupBox("File Selector");
    QHBoxLayout* flLayout =   new QHBoxLayout();
  //  QLabel* label = new QLabel("FileName:");

    // this is done solely for the reason for saving everything in the tempDir
    // this folder is deleted after this command to ensure no conflict with m_tempFolderLocation from fMainWindow
    m_txtSaveFileName = new QLineEdit(std::string(loggerFolder + "outputModel.xml").c_str());
    m_txtSaveFileName->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);

   // flLayout->addWidget(label);
    flLayout->addWidget(m_txtSaveFileName);
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
    subLayout->addWidget(algoGroup);
    //subLayout->addWidget(selectionGroup);
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
  class fSegmentationPanel : public Ui_fSegmentationPanel {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_fSegmentationPanel_H
