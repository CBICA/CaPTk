/********************************************************************************
** Form generated from reading UI file 'fImagesPanel.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FIMAGESPANEL_H
#define UI_FIMAGESPANEL_H

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
// #include <QtGui/QGroupBox>
// #include <QtGui/QCheckBox>
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
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QCheckBox>
#include "CaPTkGUIUtils.h"

QT_BEGIN_NAMESPACE
enum TAB_IMAGES_COLUMNS
{
  TAB_IMAGES_COLUMN_CLOSE, TAB_IMAGES_COLUMN_TYPE, TAB_IMAGES_COLUMN_NAME
};
class Ui_fImagesPanel
{
public:

  QTableWidget * m_imagesTable;
  QTableWidget * 	m_nonVisImagesTable;

  QSlider * m_overlaySlider;
  QCheckBox * m_overlayChkBox;
  QPushButton * m_3dViz;
  QPushButton * m_clearImagesBtn;
  QPushButton * m_clearImagesBtn_temp; // TBD
  QPushButton * HelpButton;
  QPushButton * m_CompareButton;

  void setupUi(QWidget *fImagesPanel)
  {

    const int pixelPad = 10;
    const int nameColID = 2;
    m_imagesTable = new QTableWidget();
    m_imagesTable->setColumnCount(5);
    m_imagesTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_imagesTable->setSelectionMode(QAbstractItemView::SingleSelection);
    m_imagesTable->setHorizontalHeaderItem(0, new QTableWidgetItem("Close"));
    m_imagesTable->setHorizontalHeaderItem(1, new QTableWidgetItem("Type"));
    m_imagesTable->setHorizontalHeaderItem(2, new QTableWidgetItem("Name"));
    m_imagesTable->setHorizontalHeaderItem(3, new QTableWidgetItem("Modality"));
    m_imagesTable->setHorizontalHeaderItem(4, new QTableWidgetItem("Overlay"));

    for (int i = 0; i < m_imagesTable->model()->columnCount(); i++)
    {
      if (i == nameColID)
      {
        // m_imagesTable->horizontalHeader()->setResizeMode(i, QHeaderView::Stretch);
        // NEW CHANGES
        m_imagesTable->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Stretch);
      }
      else
      {
        QString txt = m_imagesTable->model()->headerData(i, Qt::Horizontal).toString();
        int colWidth = m_imagesTable->horizontalHeader()->fontMetrics().width(txt) + pixelPad;
        // m_imagesTable->horizontalHeader()->setResizeMode(i, QHeaderView::Fixed);
        m_imagesTable->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Fixed);
        m_imagesTable->setColumnWidth(i, colWidth);
      }
    }

    m_nonVisImagesTable = new QTableWidget();
    m_nonVisImagesTable->setColumnCount(3);
    m_nonVisImagesTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_nonVisImagesTable->setSelectionMode(QAbstractItemView::SingleSelection);
    m_nonVisImagesTable->setHorizontalHeaderItem(0, new QTableWidgetItem("Close"));
    m_nonVisImagesTable->setHorizontalHeaderItem(1, new QTableWidgetItem("Type"));
    m_nonVisImagesTable->setHorizontalHeaderItem(2, new QTableWidgetItem("Name"));
    for (int i = 0; i < m_nonVisImagesTable->model()->columnCount(); i++)
    {
      if (i == nameColID)
      {
        m_nonVisImagesTable->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Stretch);
      }
      else
      {
        QString txt = m_nonVisImagesTable->model()->headerData(i, Qt::Horizontal).toString();
        int colWidth = m_nonVisImagesTable->horizontalHeader()->fontMetrics().width(txt) + pixelPad;
        m_nonVisImagesTable->horizontalHeader()->setSectionResizeMode(i, QHeaderView::Fixed);
        m_nonVisImagesTable->setColumnWidth(i, colWidth);
      }
    }


    //-------------------------------Images controls--------------------
    m_overlaySlider = new QSlider();
    m_overlaySlider->setEnabled(false);
    m_overlaySlider->setMaximum(10);
    m_overlaySlider->setSliderPosition(m_overlaySlider->maximum() / 2);
    m_overlaySlider->setOrientation(Qt::Horizontal);
    m_overlayChkBox = new QCheckBox("Change Opacity");
    m_overlayChkBox->setChecked(false);
    m_clearImagesBtn = new QPushButton("Close All");
    m_clearImagesBtn_temp = new QPushButton("Close All");
    m_3dViz = new QPushButton("3D Visualizer");
    m_CompareButton = new QPushButton("Comparison Mode");
    m_CompareButton->setCheckable(true);
#ifndef WIN32
    m_3dViz->setDisabled(true);
#endif
    HelpButton = new QPushButton();
    std::string iconDir;
    
    iconDir = getCaPTkDataDir() + "/icons/";
    if (!QDir(cbica::normPath(iconDir).c_str()).exists()) // packaged binary
    {
      if (QDir((std::string(PROJECT_SOURCE_DIR) + "/data/icons/").c_str()).exists()) // developer_mode
      {
        iconDir = std::string(PROJECT_SOURCE_DIR) + "/data/icons/";
      }
    }

    HelpButton->setIcon(QIcon((iconDir + "help.png").c_str()));
    HelpButton->setToolTip("Get Help");
    //m_HelpButton->setIconSize(QSize(30,30));

    //--------------------- Layout-------------------------------------------
    QVBoxLayout* visLayout = new QVBoxLayout();
    QHBoxLayout* visLayoutSub = new QHBoxLayout();
    visLayoutSub->addWidget(new QLabel("Images"));
    visLayoutSub->addStretch();
    visLayoutSub->addWidget(m_overlayChkBox);
    visLayoutSub->addWidget(m_overlaySlider);
    visLayoutSub->addWidget(new QLabel("|"));
    visLayoutSub->addWidget(m_CompareButton);
    visLayoutSub->addWidget(m_3dViz);
    visLayoutSub->addWidget(new QLabel("|"));
    visLayoutSub->addWidget(m_clearImagesBtn);

    visLayout->addLayout(visLayoutSub);
    visLayout->addWidget(m_imagesTable);

    // The below code is a no op now.
    // We are not showing non visualizing images anymore
    QVBoxLayout* nonVisLayout = new QVBoxLayout();
    QHBoxLayout* nonVisLayoutSub = new QHBoxLayout();
    nonVisLayoutSub->addWidget(new QLabel("Non Visualisation Images"));
    nonVisLayoutSub->addStretch();
    nonVisLayoutSub->addWidget(new QLabel("|"));
    nonVisLayoutSub->addWidget(m_clearImagesBtn_temp);
    nonVisLayoutSub->addWidget(HelpButton);

    nonVisLayout->addLayout(nonVisLayoutSub);
    nonVisLayout->addWidget(m_nonVisImagesTable);
    // The above code is a no op now.
    // We are not showing non visualizing images anymore

    QHBoxLayout	* mainLayout = new QHBoxLayout(fImagesPanel);
    mainLayout->addLayout(visLayout);
    //mainLayout->addLayout(nonVisLayout); //We are not showing non visualizing images anymore
  }

};

namespace Ui {
  class fImagesPanel : public Ui_fImagesPanel {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_fImagesPanel_H
