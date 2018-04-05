///////////////////////////////////////////////////////////////////////////////////////
// fImagesPanel.cxx
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
// License Agreement: http://www.med.upenn.edu/sbia/software/license.html
///////////////////////////////////////////////////////////////////////////////////////

#include "fImagesPanel.h"
#include "CAPTk.h"


fImagesPanel::fImagesPanel(QWidget * parent) : QWidget(parent)
{
  setupUi(this);

  connect(m_imagesTable, SIGNAL(itemSelectionChanged()), this, SLOT(ImageTableSelectionChanged()));
  connect(m_imagesTable, SIGNAL(itemClicked(QTableWidgetItem*)), this, SLOT(ImageTableSelectionChanged(QTableWidgetItem*)));
  m_overlayChkBox->setEnabled(true);

  connect(m_overlayChkBox, SIGNAL(stateChanged(int)), this, SLOT(overlayUseStateChanged(int)));
  connect(m_overlaySlider, SIGNAL(valueChanged(int)), this, SLOT(overlaySliderChanged(int)));
}

void fImagesPanel::ImageTableSelectionChanged()
{
  emit sigImageTableSelectionChanged();
  if (m_overlayChkBox->isChecked() && getSelectedOverlay() != NULL)
  {
    emit overlayChanged();
  }
}
void fImagesPanel::overlayUseStateChanged(int value)
{
  m_overlaySlider->setEnabled(value);
  emit sigOverlayCheckBoxChanged(value);
}
void fImagesPanel::overlayChanged()
{
  emit sigOverlayChanged();
}
void fImagesPanel::overlaySliderChanged(int value)
{
  emit sigOverlaySliderChanged(value);
}
void fImagesPanel::imageModalityChanged(int value)
{
  emit sigImageModalityChanged(value);
}

void fImagesPanel::ImageTableSelectionChanged(QTableWidgetItem*)
{

}

void fImagesPanel::NewImageLoaded(QString idstr, const std::string &filename, int rowIndex, const std::string &imageSubTypeStr, const int imgSubtype, const QObject* caller)
{
  m_imagesTable->setRowCount(rowIndex + 1);
  QFileInfo fileinfo(filename.c_str());

  QTableWidgetItem *item1 = new QTableWidgetItem(fileinfo.baseName());
  item1->setData(Qt::UserRole, idstr.toStdString().c_str());
  item1->setFlags(item1->flags() &  ~Qt::ItemIsEditable);

  QTablePushButton* cButton1 = new QTablePushButton;
  cButton1->setItem(item1);
  cButton1->setText(tr("X"));
  connect(cButton1, SIGNAL(clickedInto(QTableWidgetItem*)), caller, SLOT(CloseImage(QTableWidgetItem*)));


  QLabel * label1 = new QLabel;
  label1->setText(QString::fromStdString(imageSubTypeStr));
  m_imagesTable->setCellWidget(rowIndex, IMAGES_COLUMN_CLOSE, cButton1);
  m_imagesTable->setCellWidget(rowIndex, IMAGES_COLUMN_TYPE, label1);
  m_imagesTable->setItem(rowIndex, IMAGES_COLUMN_NAME, item1);

  // TBD: this needs to pick up from ImageModalityType & ImageModalityString from CAPTk.h and not be hard-coded 
  auto modalitySwitcher = new QComboBox;
  modalitySwitcher->setToolTip(QString("Select the Image Modality"));
  for (int i = 0; i < IMAGE_TYPE_FEATURES + 1; i++)
  {
    modalitySwitcher->insertItem(i, ImageModalityString[i]);
  }

  //modalitySwitcher->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
  //int buttonWidth = QLabel().fontMetrics().width("----------------");
  //modalitySwitcher->setFixedWidth(buttonWidth);
  if (imgSubtype < modalitySwitcher->count())
  {
    modalitySwitcher->setCurrentIndex(imgSubtype);
  }
  m_imagesTable->setCellWidget(rowIndex, IMAGES_COLUMN_MODALITY, modalitySwitcher);

  QRadioButton* overlayRB = new QRadioButton();
  m_imagesTable->setCellWidget(rowIndex, IMAGES_COLUMN_OVERLAY, overlayRB);

  QTableWidgetItem *item2 = new QTableWidgetItem(fileinfo.fileName());
  item2->setData(Qt::UserRole, idstr.toStdString().c_str());

  QTablePushButton* cButton2 = new QTablePushButton;
  cButton2->setItem(item2);
  cButton2->setText(tr("X"));
  cButton2->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);

  connect(cButton2, SIGNAL(clickedInto(QTableWidgetItem*)), caller, SLOT(CloseImage(QTableWidgetItem*)));
  connect(modalitySwitcher, SIGNAL(currentIndexChanged(int)), this, SLOT(imageModalityChanged(int)));
  connect(overlayRB, SIGNAL(clicked()), caller, SLOT(overlayChanged()));
  connect(m_clearImagesBtn, SIGNAL(clicked()), caller, SLOT(CloseAllImages()));//TBD fix calling every time
  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));
}

void fImagesPanel::helpClicked()
{
  emit helpClicked_Interaction("Getting_Started.html");
}