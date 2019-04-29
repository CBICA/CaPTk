#ifndef ui_fTexturePipelineDialog_H
#define ui_fTexturePipelineDialog_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QGroupBox>
#include <QtWidgets/QRadioButton>

#include "CaPTkGUIUtils.h"

QT_BEGIN_NAMESPACE

class ui_fTexturePipelineDialog
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputDirName;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QHBoxLayout * horizontalLayout;


  void setupUi(QDialog *fTexturePipelineDialog)
  {

    if (fTexturePipelineDialog->objectName().isEmpty())
      fTexturePipelineDialog->setObjectName(QString::fromUtf8("fTexturePipelineDialog"));
    //fTexturePipelineDialog->setWindowModality(Qt::NonModal);
    fTexturePipelineDialog->resize(200, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fTexturePipelineDialog->sizePolicy().hasHeightForWidth());
    fTexturePipelineDialog->setSizePolicy(sizePolicy);
    fTexturePipelineDialog->setMinimumSize(QSize(0, 0));

    //fTexturePipelineDialog->setModal(true);
    gridLayout = new QGridLayout(fTexturePipelineDialog);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // output 
    outputGroupBox = new QGroupBox(fTexturePipelineDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory Name"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputImageLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    //outputImageLabel->setSizePolicy(sizePolicy);

    outputDirName = new QLineEdit("");
    outputDirName->setObjectName(QString::fromUtf8("outputDirName"));
    sizePolicy.setHeightForWidth(outputDirName->sizePolicy().hasHeightForWidth());
    outputDirName->setSizePolicy(sizePolicy);
    outputDirName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputImageButton = new QPushButton(outputGroupBox);
    outputImageButton->setObjectName(QString::fromUtf8("outputDirButton"));
    //outputImageButton->setIcon(ButtonIcon);
    //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    outputImageButton->setText(QString("Browse"));

    //outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);

    // put the layout in perspective
    gridLayout->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fTexturePipelineDialog);
    confirmButton->setObjectName(QString::fromUtf8("Confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fTexturePipelineDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 4, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 4, 1, 1, 1);

    retranslateUi(fTexturePipelineDialog);

    QMetaObject::connectSlotsByName(fTexturePipelineDialog);
  } // setupUi

  void retranslateUi(QDialog *fTexturePipelineDialog)
  {
    fTexturePipelineDialog->setWindowTitle(QApplication::translate("fTexturePipelineDialog", "Breast Texture Feature Pipeline", 0));
    confirmButton->setText(QApplication::translate("fTexturePipelineDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fTexturePipelineDialog", "Cancel", 0));
  } // retranslateUi

};

namespace Ui {
  class fTexturePipelineDialog : public ui_fTexturePipelineDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fTexturePipelineDialog_H