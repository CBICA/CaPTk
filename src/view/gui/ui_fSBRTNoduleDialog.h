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
// License Agreement: https://www.med.upenn.edu/cbica/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ui_fSBRTNoduleDialog_H
#define ui_fSBRTNoduleDialog_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_fSBRTNoduleDialog
{
public:
  QVBoxLayout *verticalLayout;
  QHBoxLayout *horizontalLayout_2;
  QLabel *seedImageLabel;
  QLineEdit *seedImageName;
  QPushButton *seedImageBrowseButton;
  QHBoxLayout *horizontalLayout_4;
  QLabel *labelValueLabel;
  QSpinBox *labelValueSpinBox;
  QSpacerItem *horizontalSpacer_3;
  QHBoxLayout *horizontalLayout_3;
  QSpacerItem *horizontalSpacer;
  QPushButton *okButton;
  QPushButton *cancelButton;
  QSpacerItem *horizontalSpacer_2;

  void setupUi(QDialog *fSBRTNoduleDialog)
  {
    if (fSBRTNoduleDialog->objectName().isEmpty())
      fSBRTNoduleDialog->setObjectName(QStringLiteral("fSBRTNoduleDialog"));
    fSBRTNoduleDialog->resize(397, 151);
    verticalLayout = new QVBoxLayout(fSBRTNoduleDialog);
    verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
    horizontalLayout_2 = new QHBoxLayout();
    horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
    seedImageLabel = new QLabel(fSBRTNoduleDialog);
    seedImageLabel->setObjectName(QStringLiteral("seedImageLabel"));
    seedImageLabel->setMinimumSize(QSize(130, 0));

    horizontalLayout_2->addWidget(seedImageLabel);

    seedImageName = new QLineEdit(fSBRTNoduleDialog);
    seedImageName->setObjectName(QStringLiteral("seedImageName"));
    seedImageName->setMaximumSize(QSize(173, 16777215));

    horizontalLayout_2->addWidget(seedImageName);

    seedImageBrowseButton = new QPushButton(fSBRTNoduleDialog);
    seedImageBrowseButton->setObjectName(QStringLiteral("seedImageBrowseButton"));

    horizontalLayout_2->addWidget(seedImageBrowseButton);


    verticalLayout->addLayout(horizontalLayout_2);

    horizontalLayout_4 = new QHBoxLayout();
    horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
    labelValueLabel = new QLabel(fSBRTNoduleDialog);
    labelValueLabel->setObjectName(QStringLiteral("labelValueLabel"));
    labelValueLabel->setMinimumSize(QSize(130, 0));

    horizontalLayout_4->addWidget(labelValueLabel);

    labelValueSpinBox = new QSpinBox(fSBRTNoduleDialog);
    labelValueSpinBox->setObjectName(QStringLiteral("labelValueSpinBox"));
    labelValueSpinBox->setMinimumSize(QSize(50, 0));
    labelValueSpinBox->setValue(2);

    horizontalLayout_4->addWidget(labelValueSpinBox);

    horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    horizontalLayout_4->addItem(horizontalSpacer_3);


    verticalLayout->addLayout(horizontalLayout_4);

    horizontalLayout_3 = new QHBoxLayout();
    horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
    horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    horizontalLayout_3->addItem(horizontalSpacer);

    okButton = new QPushButton(fSBRTNoduleDialog);
    okButton->setObjectName(QStringLiteral("okButton"));

    horizontalLayout_3->addWidget(okButton);

    cancelButton = new QPushButton(fSBRTNoduleDialog);
    cancelButton->setObjectName(QStringLiteral("cancelButton"));

    horizontalLayout_3->addWidget(cancelButton);

    horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    horizontalLayout_3->addItem(horizontalSpacer_2);


    verticalLayout->addLayout(horizontalLayout_3);


    retranslateUi(fSBRTNoduleDialog);

    QMetaObject::connectSlotsByName(fSBRTNoduleDialog);
  } // setupUi

  void retranslateUi(QDialog *fSBRTNoduleDialog)
  {
    fSBRTNoduleDialog->setWindowTitle(QApplication::translate("fSBRTNoduleDialog", "SBRT Nodule", nullptr));
    seedImageLabel->setText(QApplication::translate("fSBRTNoduleDialog", "Seed Image (Optional)", nullptr));
    seedImageBrowseButton->setText(QApplication::translate("fSBRTNoduleDialog", "Browse", nullptr));
    labelValueLabel->setText(QApplication::translate("fSBRTNoduleDialog", "Label Value (Optional)", nullptr));
    okButton->setText(QApplication::translate("fSBRTNoduleDialog", "Ok", nullptr));
    cancelButton->setText(QApplication::translate("fSBRTNoduleDialog", "Cancel", nullptr));
  } // retranslateUi

};

namespace Ui {
  class fSBRTNoduleDialog : public Ui_fSBRTNoduleDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif 





