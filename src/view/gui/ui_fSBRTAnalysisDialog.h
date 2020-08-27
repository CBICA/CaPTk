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

#ifndef ui_fSBRTAnalysisDialog_H
#define ui_fSBRTAnalysisDialog_H

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
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_fSBRTAnalysisDialog
{
public:
  QVBoxLayout *verticalLayout;
  QHBoxLayout *horizontalLayout_5;
  QLabel *modelDirLabel;
  QLineEdit *modelDirLineEdit;
  QPushButton *modelDirBtn;
  QHBoxLayout *horizontalLayout;
  QLabel *label;
  QHBoxLayout *horizontalLayout_2;
  QLabel *survivalLabel;
  QLineEdit *survivalValue;
  QHBoxLayout *horizontalLayout_4;
  QLabel *nodalfailureLabel;
  QLineEdit *nodalFailureValue;
  QHBoxLayout *horizontalLayout_3;
  QSpacerItem *horizontalSpacer;
  QPushButton *okButton;
  QSpacerItem *horizontalSpacer_2;

  void setupUi(QDialog *fSBRTAnalysisDialog)
  {
    if (fSBRTAnalysisDialog->objectName().isEmpty())
      fSBRTAnalysisDialog->setObjectName(QStringLiteral("fSBRTAnalysisDialog"));
    fSBRTAnalysisDialog->resize(397, 176);
    verticalLayout = new QVBoxLayout(fSBRTAnalysisDialog);
    verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
    horizontalLayout_5 = new QHBoxLayout();
    horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
    horizontalLayout_5->setContentsMargins(-1, -1, 0, -1);
    modelDirLabel = new QLabel(fSBRTAnalysisDialog);
    modelDirLabel->setObjectName(QStringLiteral("modelDirLabel"));

    horizontalLayout_5->addWidget(modelDirLabel);

    modelDirLineEdit = new QLineEdit(fSBRTAnalysisDialog);
    modelDirLineEdit->setObjectName(QStringLiteral("modelDirLineEdit"));
    QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(modelDirLineEdit->sizePolicy().hasHeightForWidth());
    modelDirLineEdit->setSizePolicy(sizePolicy);

    horizontalLayout_5->addWidget(modelDirLineEdit);

    modelDirBtn = new QPushButton(fSBRTAnalysisDialog);
    modelDirBtn->setObjectName(QStringLiteral("modelDirBtn"));

    horizontalLayout_5->addWidget(modelDirBtn);


    verticalLayout->addLayout(horizontalLayout_5);

    horizontalLayout = new QHBoxLayout();
    horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
    label = new QLabel(fSBRTAnalysisDialog);
    label->setObjectName(QStringLiteral("label"));
    label->setOpenExternalLinks(false);
    label->setTextInteractionFlags(Qt::LinksAccessibleByMouse);

    horizontalLayout->addWidget(label);


    verticalLayout->addLayout(horizontalLayout);

    horizontalLayout_2 = new QHBoxLayout();
    horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
    survivalLabel = new QLabel(fSBRTAnalysisDialog);
    survivalLabel->setObjectName(QStringLiteral("survivalLabel"));
    survivalLabel->setMinimumSize(QSize(173, 0));

    horizontalLayout_2->addWidget(survivalLabel);

    survivalValue = new QLineEdit(fSBRTAnalysisDialog);
    survivalValue->setObjectName(QStringLiteral("survivalValue"));
    survivalValue->setMaximumSize(QSize(16777215, 16777215));

    horizontalLayout_2->addWidget(survivalValue);


    verticalLayout->addLayout(horizontalLayout_2);

    horizontalLayout_4 = new QHBoxLayout();
    horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
    nodalfailureLabel = new QLabel(fSBRTAnalysisDialog);
    nodalfailureLabel->setObjectName(QStringLiteral("nodalfailureLabel"));

    horizontalLayout_4->addWidget(nodalfailureLabel);

    nodalFailureValue = new QLineEdit(fSBRTAnalysisDialog);
    nodalFailureValue->setObjectName(QStringLiteral("nodalFailureValue"));

    horizontalLayout_4->addWidget(nodalFailureValue);


    verticalLayout->addLayout(horizontalLayout_4);

    horizontalLayout_3 = new QHBoxLayout();
    horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
    horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    horizontalLayout_3->addItem(horizontalSpacer);

    okButton = new QPushButton(fSBRTAnalysisDialog);
    okButton->setObjectName(QStringLiteral("okButton"));

    horizontalLayout_3->addWidget(okButton);

    horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

    horizontalLayout_3->addItem(horizontalSpacer_2);


    verticalLayout->addLayout(horizontalLayout_3);


    retranslateUi(fSBRTAnalysisDialog);

    QMetaObject::connectSlotsByName(fSBRTAnalysisDialog);
  } // setupUi

  void retranslateUi(QDialog *fSBRTAnalysisDialog)
  {
    fSBRTAnalysisDialog->setWindowTitle(QApplication::translate("fSBRTAnalysisDialog", "Prognostic Modeling", nullptr));
    modelDirLabel->setText(QApplication::translate("fSBRTAnalysisDialog", "Model Directory: ", nullptr));
    modelDirBtn->setText(QApplication::translate("fSBRTAnalysisDialog", "Browse", nullptr));
    label->setText(QApplication::translate("fSBRTAnalysisDialog", "<html><head/><body><p>Please <a href=\"www.qtcentre.org\"><span style=\" text-decoration: underline; color:#0000ff;\">download</span></a> pretrained model based on our study at PENN</p></body></html>", nullptr));
    survivalLabel->setText(QApplication::translate("fSBRTAnalysisDialog", "Predicted Risk (Survival)", nullptr));
    nodalfailureLabel->setText(QApplication::translate("fSBRTAnalysisDialog", "Predicted Risk (Nodal Failure) ", nullptr));
    okButton->setText(QApplication::translate("fSBRTAnalysisDialog", "Ok", nullptr));
  } // retranslateUi

};

namespace Ui {
  class fSBRTAnalysisDialog : public Ui_fSBRTAnalysisDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif 





