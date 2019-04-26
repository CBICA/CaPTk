#ifndef ui_fPreprocessingDialog_H
#define ui_fPreprocessingDialog_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QCheckBox>
// #include <QtGui/QDialog>
// #include <QtGui/QFrame>
// #include <QtGui/QGridLayout>
// #include <QtGui/QHBoxLayout>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QLineEdit>
// #include <QtGui/QPushButton>
// #include <QtGui/QSpacerItem>
// #include <QtGui/QSpinBox>
// NEW CHANGES
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

class ui_fPreprocessingDialog
{
public:
  QGridLayout *gridLayout_3;
  QVBoxLayout * verticalLayout;
  QGroupBox *classificationGroupBox;
  QGridLayout * classificationGridLayout;
  QGroupBox *outputGroupBox;
  QRadioButton * rdNewClassification;
  QRadioButton * rdExistingClassification;
  QRadioButton *rdCreateModel;
  QLineEdit * svmModelFileName;
  QPushButton * svmModelButton;
  QLineEdit * testSubjectsDirectoryName;
  QPushButton * testSubjectsDirectoryButton;


  QLineEdit * existingMaskDirectoryName;
  QPushButton * existingMasksButton;

  QGridLayout *outputGridLayout;
  QLabel	*outputDirectoryLabel;
  QLineEdit *outputDirectoryName;
  QPushButton *outputDirectoryButton;

  QPushButton * confirmButton;
  QPushButton * cancelButton;

  QFrame *line_3;
  QCheckBox * cbT1Data;
  QCheckBox * cbT2Data;
  QCheckBox * cbT1ceData;
  QCheckBox * cbT2FlairData;
  QCheckBox * cbDTIData;
  QCheckBox * cbPerfData;
  QCheckBox * cbDistanceData;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fPreprocessingDialog)
  {

    if (fPreprocessingDialog->objectName().isEmpty())
      fPreprocessingDialog->setObjectName(QString::fromUtf8("fPreprocessingDialog"));

    fPreprocessingDialog->setWindowModality(Qt::ApplicationModal);
    fPreprocessingDialog->resize(200, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPreprocessingDialog->sizePolicy().hasHeightForWidth());
    fPreprocessingDialog->setSizePolicy(sizePolicy);
    fPreprocessingDialog->setMinimumSize(QSize(0, 0));
    QFont font;
    font.setFamily(QString::fromUtf8("Calibri"));
    fPreprocessingDialog->setFont(font);
    //fPreprocessingDialog->setModal(true);
    gridLayout_3 = new QGridLayout(fPreprocessingDialog);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //--------------------------------------------------------------------
    classificationGroupBox = new QGroupBox(fPreprocessingDialog);
    classificationGroupBox->setTitle(QString::fromStdString("Classification"));
    classificationGridLayout = new QGridLayout(classificationGroupBox);
    classificationGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));

    rdNewClassification = new QRadioButton(classificationGroupBox);
    rdNewClassification->setObjectName(QString::fromUtf8("newClassification"));
    rdExistingClassification = new QRadioButton(classificationGroupBox);
    rdExistingClassification->setObjectName(QString::fromUtf8("existingClassification"));
    rdCreateModel = new QRadioButton(classificationGroupBox);
    rdCreateModel->setObjectName(QString::fromUtf8("createModel"));


    svmModelFileName = new QLineEdit(classificationGroupBox);
    svmModelFileName->setObjectName(QString::fromUtf8("svmModeFileName"));
    QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy13.setHorizontalStretch(0);
    sizePolicy13.setVerticalStretch(0);
    sizePolicy13.setHeightForWidth(svmModelFileName->sizePolicy().hasHeightForWidth());
    svmModelFileName->setSizePolicy(sizePolicy13);
    svmModelFileName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    testSubjectsDirectoryName = new QLineEdit(classificationGroupBox);
    testSubjectsDirectoryName->setObjectName(QString::fromUtf8("testSubjectsDirectoryName"));
    sizePolicy13.setHeightForWidth(testSubjectsDirectoryName->sizePolicy().hasHeightForWidth());
    testSubjectsDirectoryName->setSizePolicy(sizePolicy13);
    testSubjectsDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    existingMaskDirectoryName = new QLineEdit(classificationGroupBox);
    existingMaskDirectoryName->setObjectName(QString::fromUtf8("existingMaskDirectoryName"));
    sizePolicy13.setHeightForWidth(existingMaskDirectoryName->sizePolicy().hasHeightForWidth());
    existingMaskDirectoryName->setSizePolicy(sizePolicy13);
    existingMaskDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    svmModelButton = new QPushButton(classificationGroupBox);
    svmModelButton->setObjectName(QString::fromUtf8("svmModelButton"));
    QPixmap pixmap((getCaPTkDataDir() + "/icons/open.png").c_str());
    QIcon ButtonIcon(pixmap);
    svmModelButton->setIcon(ButtonIcon);
    svmModelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    testSubjectsDirectoryButton = new QPushButton(classificationGroupBox);
    testSubjectsDirectoryButton->setObjectName(QString::fromUtf8("testSubjectsDirectoryButton"));
    testSubjectsDirectoryButton->setIcon(ButtonIcon);
    testSubjectsDirectoryButton->setIconSize(QSize(20, 20));


    existingMasksButton = new QPushButton(classificationGroupBox);
    existingMasksButton->setObjectName(QString::fromUtf8("existingMasksButton"));
    existingMasksButton->setIcon(ButtonIcon);
    existingMasksButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent


    classificationGridLayout->addWidget(rdNewClassification, 0, 0, 1, 1);
    classificationGridLayout->addWidget(rdExistingClassification, 1, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName, 1, 1, 1, 1);
    classificationGridLayout->addWidget(svmModelButton, 1, 2, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryName, 2, 1, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryButton, 2, 2, 1, 1);

    classificationGridLayout->addWidget(rdCreateModel, 3, 0, 1, 1);
    classificationGridLayout->addWidget(existingMaskDirectoryName, 3, 1, 1, 1);
    classificationGridLayout->addWidget(existingMasksButton, 3, 2, 1, 1);

    line_3 = new QFrame(classificationGroupBox);
    line_3->setFrameStyle(QFrame::HLine);
    line_3->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    line_3->setFrameShadow(QFrame::Sunken);

    classificationGridLayout->addWidget(line_3, 4, 1, 1, 2);

    cbT1Data = new QCheckBox(classificationGroupBox);
    cbT1Data->setObjectName(QString::fromUtf8("cbT1Data"));
    cbT2Data = new QCheckBox(classificationGroupBox);
    cbT2Data->setObjectName(QString::fromUtf8("cbT2Data"));
    cbT1ceData = new QCheckBox(classificationGroupBox);
    cbT1ceData->setObjectName(QString::fromUtf8("cbT1ceData"));
    cbT2FlairData = new QCheckBox(classificationGroupBox);
    cbT2FlairData->setObjectName(QString::fromUtf8("cbT2FlairData"));
    cbDTIData = new QCheckBox(classificationGroupBox);
    cbDTIData->setObjectName(QString::fromUtf8("cbDTIData"));
    cbPerfData = new QCheckBox(classificationGroupBox);
    cbPerfData->setObjectName(QString::fromUtf8("cbPerfData"));
    cbDistanceData = new QCheckBox(classificationGroupBox);
    cbDistanceData->setObjectName(QString::fromUtf8("cbDistanceData"));


    horizontalLayout = new QHBoxLayout();
    horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));

    horizontalLayout->addWidget(cbT1Data);
    horizontalLayout->addWidget(cbT2Data);
    horizontalLayout->addWidget(cbT1ceData);
    horizontalLayout->addWidget(cbT2FlairData);
    horizontalLayout->addWidget(cbDTIData);
    horizontalLayout->addWidget(cbPerfData);
    horizontalLayout->addWidget(cbDistanceData);

    classificationGridLayout->addLayout(horizontalLayout, 5, 0, 1, 3);
    //--------------------------output-------------------------------------------------------


    outputGroupBox = new QGroupBox(fPreprocessingDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputDirectoryLabel = new QLabel(outputGroupBox);
    sizePolicy13.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    outputDirectoryLabel->setSizePolicy(sizePolicy13);

    outputDirectoryName = new QLineEdit(outputGroupBox);
    outputDirectoryName->setObjectName(QString::fromUtf8("outputDirectoryName"));
    sizePolicy13.setHeightForWidth(outputDirectoryName->sizePolicy().hasHeightForWidth());
    outputDirectoryName->setSizePolicy(sizePolicy13);
    outputDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputDirectoryButton = new QPushButton(outputGroupBox);
    outputDirectoryButton->setObjectName(QString::fromUtf8("outputDirectoryButton"));
    outputDirectoryButton->setIcon(ButtonIcon);
    outputDirectoryButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    outputGridLayout->addWidget(outputDirectoryLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirectoryName, 1, 0, 1, 1);
    outputGridLayout->addWidget(outputDirectoryButton, 1, 1, 1, 1);

    gridLayout_3->addWidget(classificationGroupBox, 0, 0, 1, 2);
    gridLayout_3->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fPreprocessingDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    confirmButton->setIcon(ButtonIcon);
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fPreprocessingDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    cancelButton->setIcon(ButtonIcon);
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout_3->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout_3->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fPreprocessingDialog);

    QMetaObject::connectSlotsByName(fPreprocessingDialog);
  } // setupUi

  void retranslateUi(QDialog *fPreprocessingDialog)
  {
    // fPreprocessingDialog->setWindowTitle(QApplication::translate("fPreprocessingDialog", "Glioblastoma Infiltration Index", 0, QApplication::UnicodeUTF8));
    // rdNewClassification->setText(QApplication::translate("fPreprocessingDialog", "Subject based", 0, QApplication::UnicodeUTF8));
    // rdExistingClassification->setText(QApplication::translate("fPreprocessingDialog", "Use existing model", 0, QApplication::UnicodeUTF8));
    // rdCreateModel->setText(QApplication::translate("fPreprocessingDialog", "Train new model", 0, QApplication::UnicodeUTF8));
    // confirmButton->setText(QApplication::translate("fPreprocessingDialog", "Confirm", 0, QApplication::UnicodeUTF8));
    // cancelButton->setText(QApplication::translate("fPreprocessingDialog", "Cancel", 0, QApplication::UnicodeUTF8));
    //
    // cbT1Data->setText(QApplication::translate("fPreprocessingDialog", "T1", 0, QApplication::UnicodeUTF8));
    // cbT1ceData->setText(QApplication::translate("fPreprocessingDialog", "T1ce", 0, QApplication::UnicodeUTF8));
    // cbT2Data->setText(QApplication::translate("fPreprocessingDialog", "T2", 0, QApplication::UnicodeUTF8));
    // cbT2FlairData->setText(QApplication::translate("fPreprocessingDialog", "Flair", 0, QApplication::UnicodeUTF8));
    // cbDTIData->setText(QApplication::translate("fPreprocessingDialog", "DTI", 0, QApplication::UnicodeUTF8));
    // cbPerfData->setText(QApplication::translate("fPreprocessingDialog", "Perfusion", 0, QApplication::UnicodeUTF8));
    // cbDistanceData->setText(QApplication::translate("fPreprocessingDialog", "Distance", 0, QApplication::UnicodeUTF8));
    // NEW CHANGES
    fPreprocessingDialog->setWindowTitle(QApplication::translate("fPreprocessingDialog", "Glioblastoma Infiltration Index", 0));
    rdNewClassification->setText(QApplication::translate("fPreprocessingDialog", "Subject based", 0));
    rdExistingClassification->setText(QApplication::translate("fPreprocessingDialog", "Use existing model", 0));
    rdCreateModel->setText(QApplication::translate("fPreprocessingDialog", "Train new model", 0));
    confirmButton->setText(QApplication::translate("fPreprocessingDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fPreprocessingDialog", "Cancel", 0));

    cbT1Data->setText(QApplication::translate("fPreprocessingDialog", "T1", 0));
    cbT1ceData->setText(QApplication::translate("fPreprocessingDialog", "T1-Gd", 0));
    cbT2Data->setText(QApplication::translate("fPreprocessingDialog", "T2", 0));
    cbT2FlairData->setText(QApplication::translate("fPreprocessingDialog", "Flair", 0));
    cbDTIData->setText(QApplication::translate("fPreprocessingDialog", "DTI", 0));
    cbPerfData->setText(QApplication::translate("fPreprocessingDialog", "Perfusion", 0));
    cbDistanceData->setText(QApplication::translate("fPreprocessingDialog", "Distance", 0));
  } // retranslateUi

};

namespace Ui {
  class fPreprocessingDialog : public ui_fPreprocessingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPreprocessingDialog_H
