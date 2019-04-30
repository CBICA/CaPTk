#ifndef ui_fDeepMedicDialog_H
#define ui_fDeepMedicDialog_H

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

class ui_fDeepMedicDialog
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *modelGroupBox;
  QGridLayout *modelGridLayout;
  QLineEdit *modelDirName;
  QPushButton *modelImageButton;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputDirName;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QHBoxLayout * horizontalLayout;

  QGroupBox *modelSelectionGroupBox;
  QGridLayout *modelSelectionGridLayout;
  QRadioButton *brainTumorSegmentationButton;
  QRadioButton *skullStrippingButton;
  QRadioButton *customButton;

  void setupUi(QDialog *fDeepMedicDialog)
  {

    if (fDeepMedicDialog->objectName().isEmpty())
      fDeepMedicDialog->setObjectName(QString::fromUtf8("fDeepMedicDialog"));
    //fDeepMedicDialog->setWindowModality(Qt::NonModal);
    fDeepMedicDialog->resize(200, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fDeepMedicDialog->sizePolicy().hasHeightForWidth());
    fDeepMedicDialog->setSizePolicy(sizePolicy);
    fDeepMedicDialog->setMinimumSize(QSize(0, 0));

    //--------------------------------------------------------------------
    modelSelectionGroupBox = new QGroupBox(fDeepMedicDialog);
    modelSelectionGroupBox->setTitle(QString::fromStdString("Model Selector"));
    modelSelectionGridLayout = new QGridLayout(modelSelectionGroupBox);
    modelSelectionGridLayout->setObjectName(QString::fromUtf8("modelSelectionGridLayout"));

    //rdNewClassification = new QRadioButton(classificationGroupBox);
    //rdNewClassification->setObjectName(QString::fromUtf8("newClassification"));
    brainTumorSegmentationButton = new QRadioButton(modelSelectionGroupBox);
    brainTumorSegmentationButton->setObjectName(QString::fromUtf8("brainTumorSegmentationButton"));
    brainTumorSegmentationButton->setText("Brain Tumor Segmentation");
    skullStrippingButton = new QRadioButton(modelSelectionGroupBox);
    skullStrippingButton->setObjectName(QString::fromUtf8("skullStrippingButton"));
    skullStrippingButton->setText("Skull Stripping");
    customButton = new QRadioButton(modelSelectionGroupBox);
    customButton->setObjectName(QString::fromUtf8("customButton"));
    customButton->setText("Custom");

    modelSelectionGridLayout->addWidget(brainTumorSegmentationButton, 1, 1, 1, 1);
    modelSelectionGridLayout->addWidget(skullStrippingButton, 2, 1, 1, 1);
    modelSelectionGridLayout->addWidget(customButton, 3, 1, 1, 1);

    //fDeepMedicDialog->setModal(true);
    gridLayout = new QGridLayout(fDeepMedicDialog);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // existing model 
    modelGroupBox = new QGroupBox(fDeepMedicDialog);
    modelGroupBox->setTitle(QString::fromStdString("Existing Model Directory Name"));

    modelGridLayout = new QGridLayout(modelGroupBox);
    modelGridLayout->setObjectName(QString::fromUtf8("modelGridLayout"));

    auto currentModelDir = cbica::normPath(getCaPTkDataDir() + "/deepMedic/saved_models/brainSegmentation/");
    modelDirName = new QLineEdit(currentModelDir.c_str());
    modelDirName->setObjectName(QString::fromUtf8("modelDirName"));
    sizePolicy.setHeightForWidth(modelDirName->sizePolicy().hasHeightForWidth());
    modelDirName->setSizePolicy(sizePolicy);
    modelDirName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    modelDirName->setReadOnly(true);

    modelImageButton = new QPushButton(modelGroupBox);
    modelImageButton->setObjectName(QString::fromUtf8("modelDirButton"));
    modelImageButton->setText("Browse");
    modelImageButton->setToolTip("Location of 'modelConfig.txt' and 'model.ckpt'");

    modelGridLayout->addWidget(modelDirName, 0, 0, 1, 1);
    modelGridLayout->addWidget(modelImageButton, 0, 1, 1, 1);

    // output 
    outputGroupBox = new QGroupBox(fDeepMedicDialog);
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
    gridLayout->addWidget(modelSelectionGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(modelGroupBox, 2, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 3, 0, 1, 2);


    confirmButton = new QPushButton(fDeepMedicDialog);
    confirmButton->setObjectName(QString::fromUtf8("Confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fDeepMedicDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 4, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 4, 1, 1, 1);

    retranslateUi(fDeepMedicDialog);

    QMetaObject::connectSlotsByName(fDeepMedicDialog);
  } // setupUi

  void retranslateUi(QDialog *fDeepMedicDialog)
  {
    fDeepMedicDialog->setWindowTitle(QApplication::translate("fDeepMedicDialog", "DeepLearning-based Segmentation", 0));
    confirmButton->setText(QApplication::translate("fDeepMedicDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fDeepMedicDialog", "Cancel", 0));
  } // retranslateUi

};

namespace Ui {
  class fDeepMedicDialog : public ui_fDeepMedicDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDeepMedicDialog_H