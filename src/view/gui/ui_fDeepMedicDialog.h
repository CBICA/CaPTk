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

QT_BEGIN_NAMESPACE

class ui_fDeepMedicDialog
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputDirName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fDeepMedicDialog)
  {

    if (fDeepMedicDialog->objectName().isEmpty())
      fDeepMedicDialog->setObjectName(QString::fromUtf8("fDeepMedicDialog"));
    //fDeepMedicDialog->setWindowModality(Qt::NonModal);
    //fDeepMedicDialog->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fDeepMedicDialog->sizePolicy().hasHeightForWidth());
    fDeepMedicDialog->setSizePolicy(sizePolicy);
    fDeepMedicDialog->setMinimumSize(QSize(0, 0));

    //fDeepMedicDialog->setModal(true);
    gridLayout = new QGridLayout(fDeepMedicDialog);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // output 
    outputGroupBox = new QGroupBox(fDeepMedicDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory Name"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputImageLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    //outputImageLabel->setSizePolicy(sizePolicy);

    outputDirName = new QLineEdit(" ");
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
    //gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fDeepMedicDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fDeepMedicDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fDeepMedicDialog);

    QMetaObject::connectSlotsByName(fDeepMedicDialog);
  } // setupUi

  void retranslateUi(QDialog *fDeepMedicDialog)
  {
    fDeepMedicDialog->setWindowTitle(QApplication::translate("fDeepMedicDialog", "DeepMedic Segmentation", 0));
    confirmButton->setText(QApplication::translate("fDeepMedicDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fDeepMedicDialog", "Cancel", 0));
  } // retranslateUi

};

namespace Ui {
  class fDeepMedicDialog : public ui_fDeepMedicDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDeepMedicDialog_H