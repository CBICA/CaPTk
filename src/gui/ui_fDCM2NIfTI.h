#ifndef ui_fDCM2NIfTIConverter_H
#define ui_fDCM2NIfTIConverter_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDialog>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QGroupBox>
#include <QtGui/QRadioButton>

QT_BEGIN_NAMESPACE

class ui_fDCM2NIfTIConverter
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;
  QLabel *inputImageLabel;
  QLineEdit *inputImageName;
  QPushButton *inputImageButton;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QFrame *line_3;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fDCM2NIfTIConverter)
  {

    if (fDCM2NIfTIConverter->objectName().isEmpty())
      fDCM2NIfTIConverter->setObjectName(QString::fromUtf8("fDCM2NIfTIConverter"));
    //fDCM2NIfTIConverter->setWindowModality(Qt::NonModal);
    fDCM2NIfTIConverter->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fDCM2NIfTIConverter->sizePolicy().hasHeightForWidth());
    fDCM2NIfTIConverter->setSizePolicy(sizePolicy);
    fDCM2NIfTIConverter->setMinimumSize(QSize(0, 0));

    //fDCM2NIfTIConverter->setModal(true);
    gridLayout = new QGridLayout(fDCM2NIfTIConverter);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // input 
    inputGroupBox = new QGroupBox(fDCM2NIfTIConverter);
    inputGroupBox->setTitle(QString::fromStdString("First Image in Series"));

    inputGridLayout = new QGridLayout(inputGroupBox);
    inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputImageLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputImageLabel->sizePolicy().hasHeightForWidth());
    inputImageLabel->setSizePolicy(sizePolicy);

    inputImageName = new QLineEdit(" ");
    inputImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputImageName->sizePolicy().hasHeightForWidth());
    inputImageName->setSizePolicy(sizePolicy);
    inputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputImageButton = new QPushButton(inputGroupBox);
    inputImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputImageButton->setText(QString("Browse"));

    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 1, 1, 1, 1);

    // output 
    outputGroupBox = new QGroupBox(fDCM2NIfTIConverter);
    outputGroupBox->setTitle(QString::fromStdString("Output File Name"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputImageLabel = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    outputImageLabel->setSizePolicy(sizePolicy);

    outputImageName = new QLineEdit(" ");
    outputImageName->setObjectName(QString::fromUtf8("outputImageName"));
    sizePolicy.setHeightForWidth(outputImageName->sizePolicy().hasHeightForWidth());
    outputImageName->setSizePolicy(sizePolicy);
    outputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputImageButton = new QPushButton(outputGroupBox);
    outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
    //outputImageButton->setIcon(ButtonIcon);
    //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    outputImageButton->setText(QString("Browse"));

    outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageName, 1, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 1, 1, 1, 1);

    // put the layout in perspective
    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fDCM2NIfTIConverter);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fDCM2NIfTIConverter);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fDCM2NIfTIConverter);

    QMetaObject::connectSlotsByName(fDCM2NIfTIConverter);
  } // setupUi

  void retranslateUi(QDialog *fDCM2NIfTIConverter)
  {
    fDCM2NIfTIConverter->setWindowTitle(QApplication::translate("fDCM2NIfTIConverter", "DICOM to NIfTI", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fDCM2NIfTIConverter", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fDCM2NIfTIConverter", "Cancel", 0, QApplication::UnicodeUTF8));
  } // retranslateUi

};

namespace Ui {
  class fDCM2NIfTIConverter : public ui_fDCM2NIfTIConverter {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDCM2NIfTIConverter_H