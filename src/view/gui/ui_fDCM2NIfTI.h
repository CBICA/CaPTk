#ifndef ui_fDCM2NIfTIConverter_H
#define ui_fDCM2NIfTIConverter_H

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

QT_BEGIN_NAMESPACE

class ui_fDCM2NIfTIConverter
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;
  QLabel *inputDirLabel;
  QLineEdit *inputDirName;
  QPushButton *inputDirButton;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputDirName;
  QLabel *outputDirLabel;
  QPushButton *outputDirButton;

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
    inputGroupBox->setTitle(QString::fromStdString("Input Directory"));

    inputGridLayout = new QGridLayout(inputGroupBox);
    inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputDirLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputDirLabel->sizePolicy().hasHeightForWidth());
    inputDirLabel->setSizePolicy(sizePolicy);

    inputDirName = new QLineEdit("");
    inputDirName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputDirName->sizePolicy().hasHeightForWidth());
    inputDirName->setSizePolicy(sizePolicy);
    inputDirName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputDirButton = new QPushButton(inputGroupBox);
    inputDirButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputDirButton->setText(QString("Browse"));

    inputGridLayout->addWidget(inputDirLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputDirName, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputDirButton, 1, 1, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fDCM2NIfTIConverter);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputDirLabel = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(outputDirLabel->sizePolicy().hasHeightForWidth());
    outputDirLabel->setSizePolicy(sizePolicy);

    outputDirName = new QLineEdit("");
    outputDirName->setObjectName(QString::fromUtf8("outputDirName"));
    sizePolicy.setHeightForWidth(outputDirName->sizePolicy().hasHeightForWidth());
    outputDirName->setSizePolicy(sizePolicy);
    outputDirName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputDirButton = new QPushButton(outputGroupBox);
    outputDirButton->setObjectName(QString::fromUtf8("outputImageButton"));
    //outputImageButton->setIcon(ButtonIcon);
    //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent
    outputDirButton->setText(QString("Browse"));

    outputGridLayout->addWidget(outputDirLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirName, 1, 0, 1, 1);
    outputGridLayout->addWidget(outputDirButton, 1, 1, 1, 1);

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
    fDCM2NIfTIConverter->setWindowTitle(QApplication::translate("fDCM2NIfTIConverter", "DICOM to NIfTI", 0));
    confirmButton->setText(QApplication::translate("fDCM2NIfTIConverter", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fDCM2NIfTIConverter", "Cancel", 0));
  } // retranslateUi

};

namespace Ui {
  class fDCM2NIfTIConverter : public ui_fDCM2NIfTIConverter {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDCM2NIfTIConverter_H