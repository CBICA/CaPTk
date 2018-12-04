#ifndef ui_fPCAEstimator_H
#define ui_fPCAEstimator_H

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
#include "CapTkGUIUtils.h"
#include <QDir>

QT_BEGIN_NAMESPACE

class ui_fPCAEstimator
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *referenceGroupBox;
  QGridLayout *referenceGridLayout;
  QLabel *referenceImageLabel;
  QLineEdit *referenceImageName;
  QPushButton *referenceImageButton;
  QLineEdit *referenceMaskName;
  QPushButton *referenceMaskButton;

  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;
  QLabel *inputImageLabel;
  QLineEdit *inputImageName;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QFrame *line_3;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fPCAEstimator)
  {
    std::string dataDir;

    if (QDir((getCaPTkDataDir() + "/sri24/").c_str()).exists()) // packaged binary
    {
      dataDir = getCaPTkDataDir() + "/sri24/";
    }
    else if (QDir(QApplication::applicationDirPath() + "/../../data/sri24/").exists()) // developer_mode
    {
      dataDir = captk_currentApplicationPath + "/../../data/sri24/";
    }


    if (fPCAEstimator->objectName().isEmpty())
      fPCAEstimator->setObjectName(QString::fromUtf8("fPCAEstimator"));
    //fPCAEstimator->setWindowModality(Qt::NonModal);
    fPCAEstimator->resize(400, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPCAEstimator->sizePolicy().hasHeightForWidth());
    fPCAEstimator->setSizePolicy(sizePolicy);
    fPCAEstimator->setMinimumSize(QSize(0, 0));

    //fPCAEstimator->setModal(true);
    gridLayout = new QGridLayout(fPCAEstimator);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // input
    inputGroupBox = new QGroupBox(fPCAEstimator);
    inputGroupBox->setTitle(QString::fromStdString("Input Data"));

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


    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 1, 0, 1, 1);
    //inputGridLayout->addWidget(inputImageButton, 1, 1, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fPCAEstimator);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputImageLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    //outputImageLabel->setSizePolicy(sizePolicy);

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

    //outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);

    // put the layout in perspective
    gridLayout->addWidget(inputGroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fPCAEstimator);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fPCAEstimator);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fPCAEstimator);

    QMetaObject::connectSlotsByName(fPCAEstimator);
  } // setupUi

  void retranslateUi(QDialog *fPCAEstimator)
  {
   //  fPCAEstimator->setWindowTitle(QApplication::translate("fPCAEstimator", "Principal component analysis", 0, QApplication::UnicodeUTF8));
   //  confirmButton->setText(QApplication::translate("fPCAEstimator", "Confirm", 0, QApplication::UnicodeUTF8));
	// inputImageLabel->setText(QApplication::translate("fPCAEstimator", "Number of PCAs", 0, QApplication::UnicodeUTF8));
   //  cancelButton->setText(QApplication::translate("fPCAEstimator", "Cancel", 0, QApplication::UnicodeUTF8));
   // NEW CHANGES
   fPCAEstimator->setWindowTitle(QApplication::translate("fPCAEstimator", "Principal component analysis", 0));
   confirmButton->setText(QApplication::translate("fPCAEstimator", "Confirm", 0));
    inputImageLabel->setText(QApplication::translate("fPCAEstimator", "Number of Principal Components", 0));
   cancelButton->setText(QApplication::translate("fPCAEstimator", "Cancel", 0));
  } // retranslateUi

};

namespace Ui {
  class fPCAEstimator : public ui_fPCAEstimator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPCAEstimator_H
