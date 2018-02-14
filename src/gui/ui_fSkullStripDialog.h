#ifndef ui_fSkullStripper_H
#define ui_fSkullStripper_H

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

class ui_fSkullStripper
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

  void setupUi(QDialog *fSkullStripper)
  {
    std::string dataDir;

    if (QDir(QApplication::applicationDirPath() + "/../data/sri24/").exists()) // packaged binary
    {
      dataDir = QApplication::applicationDirPath().toStdString() + "/../data/sri24/";
    }
    else if (QDir(QApplication::applicationDirPath() + "/../../data/sri24/").exists()) // developer_mode
    {
      dataDir = QApplication::applicationDirPath().toStdString() + "/../../data/sri24/";
    }


    if (fSkullStripper->objectName().isEmpty())
      fSkullStripper->setObjectName(QString::fromUtf8("fSkullStripper"));
    //fSkullStripper->setWindowModality(Qt::NonModal);
    fSkullStripper->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fSkullStripper->sizePolicy().hasHeightForWidth());
    fSkullStripper->setSizePolicy(sizePolicy);
    fSkullStripper->setMinimumSize(QSize(0, 0));

    fSkullStripper->setModal(true);
    gridLayout = new QGridLayout(fSkullStripper);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // reference 
    referenceGroupBox = new QGroupBox(fSkullStripper);
    referenceGroupBox->setTitle(QString::fromStdString("Reference Image"));

    referenceGridLayout = new QGridLayout(referenceGroupBox);
    referenceGridLayout->setObjectName(QString::fromUtf8("referenceGridLayout"));

    referenceImageLabel = new QLabel(referenceGroupBox);
    sizePolicy.setHeightForWidth(referenceImageLabel->sizePolicy().hasHeightForWidth());
    referenceImageLabel->setSizePolicy(sizePolicy);

    referenceImageName = new QLineEdit(" ");
    referenceImageName->setObjectName(QString::fromUtf8("referenceImageName"));
    sizePolicy.setHeightForWidth(referenceImageName->sizePolicy().hasHeightForWidth());
    referenceImageName->setSizePolicy(sizePolicy);
    referenceImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    referenceImageName->setText((dataDir + "atlastImage.nii.gz").c_str());

    referenceImageButton = new QPushButton(referenceGroupBox);
    referenceImageButton->setObjectName(QString::fromUtf8("referenceImageButton"));
    referenceImageButton->setText(QString("Reference Atlas"));
    //referenceMaskButton->setToolTip(QString("Atlas based on which stripping is done"));

    referenceMaskName = new QLineEdit(" ");
    referenceMaskName->setObjectName(QString::fromUtf8("referenceImageName"));
    sizePolicy.setHeightForWidth(referenceMaskName->sizePolicy().hasHeightForWidth());
    referenceMaskName->setSizePolicy(sizePolicy);
    referenceMaskName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    referenceMaskName->setText((dataDir + "atlasMask.nii.gz").c_str());
    
    referenceMaskButton = new QPushButton(referenceGroupBox);
    referenceMaskButton->setObjectName(QString::fromUtf8("testSubjectsDirectoryButton"));
    referenceMaskButton->setText(QString("Reference Mask"));
    referenceMaskButton->setToolTip(QString("Atlas Mask based on which stripping is done"));

    referenceGridLayout->addWidget(referenceImageLabel, 0, 0, 1, 1);
    referenceGridLayout->addWidget(referenceImageName, 1, 0, 1, 1);
    referenceGridLayout->addWidget(referenceImageButton, 1, 1, 1, 1);
    referenceGridLayout->addWidget(referenceMaskName, 2, 0, 1, 1);
    referenceGridLayout->addWidget(referenceMaskButton, 2, 1, 1, 1);

    // input 
    inputGroupBox = new QGroupBox(fSkullStripper);
    inputGroupBox->setTitle(QString::fromStdString("Input Image"));

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
    outputGroupBox = new QGroupBox(fSkullStripper);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

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
    gridLayout->addWidget(referenceGroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fSkullStripper);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fSkullStripper);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fSkullStripper);

    QMetaObject::connectSlotsByName(fSkullStripper);
  } // setupUi

  void retranslateUi(QDialog *fSkullStripper)
  {
    fSkullStripper->setWindowTitle(QApplication::translate("fSkullStripper", "Skull Stripping", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fSkullStripper", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fSkullStripper", "Cancel", 0, QApplication::UnicodeUTF8));


  } // retranslateUi

};

namespace Ui {
  class fSkullStripper : public ui_fSkullStripper {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fSkullStripper_H






