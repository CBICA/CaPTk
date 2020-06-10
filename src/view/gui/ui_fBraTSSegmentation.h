#ifndef ui_fBraTSSegmentation_H
#define ui_fBraTSSegmentation_H

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

class ui_fBraTSSegmentation
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *inputT1CEGroupBox;
  QGridLayout *inputT1CEGridLayout;

  QLabel *inputT1CEImageLabel;
  QLineEdit *inputT1CEImageName;
  QPushButton *inputT1CEImageButton;

  QGroupBox *inputT1GroupBox;
  QGridLayout *inputT1GridLayout;

  QLabel *inputT1ImageLabel;
  QLineEdit *inputT1ImageName;
  QPushButton *inputT1ImageButton;

  QGroupBox *inputT2GroupBox;
  QGridLayout *inputT2GridLayout;

  QLabel *inputT2ImageLabel;
  QLineEdit *inputT2ImageName;
  QPushButton *inputT2ImageButton;

  QGroupBox *inputFLGroupBox;
  QGridLayout *inputFLGridLayout;

  QLabel *inputFLImageLabel;
  QLineEdit *inputFLImageName;
  QPushButton *inputFLImageButton;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QFrame *line_3;

  QHBoxLayout * horizontalLayout;

  //QLabel *longRunningWarning;

  void setupUi(QDialog *fBraTSSegmentation)
  {

    if (fBraTSSegmentation->objectName().isEmpty())
      fBraTSSegmentation->setObjectName(QString::fromUtf8("fBraTSSegmentation"));
    //fBraTSSegmentation->setWindowModality(Qt::NonModal);
    fBraTSSegmentation->resize(400, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fBraTSSegmentation->sizePolicy().hasHeightForWidth());
    fBraTSSegmentation->setSizePolicy(sizePolicy);
    fBraTSSegmentation->setMinimumSize(QSize(0, 0));

    //fBraTSSegmentation->setModal(true);
    gridLayout = new QGridLayout(fBraTSSegmentation);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // input t1ce
    inputT1CEGroupBox = new QGroupBox(fBraTSSegmentation);
    inputT1CEGroupBox->setTitle(QString::fromStdString("Input T1CE Image"));

    inputT1CEGridLayout = new QGridLayout(inputT1CEGroupBox);
    inputT1CEGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputT1CEImageLabel = new QLabel(inputT1CEGroupBox);
    sizePolicy.setHeightForWidth(inputT1CEImageLabel->sizePolicy().hasHeightForWidth());
    inputT1CEImageLabel->setSizePolicy(sizePolicy);

    inputT1CEImageName = new QLineEdit("");
    inputT1CEImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputT1CEImageName->sizePolicy().hasHeightForWidth());
    inputT1CEImageName->setSizePolicy(sizePolicy);
    inputT1CEImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputT1CEImageButton = new QPushButton(inputT1CEGroupBox);
    inputT1CEImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputT1CEImageButton->setText(QString("Browse"));

    inputT1CEGridLayout->addWidget(inputT1CEImageLabel, 0, 0, 1, 1);
    inputT1CEGridLayout->addWidget(inputT1CEImageName, 1, 0, 1, 1);
    inputT1CEGridLayout->addWidget(inputT1CEImageButton, 1, 1, 1, 1);

    // input t1
    inputT1GroupBox = new QGroupBox(fBraTSSegmentation);
    inputT1GroupBox->setTitle(QString::fromStdString("Input T1 Image"));

    inputT1GridLayout = new QGridLayout(inputT1GroupBox);
    inputT1GridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputT1ImageLabel = new QLabel(inputT1GroupBox);
    sizePolicy.setHeightForWidth(inputT1ImageLabel->sizePolicy().hasHeightForWidth());
    inputT1ImageLabel->setSizePolicy(sizePolicy);

    inputT1ImageName = new QLineEdit("");
    inputT1ImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputT1ImageName->sizePolicy().hasHeightForWidth());
    inputT1ImageName->setSizePolicy(sizePolicy);
    inputT1ImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputT1ImageButton = new QPushButton(inputT1GroupBox);
    inputT1ImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputT1ImageButton->setText(QString("Browse"));

    inputT1GridLayout->addWidget(inputT1ImageLabel, 0, 0, 1, 1);
    inputT1GridLayout->addWidget(inputT1ImageName, 1, 0, 1, 1);
    inputT1GridLayout->addWidget(inputT1ImageButton, 1, 1, 1, 1);

    // input t1
    inputT2GroupBox = new QGroupBox(fBraTSSegmentation);
    inputT2GroupBox->setTitle(QString::fromStdString("Input T2 Image"));

    inputT2GridLayout = new QGridLayout(inputT2GroupBox);
    inputT2GridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputT2ImageLabel = new QLabel(inputT2GroupBox);
    sizePolicy.setHeightForWidth(inputT2ImageLabel->sizePolicy().hasHeightForWidth());
    inputT2ImageLabel->setSizePolicy(sizePolicy);

    inputT2ImageName = new QLineEdit("");
    inputT2ImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputT2ImageName->sizePolicy().hasHeightForWidth());
    inputT2ImageName->setSizePolicy(sizePolicy);
    inputT2ImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputT2ImageButton = new QPushButton(inputT2GroupBox);
    inputT2ImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputT2ImageButton->setText(QString("Browse"));

    inputT2GridLayout->addWidget(inputT2ImageLabel, 0, 0, 1, 1);
    inputT2GridLayout->addWidget(inputT2ImageName, 1, 0, 1, 1);
    inputT2GridLayout->addWidget(inputT2ImageButton, 1, 1, 1, 1);

    // input t1
    inputFLGroupBox = new QGroupBox(fBraTSSegmentation);
    inputFLGroupBox->setTitle(QString::fromStdString("Input FL Image"));

    inputFLGridLayout = new QGridLayout(inputFLGroupBox);
    inputFLGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputFLImageLabel = new QLabel(inputFLGroupBox);
    sizePolicy.setHeightForWidth(inputFLImageLabel->sizePolicy().hasHeightForWidth());
    inputFLImageLabel->setSizePolicy(sizePolicy);

    inputFLImageName = new QLineEdit("");
    inputFLImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputFLImageName->sizePolicy().hasHeightForWidth());
    inputFLImageName->setSizePolicy(sizePolicy);
    inputFLImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputFLImageButton = new QPushButton(inputFLGroupBox);
    inputFLImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputFLImageButton->setText(QString("Browse"));

    inputFLGridLayout->addWidget(inputFLImageLabel, 0, 0, 1, 1);
    inputFLGridLayout->addWidget(inputFLImageName, 1, 0, 1, 1);
    inputFLGridLayout->addWidget(inputFLImageButton, 1, 1, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fBraTSSegmentation);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputImageLabel = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    outputImageLabel->setSizePolicy(sizePolicy);

    outputImageName = new QLineEdit("");
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
    gridLayout->addWidget(inputT1CEGroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(inputT1GroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(inputT2GroupBox, 2, 0, 1, 2);
    gridLayout->addWidget(inputFLGroupBox, 3, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 4, 0, 1, 2);


    confirmButton = new QPushButton(fBraTSSegmentation);
    confirmButton->setObjectName(QString::fromUtf8("Confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fBraTSSegmentation);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout->addWidget(confirmButton, 5, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 5, 1, 1, 1);

    retranslateUi(fBraTSSegmentation);

    QMetaObject::connectSlotsByName(fBraTSSegmentation);
  } // setupUi

  void retranslateUi(QDialog *fBraTSSegmentation)
  {
    fBraTSSegmentation->setWindowTitle(QApplication::translate("fBraTSSegmentation", "BraTS Pipeline", 0));
    confirmButton->setText(QApplication::translate("fBraTSSegmentation", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fBraTSSegmentation", "Cancel", 0));

  } // retranslateUi

};

namespace Ui {
  class fBraTSSegmentation : public ui_fBraTSSegmentation {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fBraTSSegmentation_H
