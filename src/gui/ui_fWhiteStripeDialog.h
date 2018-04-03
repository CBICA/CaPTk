#ifndef ui_fWhiteStripe_H
#define ui_fWhiteStripe_H

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
#include "QtGui/qcombobox.h"
#include <qvalidator.h>

QT_BEGIN_NAMESPACE

class ui_fWhiteStripe
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *optionsGroupBox;
  QGridLayout *optionsGridLayout;

  QLineEdit *options_twsWidth;
  QLabel *options_twsWidth_label;
  QLineEdit *options_sliceStartZ;
  QLabel *options_sliceStartZ_label;
  QLineEdit *options_sliceStopZ;
  QLabel *options_sliceStopZ_label;
  QLineEdit *options_tissuesMax;
  QLabel *options_tissuesMax_label;
  QLineEdit *options_smoothMax;
  QLabel *options_smoothMax_label;
  QLineEdit *options_smoothDelta;
  QLabel *options_smoothDelta_label;
  QLineEdit *options_histSize;
  QLabel *options_histSize_label;

  QLabel *options_T1Selector_label;
  QComboBox *options_T1Selector;
  //QRadioButton *options_T1selected;
  //QRadioButton *options_T2selected;

  QRadioButton *options_skullStrippedImage;
  QRadioButton *options_axialSlicing;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;


  void setupUi(QDialog *fWhiteStripeObj)

  {

    if (fWhiteStripeObj->objectName().isEmpty())
      fWhiteStripeObj->setObjectName(QString::fromUtf8("fWhiteStripeObj"));
    //fWhiteStripeObj->setWindowModality(Qt::NonModal);
    fWhiteStripeObj->resize(200, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fWhiteStripeObj->sizePolicy().hasHeightForWidth());
    fWhiteStripeObj->setSizePolicy(sizePolicy);
    fWhiteStripeObj->setMinimumSize(QSize(0, 0));

    //fWhiteStripeObj->setModal(true);
    gridLayout = new QGridLayout(fWhiteStripeObj);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    optionsGroupBox = new QGroupBox(fWhiteStripeObj);
    optionsGroupBox->setTitle(QString::fromStdString("Options"));

    optionsGridLayout = new QGridLayout(optionsGroupBox);
    optionsGridLayout->setObjectName(QString::fromUtf8("optionsGridLayout"));
    int gridRowCounter = 0;

    options_T1Selector_label = new QLabel(optionsGroupBox);
    options_T1Selector_label->setObjectName("Image Modality");
    options_T1Selector_label->setSizePolicy(sizePolicy);
    options_T1Selector = new QComboBox(optionsGroupBox);
    QSizePolicy sizePolicy8(QSizePolicy::Minimum, QSizePolicy::Fixed);
    options_T1Selector->setSizePolicy(sizePolicy);
    //options_T1Selector->setMaximumSize(QSize(100, 1000));
    options_T1Selector->insertItem(0, "T1 Image");
    options_T1Selector->insertItem(1, "T2 Image");
    options_T1Selector->setToolTip("Set the type of image being processed");

    //options_T1selected = new QRadioButton(optionsGroupBox);
    //options_T1selected->setText("T1 Image");

    //options_T2selected = new QRadioButton(optionsGroupBox);
    //options_T2selected->setText("T2 Image");

    optionsGridLayout->addWidget(options_T1Selector_label, gridRowCounter, 0);
    optionsGridLayout->addWidget(options_T1Selector, gridRowCounter, 1);

    options_twsWidth_label = new QLabel(optionsGroupBox);
    options_twsWidth_label->setObjectName("Wstripe_Radius");
    options_twsWidth_label->setSizePolicy(sizePolicy);
    options_twsWidth = new QLineEdit("0.05");
    options_twsWidth->setObjectName("twsWidth");
    options_twsWidth->setToolTip("WhiteStripe Radius (0.00 to 4.99)");
    options_twsWidth->setValidator(new QDoubleValidator(0.00, 4.99, 2, optionsGroupBox));
    options_twsWidth->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_twsWidth_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_twsWidth, gridRowCounter, 1, 1, 1);

    options_tissuesMax_label = new QLabel(optionsGroupBox);
    options_tissuesMax_label->setObjectName("Wstripe_tissuesMax");
    options_tissuesMax_label->setSizePolicy(sizePolicy);
    options_tissuesMax = new QLineEdit("5");
    options_tissuesMax->setObjectName("sliceStartZ");
    options_tissuesMax->setToolTip("Max Tissues (5 or 10 or 20)");
    options_tissuesMax->setValidator(new QIntValidator(5, 20, optionsGroupBox));
    options_tissuesMax->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_tissuesMax_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_tissuesMax, gridRowCounter, 1, 1, 1);

    options_smoothMax_label = new QLabel(optionsGroupBox);
    options_smoothMax_label->setObjectName("Wstripe_smoothMax");
    options_smoothMax_label->setSizePolicy(sizePolicy);
    options_smoothMax = new QLineEdit("10.0");
    options_smoothMax->setObjectName("smoothMax");
    options_smoothMax->setToolTip("Max smoothing (0.0 to 50.0)");
    options_smoothMax->setValidator(new QDoubleValidator(0.0, 50.0, 1, optionsGroupBox));
    options_smoothMax->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_smoothMax_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_smoothMax, gridRowCounter, 1, 1, 1);

    options_smoothDelta_label = new QLabel(optionsGroupBox);
    options_smoothDelta_label->setObjectName("Wstripe_smoothDelta");
    options_smoothDelta_label->setSizePolicy(sizePolicy);
    options_smoothDelta = new QLineEdit("0.5");
    options_smoothDelta->setObjectName("smoothDelta");
    options_smoothDelta->setToolTip("Smoothing Delta (0.0 to 10.0)");
    options_smoothDelta->setValidator(new QDoubleValidator(0.0, 10.0, 1, optionsGroupBox));
    options_smoothDelta->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_smoothDelta_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_smoothDelta, gridRowCounter, 1, 1, 1);

    options_histSize_label = new QLabel(optionsGroupBox);
    options_histSize_label->setObjectName("Wstripe_histSize");
    options_histSize_label->setSizePolicy(sizePolicy);
    options_histSize = new QLineEdit("2000");
    options_histSize->setObjectName("histSize");
    options_histSize->setToolTip("Histogram Bins (100 to 3000)");
    options_histSize->setValidator(new QIntValidator(100, 3000, optionsGroupBox));
    options_histSize->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_histSize_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_histSize, gridRowCounter, 1, 1, 1);

    options_skullStrippedImage = new QRadioButton(optionsGroupBox);
    options_skullStrippedImage->setText("Skull Stripped Image");

    options_axialSlicing = new QRadioButton(optionsGroupBox);
    options_axialSlicing->setText("Enable Axial Slicing");
    options_axialSlicing->setToolTip("Required to ensure WhiteStripe doesn't pick up too much bone");

    gridRowCounter++;
    optionsGridLayout->addWidget(options_skullStrippedImage, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_axialSlicing, gridRowCounter, 1, 1, 1);

    options_sliceStartZ_label = new QLabel(optionsGroupBox);
    options_sliceStartZ_label->setObjectName("Wstripe_sliceStartZ");
    options_sliceStartZ_label->setSizePolicy(sizePolicy);
    options_sliceStartZ = new QLineEdit("-1");
    options_sliceStartZ->setDisabled(true);
    options_sliceStartZ->setObjectName("sliceStartZ");
    options_sliceStartZ->setToolTip("Start Z-slice for cropping");
    options_sliceStartZ->setValidator(new QIntValidator(50, 100, optionsGroupBox));
    options_sliceStartZ->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_sliceStartZ_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_sliceStartZ, gridRowCounter, 1, 1, 1);

    options_sliceStopZ_label = new QLabel(optionsGroupBox);
    options_sliceStopZ_label->setObjectName("Wstripe_sliceStopZ");
    options_sliceStopZ_label->setSizePolicy(sizePolicy);
    options_sliceStopZ = new QLineEdit("-1");
    options_sliceStopZ->setDisabled(true);
    options_sliceStopZ->setObjectName("sliceStartZ");
    options_sliceStopZ->setToolTip("Start Z-slice for cropping");
    options_sliceStopZ->setValidator(new QIntValidator(100, 150, optionsGroupBox));
    options_sliceStopZ->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    gridRowCounter++;
    optionsGridLayout->addWidget(options_sliceStopZ_label, gridRowCounter, 0, 1, 1);
    optionsGridLayout->addWidget(options_sliceStopZ, gridRowCounter, 1, 1, 1);

    // output 
    outputGroupBox = new QGroupBox(fWhiteStripeObj);
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
    gridLayout->addWidget(optionsGroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fWhiteStripeObj);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    confirmButton->setText(QApplication::translate("fWhiteStripeObj", "Confirm", 0, QApplication::UnicodeUTF8));

    cancelButton = new QPushButton(fWhiteStripeObj);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    cancelButton->setText(QApplication::translate("fWhiteStripeObj", "Cancel", 0, QApplication::UnicodeUTF8));

    gridLayout->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fWhiteStripeObj);

    QMetaObject::connectSlotsByName(fWhiteStripeObj);
  } // setupUi

  void retranslateUi(QDialog *fWhiteStripeObj)
  {
    fWhiteStripeObj->setWindowTitle(QApplication::translate("fWhiteStripeObj", "WhiteStripe Normalization", 0, QApplication::UnicodeUTF8));
    options_twsWidth_label->setText("Wstripe_Radius");
    options_sliceStartZ_label->setText("Wstripe_sliceStartZ");
    options_sliceStopZ_label->setText("Wstripe_sliceStopZ");
    options_tissuesMax_label->setText("Wstripe_tissuesMax");
    options_smoothMax_label->setText("Wstripe_smoothMax");
    options_smoothDelta_label->setText("Wstripe_smoothDelta");
    options_histSize_label->setText("Wstripe_histSize");
  } // retranslateUi

};

namespace Ui {
  class fWhiteStripeObj : public ui_fWhiteStripe {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fWhiteStripeObj_H






