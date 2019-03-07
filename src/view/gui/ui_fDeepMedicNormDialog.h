#ifndef ui_fDeepMedicNormalizer_H
#define ui_fDeepMedicNormalizer_H

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

class ui_fDeepMedicNormalizer
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *maskGroupBox;
  QGridLayout *maskGridLayout;
  QLabel *maskImageLabel;
  QLineEdit *maskImageName;
  QPushButton *maskImageButton;

  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;
  QLabel *inputImageLabel;
  QLineEdit *inputImageName;
  QPushButton *inputImageButton;

  QGroupBox *optionsGroupBox;
  QGridLayout *optionsGridLayout;
  QLabel *options_quantileLowerLabel;
  QLineEdit *options_quantileLowerName;
  QLabel *options_quantileUpperLabel;
  QLineEdit *options_quantileUpperName;
  QLabel *options_cutoffLowerLabel;
  QLineEdit *options_cutoffLowerName;
  QLabel *options_cutoffUpperLabel;
  QLineEdit *options_cutoffUpperName;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QFrame *line_3;

  QHBoxLayout * horizontalLayout;

  QCheckBox *wholeImageThresholdCheckBoxBox;
  
  //QLabel *longRunningWarning;

  void setupUi(QDialog *fDeepMedicNormalizer)
  {

    if (fDeepMedicNormalizer->objectName().isEmpty())
      fDeepMedicNormalizer->setObjectName(QString::fromUtf8("fDeepMedicNormalizer"));
    //fDeepMedicNormalizer->setWindowModality(Qt::NonModal);
    fDeepMedicNormalizer->resize(256, 256); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fDeepMedicNormalizer->sizePolicy().hasHeightForWidth());
    fDeepMedicNormalizer->setSizePolicy(sizePolicy);
    fDeepMedicNormalizer->setMinimumSize(QSize(0, 0));

    //fDeepMedicNormalizer->setModal(true);
    gridLayout = new QGridLayout(fDeepMedicNormalizer);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // input
    inputGroupBox = new QGroupBox(fDeepMedicNormalizer);
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

    // mask
    maskGroupBox = new QGroupBox(fDeepMedicNormalizer);
    maskGroupBox->setTitle(QString::fromStdString("Mask Image"));

    maskGridLayout = new QGridLayout(maskGroupBox);
    maskGridLayout->setObjectName(QString::fromUtf8("maskGridLayout"));

    maskImageLabel = new QLabel(maskGroupBox);
    sizePolicy.setHeightForWidth(maskImageLabel->sizePolicy().hasHeightForWidth());
    maskImageLabel->setSizePolicy(sizePolicy);

    maskImageName = new QLineEdit(" ");
    maskImageName->setObjectName(QString::fromUtf8("maskImageName"));
    sizePolicy.setHeightForWidth(maskImageName->sizePolicy().hasHeightForWidth());
    maskImageName->setSizePolicy(sizePolicy);
    maskImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    maskImageButton = new QPushButton(maskGroupBox);
    maskImageButton->setObjectName(QString::fromUtf8("maskImageButton"));
    maskImageButton->setText(QString("Browse"));
    //maskMaskButton->setToolTip(QString("Atlas based on which stripping is done"));

    maskGridLayout->addWidget(maskImageLabel, 0, 0, 1, 1);
    maskGridLayout->addWidget(maskImageName, 1, 0, 1, 1);
    maskGridLayout->addWidget(maskImageButton, 1, 1, 1, 1);

    // options
    optionsGroupBox = new QGroupBox(fDeepMedicNormalizer);
    optionsGroupBox->setTitle(QString::fromStdString("Options"));

    optionsGridLayout = new QGridLayout(optionsGroupBox);
    optionsGridLayout->setObjectName(QString::fromUtf8("optionsGridLayout"));
    int gridRowCounter = 0;

    options_quantileLowerLabel = new QLabel(optionsGroupBox);
    sizePolicy.setHeightForWidth(options_quantileLowerLabel->sizePolicy().hasHeightForWidth());
    options_quantileLowerLabel->setSizePolicy(sizePolicy);
    options_quantileLowerLabel->setText("Quantalization Lower");
    options_quantileLowerName = new QLineEdit(" ");
    options_quantileLowerName->setObjectName(QString::fromUtf8("options_quantileLowerName"));
    sizePolicy.setHeightForWidth(options_quantileLowerName->sizePolicy().hasHeightForWidth());
    options_quantileLowerName->setSizePolicy(sizePolicy);
    options_quantileLowerName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    options_quantileLowerName->setText("5");
    options_quantileLowerName->setValidator(new QIntValidator(0, 100, optionsGroupBox));

    optionsGridLayout->addWidget(options_quantileLowerLabel, gridRowCounter, 0);
    optionsGridLayout->addWidget(options_quantileLowerName, gridRowCounter, 1);

    options_quantileUpperLabel = new QLabel(optionsGroupBox);
    sizePolicy.setHeightForWidth(options_quantileUpperLabel->sizePolicy().hasHeightForWidth());
    options_quantileUpperLabel->setSizePolicy(sizePolicy);
    options_quantileUpperLabel->setText("Quantalization Upper");
    options_quantileUpperName = new QLineEdit(" ");
    options_quantileUpperName->setObjectName(QString::fromUtf8("options_quantileUpperName"));
    sizePolicy.setHeightForWidth(options_quantileUpperName->sizePolicy().hasHeightForWidth());
    options_quantileUpperName->setSizePolicy(sizePolicy);
    options_quantileUpperName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    options_quantileUpperName->setText("95");
    options_quantileUpperName->setValidator(new QIntValidator(0, 100, optionsGroupBox));

    gridRowCounter++;
    optionsGridLayout->addWidget(options_quantileUpperLabel, gridRowCounter, 0);
    optionsGridLayout->addWidget(options_quantileUpperName, gridRowCounter, 1);

    options_cutoffLowerLabel = new QLabel(optionsGroupBox);
    sizePolicy.setHeightForWidth(options_cutoffLowerLabel->sizePolicy().hasHeightForWidth());
    options_cutoffLowerLabel->setSizePolicy(sizePolicy);
    options_cutoffLowerLabel->setText("Cut-Off Lower");
    options_cutoffLowerName = new QLineEdit(" ");
    options_cutoffLowerName->setObjectName(QString::fromUtf8("options_cutoffLowerName"));
    sizePolicy.setHeightForWidth(options_cutoffLowerName->sizePolicy().hasHeightForWidth());
    options_cutoffLowerName->setSizePolicy(sizePolicy);
    options_cutoffLowerName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    options_cutoffLowerName->setText("3");
    options_cutoffLowerName->setValidator(new QIntValidator(0, 100, optionsGroupBox));

    gridRowCounter++;
    optionsGridLayout->addWidget(options_cutoffLowerLabel, gridRowCounter, 0);
    optionsGridLayout->addWidget(options_cutoffLowerName, gridRowCounter, 1);

    options_cutoffUpperLabel = new QLabel(optionsGroupBox);
    sizePolicy.setHeightForWidth(options_cutoffUpperLabel->sizePolicy().hasHeightForWidth());
    options_cutoffUpperLabel->setSizePolicy(sizePolicy);
    options_cutoffUpperLabel->setText("Cut-Off Upper");
    options_cutoffUpperName = new QLineEdit(" ");
    options_cutoffUpperName->setObjectName(QString::fromUtf8("options_cutoffUpperName"));
    sizePolicy.setHeightForWidth(options_cutoffUpperName->sizePolicy().hasHeightForWidth());
    options_cutoffUpperName->setSizePolicy(sizePolicy);
    options_cutoffUpperName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    options_cutoffUpperName->setText("3");
    options_cutoffUpperName->setValidator(new QIntValidator(0, 100, optionsGroupBox));

    gridRowCounter++;
    optionsGridLayout->addWidget(options_cutoffUpperLabel, gridRowCounter, 0);
    optionsGridLayout->addWidget(options_cutoffUpperName, gridRowCounter, 1);

    wholeImageThresholdCheckBoxBox = new QCheckBox("Whole Image Threshold");
    wholeImageThresholdCheckBoxBox->setEnabled(true);

    gridRowCounter++;
    optionsGridLayout->addWidget(wholeImageThresholdCheckBoxBox, gridRowCounter, 0, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fDeepMedicNormalizer);
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
    gridLayout->addWidget(inputGroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(maskGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(optionsGroupBox, 2, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 3, 0, 1, 2);


    confirmButton = new QPushButton(fDeepMedicNormalizer);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fDeepMedicNormalizer);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout->addWidget(confirmButton, 4, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 4, 1, 1, 1);

    retranslateUi(fDeepMedicNormalizer);

    QMetaObject::connectSlotsByName(fDeepMedicNormalizer);
  } // setupUi

  void retranslateUi(QDialog *fDeepMedicNormalizer)
  {
    fDeepMedicNormalizer->setWindowTitle(QApplication::translate("fDeepMedicNormalizer", "Z-Scoring Normalizer", 0));
    confirmButton->setText(QApplication::translate("fDeepMedicNormalizer", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fDeepMedicNormalizer", "Cancel", 0));

  } // retranslateUi

};

namespace Ui {
  class fDeepMedicNormalizer : public ui_fDeepMedicNormalizer {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDeepMedicNormalizer_H
