#ifndef ui_fPerfusionAligner_H
#define ui_fPerfusionAligner_H

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

class ui_fPerfusionAligner
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;


  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;
  
  QLabel *inputImageLabel;
  QLineEdit *inputImageName;
  QPushButton *inputImageButton;

  QLabel *inputMaskLabel;
  QLineEdit *inputMaskName;
  QPushButton *inputMaskButton;

  QLabel *inputAfterPointsLabel;
  QLineEdit *inputAfterPointsName;

  QLabel *inputBeforePointsLabel;
  QLineEdit *inputBeforePointsName;

  QLabel *inputEchoTimeLabel;
  QLineEdit *inputEchoTimeName;

  QLabel *inputStdDevLabel;
  QLineEdit *inputScaleDropBeforeMeanName;

  QLabel *inputBaselineLabel;
  QLineEdit *inputBaselineLine;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QLabel *longRunningWarning;
  QPushButton *outputImageButton;

  QLabel *outputEchoTimeLabel;
  QLineEdit *outputEchoTimeName;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fPerfusionAlignmentDialog)
  {
    if (fPerfusionAlignmentDialog->objectName().isEmpty())
      fPerfusionAlignmentDialog->setObjectName(QString::fromUtf8("fPerfusionAlignmentDialog"));
    //fPerfusionAlignmentDialog->setWindowModality(Qt::NonModal);
    fPerfusionAlignmentDialog->resize(400, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPerfusionAlignmentDialog->sizePolicy().hasHeightForWidth());
    fPerfusionAlignmentDialog->setSizePolicy(sizePolicy);
    fPerfusionAlignmentDialog->setMinimumSize(QSize(0, 0));

    //fPerfusionAlignmentDialog->setModal(true);
    gridLayout = new QGridLayout(fPerfusionAlignmentDialog);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    //input
    inputGroupBox = new QGroupBox(fPerfusionAlignmentDialog);
    inputGroupBox->setTitle(QString::fromStdString("Input Data"));

    inputGridLayout = new QGridLayout(inputGroupBox);
    inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputImageLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputImageLabel->sizePolicy().hasHeightForWidth());
    inputImageLabel->setSizePolicy(sizePolicy);

    inputImageName = new QLineEdit("");
    inputImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputImageName->sizePolicy().hasHeightForWidth());
    inputImageName->setSizePolicy(sizePolicy);
    inputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputImageButton = new QPushButton(inputGroupBox);
    inputImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputImageButton->setText("...");
    inputImageButton->setToolTip("Browse and select file");

    inputMaskLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputMaskLabel->sizePolicy().hasHeightForWidth());
    inputMaskLabel->setSizePolicy(sizePolicy);

    inputMaskName = new QLineEdit("Select Mask File");
    inputMaskName->setObjectName(QString::fromUtf8("inputMaskName"));
    sizePolicy.setHeightForWidth(inputMaskName->sizePolicy().hasHeightForWidth());
    inputMaskName->setSizePolicy(sizePolicy);
    inputMaskName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputMaskButton = new QPushButton(inputGroupBox);
    inputMaskButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputMaskButton->setText("...");
    inputMaskButton->setToolTip("Browse and select Mask file");

    inputBeforePointsLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputBeforePointsLabel->sizePolicy().hasHeightForWidth());
    inputBeforePointsLabel->setSizePolicy(sizePolicy);

    inputAfterPointsLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputAfterPointsLabel->sizePolicy().hasHeightForWidth());
    inputAfterPointsLabel->setSizePolicy(sizePolicy);

    inputEchoTimeLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputEchoTimeLabel->sizePolicy().hasHeightForWidth());
    inputEchoTimeLabel->setSizePolicy(sizePolicy);
    
    inputBeforePointsName = new QLineEdit("5");
    inputBeforePointsName->setObjectName(QString::fromUtf8("inputBeforePointsName"));
    sizePolicy.setHeightForWidth(inputBeforePointsName->sizePolicy().hasHeightForWidth());
    inputBeforePointsName->setSizePolicy(sizePolicy);
    inputBeforePointsName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    inputBeforePointsName->setValidator(new QIntValidator(0, 100, inputGroupBox));

    inputAfterPointsName = new QLineEdit("15");
    inputAfterPointsName->setObjectName(QString::fromUtf8("inputAfterPointsName"));
    sizePolicy.setHeightForWidth(inputAfterPointsName->sizePolicy().hasHeightForWidth());
    inputAfterPointsName->setSizePolicy(sizePolicy);
    inputAfterPointsName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    inputAfterPointsName->setValidator(new QIntValidator(0, 100, inputGroupBox));


    inputEchoTimeName = new QLineEdit("2");
    inputEchoTimeName->setObjectName(QString::fromUtf8("inputEchoTimeName"));
    sizePolicy.setHeightForWidth(inputEchoTimeName->sizePolicy().hasHeightForWidth());
    inputEchoTimeName->setSizePolicy(sizePolicy);
    inputEchoTimeName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    inputEchoTimeName->setValidator(new QIntValidator(0, 100, inputGroupBox));

    inputStdDevLabel = new QLabel(inputGroupBox);
    inputStdDevLabel->setText("Intensity drop in mean curve inside ROI");
    inputScaleDropBeforeMeanName = new QLineEdit("100");
    inputScaleDropBeforeMeanName->setObjectName(QString::fromUtf8("inputEchoTimeName"));
    sizePolicy.setHeightForWidth(inputScaleDropBeforeMeanName->sizePolicy().hasHeightForWidth());
    inputScaleDropBeforeMeanName->setSizePolicy(sizePolicy);
    inputScaleDropBeforeMeanName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    inputScaleDropBeforeMeanName->setValidator(new QIntValidator(0, 10000, inputGroupBox));

    inputBaselineLabel = new QLabel(inputGroupBox);
    inputBaselineLabel->setText("Maximum intensity before drop in mean curve inside ROI");
    inputBaselineLine = new QLineEdit("300");
    inputBaselineLine->setObjectName(QString::fromUtf8("inputEchoTimeName"));
    sizePolicy.setHeightForWidth(inputBaselineLine->sizePolicy().hasHeightForWidth());
    inputBaselineLine->setSizePolicy(sizePolicy);
    inputBaselineLine->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    inputBaselineLine->setValidator(new QIntValidator(0, 20000, inputGroupBox));

    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 0, 1, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 0, 2, 1, 1);

    inputGridLayout->addWidget(inputMaskLabel, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputMaskName, 1, 1, 1, 1);
    inputGridLayout->addWidget(inputMaskButton, 1, 2, 1, 1);

    inputGridLayout->addWidget(inputBeforePointsLabel, 2, 0, 1, 1);
    inputGridLayout->addWidget(inputBeforePointsName, 2, 1, 1, 1);

    inputGridLayout->addWidget(inputAfterPointsLabel, 3, 0, 1, 1);
    inputGridLayout->addWidget(inputAfterPointsName, 3, 1, 1, 1);

    inputGridLayout->addWidget(inputEchoTimeLabel, 4, 0, 1, 1);
    inputGridLayout->addWidget(inputEchoTimeName, 4, 1, 1, 1);
    inputGridLayout->addWidget(inputStdDevLabel, 5, 0, 1, 1);
    inputGridLayout->addWidget(inputScaleDropBeforeMeanName, 5, 1, 1, 1);
    inputGridLayout->addWidget(inputBaselineLabel, 6, 0, 1, 1);
    inputGridLayout->addWidget(inputBaselineLine, 6, 1, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fPerfusionAlignmentDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output"));
    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputEchoTimeLabel = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(outputEchoTimeLabel->sizePolicy().hasHeightForWidth());
    outputEchoTimeLabel->setSizePolicy(sizePolicy);
    outputEchoTimeLabel->setText("Desired Time-resolution (seconds)");

    outputEchoTimeName = new QLineEdit("1");
    outputEchoTimeName->setObjectName(QString::fromUtf8("outputEchoTimeName"));
    sizePolicy.setHeightForWidth(outputEchoTimeName->sizePolicy().hasHeightForWidth());
    outputEchoTimeName->setSizePolicy(sizePolicy);
    outputEchoTimeName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    outputEchoTimeName->setValidator(new QIntValidator(0, 100, outputGroupBox));

    outputImageName = new QLineEdit("");
    outputImageName->setObjectName(QString::fromUtf8("outputImageName"));
    sizePolicy.setHeightForWidth(outputImageName->sizePolicy().hasHeightForWidth());
    outputImageName->setSizePolicy(sizePolicy);
    outputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputImageButton = new QPushButton(outputGroupBox);
    outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
    outputImageButton->setText(QString("Browse"));

    longRunningWarning = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    longRunningWarning->setSizePolicy(sizePolicy);
    longRunningWarning->setAlignment(Qt::AlignRight);
    longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

    outputGridLayout->addWidget(outputEchoTimeLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputEchoTimeName, 0, 1, 1, 1);
    outputGridLayout->addWidget(outputImageName, 1, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 1, 1, 1, 1);
    outputGridLayout->addWidget(longRunningWarning, 2, 0, 1, 2);

    // put the layout in perspective
    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fPerfusionAlignmentDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fPerfusionAlignmentDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fPerfusionAlignmentDialog);

    QMetaObject::connectSlotsByName(fPerfusionAlignmentDialog);
  } // setupUi

  void retranslateUi(QDialog *fPerfusionAlignmentDialog)
  {
    // fPerfusionAlignmentDialog->setWindowTitle(QApplication::translate("fPerfusionAlignmentDialog", "Perfusion Derivatives", 0, QApplication::UnicodeUTF8));
    // confirmButton->setText(QApplication::translate("fPerfusionAlignmentDialog", "Confirm", 0, QApplication::UnicodeUTF8));
    // cancelButton->setText(QApplication::translate("fPerfusionAlignmentDialog", "Cancel", 0, QApplication::UnicodeUTF8));
    // m_rcbv->setText(QApplication::translate("fPerfusionAlignmentDialog", "Relative Cerebral Blood Volume", 0, QApplication::UnicodeUTF8));
    // m_psr->setText(QApplication::translate("fPerfusionAlignmentDialog", "Percent Signal Recovery", 0, QApplication::UnicodeUTF8));
    // m_ph->setText(QApplication::translate("fPerfusionAlignmentDialog", "Peak Height", 0, QApplication::UnicodeUTF8));
    // NEW CHANGES
    fPerfusionAlignmentDialog->setWindowTitle(QApplication::translate("fPerfusionAlignmentDialog", "Perfusion Derivatives", 0));
    inputBeforePointsLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "# of Points before drop", 0));
    inputAfterPointsLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "# of Points after drop", 0));
    inputEchoTimeLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "Time-resolution (seconds)", 0));
    inputImageLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "DSC-MRI Image", 0));

    confirmButton->setText(QApplication::translate("fPerfusionAlignmentDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fPerfusionAlignmentDialog", "Cancel", 0));
  } // retranslateUi
};

namespace Ui {
  class fPerfusionAligner : public ui_fPerfusionAligner {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPerfusionAligner_H