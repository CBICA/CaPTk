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

  QLabel *inputT1ceImageLabel;
  QLineEdit *inputT1ceImageName;
  QPushButton *inputT1ceImageButton;

  QLabel *inputAfterPointsLabel;
  QLineEdit *inputAfterPointsName;

  QLabel *inputBeforePointsLabel;
  QLineEdit *inputBeforePointsName;

  QLabel *inputEchoTimeLabel;
  QLineEdit *inputEchoTimeName;


  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QLabel *longRunningWarning;
  QPushButton *outputImageButton;


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
    inputImageButton->setText(QString("Browse"));






    inputT1ceImageLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputT1ceImageLabel->sizePolicy().hasHeightForWidth());
    inputT1ceImageLabel->setSizePolicy(sizePolicy);

    inputT1ceImageName = new QLineEdit("");
    inputT1ceImageName->setObjectName(QString::fromUtf8("inputT1ceImageName"));
    sizePolicy.setHeightForWidth(inputT1ceImageName->sizePolicy().hasHeightForWidth());
    inputT1ceImageName->setSizePolicy(sizePolicy);
    inputT1ceImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputT1ceImageButton = new QPushButton(inputGroupBox);
    inputT1ceImageButton->setObjectName(QString::fromUtf8("inputT1ceImageButton"));
    inputT1ceImageButton->setText(QString("Browse"));


    inputBeforePointsLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputBeforePointsLabel->sizePolicy().hasHeightForWidth());
    inputBeforePointsLabel->setSizePolicy(sizePolicy);

    inputAfterPointsLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputAfterPointsLabel->sizePolicy().hasHeightForWidth());
    inputAfterPointsLabel->setSizePolicy(sizePolicy);

    inputEchoTimeLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputEchoTimeLabel->sizePolicy().hasHeightForWidth());
    inputEchoTimeLabel->setSizePolicy(sizePolicy);



    inputBeforePointsName = new QLineEdit("");
    inputBeforePointsName->setObjectName(QString::fromUtf8("inputBeforePointsName"));
    sizePolicy.setHeightForWidth(inputBeforePointsName->sizePolicy().hasHeightForWidth());
    inputBeforePointsName->setSizePolicy(sizePolicy);
    inputBeforePointsName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputAfterPointsName = new QLineEdit("");
    inputAfterPointsName->setObjectName(QString::fromUtf8("inputAfterPointsName"));
    sizePolicy.setHeightForWidth(inputAfterPointsName->sizePolicy().hasHeightForWidth());
    inputAfterPointsName->setSizePolicy(sizePolicy);
    inputAfterPointsName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    inputEchoTimeName = new QLineEdit("");
    inputEchoTimeName->setObjectName(QString::fromUtf8("inputEchoTimeName"));
    sizePolicy.setHeightForWidth(inputEchoTimeName->sizePolicy().hasHeightForWidth());
    inputEchoTimeName->setSizePolicy(sizePolicy);
    inputEchoTimeName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);





    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 0, 1, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 0, 2, 1, 1);

    inputGridLayout->addWidget(inputT1ceImageLabel, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputT1ceImageName, 1, 1, 1, 1);
    inputGridLayout->addWidget(inputT1ceImageButton, 1, 2, 1, 1);

    inputGridLayout->addWidget(inputBeforePointsLabel, 2, 0, 1, 1);
    inputGridLayout->addWidget(inputBeforePointsName, 2, 1, 1, 1);

    inputGridLayout->addWidget(inputAfterPointsLabel, 3, 0, 1, 1);
    inputGridLayout->addWidget(inputAfterPointsName, 3, 1, 1, 1);

    inputGridLayout->addWidget(inputEchoTimeLabel, 4, 0, 1, 1);
    inputGridLayout->addWidget(inputEchoTimeName, 4, 1, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fPerfusionAlignmentDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output"));
    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

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

    outputGridLayout->addWidget(outputImageName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);
    outputGridLayout->addWidget(longRunningWarning, 1, 0, 1, 2);

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
    inputEchoTimeLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "Echo time (seconds)", 0));
    inputImageLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "DSC-MRI Image", 0));
    inputT1ceImageLabel->setText(QApplication::translate("fPerfusionAlignmentDialog", "T1 PostContrast Image", 0));

    confirmButton->setText(QApplication::translate("fPerfusionAlignmentDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fPerfusionAlignmentDialog", "Cancel", 0));
  } // retranslateUi
};

namespace Ui {
  class fPerfusionAligner : public ui_fPerfusionAligner {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPerfusionAligner_H