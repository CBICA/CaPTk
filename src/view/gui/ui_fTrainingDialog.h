#ifndef ui_fTrainingSimulator_H
#define ui_fTrainingSimulator_H

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

class ui_fTrainingSimulator
{
public:
  QGridLayout * gridLayout;
  QVBoxLayout * verticalLayout;


  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;

  QLabel		*inputImageLabel;
  QLineEdit	*inputImageName;
  QPushButton *inputImageButton;

  QLabel		*inputMaskLabel;
  QLineEdit	*inputMaskName;
  QPushButton *inputMaskButton;


  QLabel		*inputFoldsLabel;
  QLineEdit	*inputFoldsName;
  
  
  QGroupBox	*outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit	*outputImageName;
  QLabel		*outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QRadioButton *mLinearKernel;
  QRadioButton *mRBFKernel;

  
  QHBoxLayout * horizontalLayout;

  QLabel *longRunningWarning;

  void setupUi(QDialog *fTrainingSimulator)
  {
    if (fTrainingSimulator->objectName().isEmpty())
      fTrainingSimulator->setObjectName(QString::fromUtf8("fTrainingSimulator"));

    fTrainingSimulator->resize(400, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fTrainingSimulator->sizePolicy().hasHeightForWidth());
    fTrainingSimulator->setSizePolicy(sizePolicy);
    fTrainingSimulator->setMinimumSize(QSize(0, 0));

    fTrainingSimulator->setModal(true);
    gridLayout = new QGridLayout(fTrainingSimulator);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    //input
    inputGroupBox = new QGroupBox(fTrainingSimulator);
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


    inputMaskLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputMaskLabel->sizePolicy().hasHeightForWidth());
    inputMaskLabel->setSizePolicy(sizePolicy);

    inputFoldsLabel = new QLabel(inputGroupBox);
    inputFoldsLabel->setSizePolicy(sizePolicy);

    inputMaskName = new QLineEdit("");
    inputMaskName->setObjectName(QString::fromUtf8("inputMaskName"));
    sizePolicy.setHeightForWidth(inputMaskName->sizePolicy().hasHeightForWidth());
    inputMaskName->setSizePolicy(sizePolicy);
    inputMaskName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputFoldsName = new QLineEdit("");
    inputFoldsName->setObjectName(QString::fromUtf8("inputFoldsName"));
    sizePolicy.setHeightForWidth(inputFoldsName->sizePolicy().hasHeightForWidth());
    inputFoldsName->setSizePolicy(sizePolicy);
    inputFoldsName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    inputMaskButton = new QPushButton(inputGroupBox);
    inputMaskButton->setObjectName(QString::fromUtf8("inputMaskButton"));
    inputMaskButton->setText(QString("Browse"));

    mLinearKernel = new QRadioButton("SVM: Linear");
    mLinearKernel->setEnabled(true);
    mRBFKernel = new QRadioButton("SVM: RBF");
    mRBFKernel->setEnabled(true);



    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 0, 1, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 0, 2, 1, 1);

    inputGridLayout->addWidget(inputMaskLabel, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputMaskName, 1, 1, 1, 1);
    inputGridLayout->addWidget(inputMaskButton, 1, 2, 1, 1);


    inputGridLayout->addWidget(inputFoldsLabel, 2, 0, 1, 1);
    inputGridLayout->addWidget(inputFoldsName, 2, 1, 1, 1);

    inputGridLayout->addWidget(mLinearKernel, 3, 0, 1, 1);
    inputGridLayout->addWidget(mRBFKernel, 3, 1, 1, 1);


    // output
    outputGroupBox = new QGroupBox(fTrainingSimulator);
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
    outputImageButton->setText(QString("Browse"));

    longRunningWarning = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    longRunningWarning->setSizePolicy(sizePolicy);
    longRunningWarning->setAlignment(Qt::AlignRight);
    longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

    //outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);
    outputGridLayout->addWidget(longRunningWarning, 1, 0, 1, 2);

    // put the layout in perspective
    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fTrainingSimulator);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fTrainingSimulator);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fTrainingSimulator);

    QMetaObject::connectSlotsByName(fTrainingSimulator);
  } // setupUi

  void retranslateUi(QDialog *fTrainingSimulator)
  {
    // NEW CHANGES
    fTrainingSimulator->setWindowTitle(QApplication::translate("fTrainingSimulator", "Training Module", 0));
    confirmButton->setText(QApplication::translate("fTrainingSimulator", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fTrainingSimulator", "Cancel", 0));
    inputImageLabel->setText(QApplication::translate("fTrainingSimulator", "Features File:", 0));
    inputMaskLabel->setText(QApplication::translate("fTrainingSimulator", "Target File:", 0));
    inputFoldsLabel->setText(QApplication::translate("fTrainingSimulator", "No. of folds:", 0));
  } // retranslateUi
};

namespace Ui {
  class fTrainingSimulator : public ui_fTrainingSimulator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fTrainingSimulator_H
