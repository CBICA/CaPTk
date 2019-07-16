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


  QGroupBox	*outputGroupBox;
  QGridLayout *outputGridLayout;

  QGroupBox	*confirmGroupBox;
  QGridLayout *confirmGridLayout;


  QLineEdit	*outputImageName;
  QLabel		*outputImageLabel;
  QPushButton *outputImageButton;
  QPushButton *mSplitModelDirectoryButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QRadioButton *mLinearKernel;
  QRadioButton *mRBFKernel;
  QRadioButton *mCrossValidation;
  QRadioButton *mSplitTrainTest;
  QRadioButton *mSplitTrain;
  QRadioButton *mSplitTest;

  
  QHBoxLayout * horizontalLayout;

  QLabel *longRunningWarning;
  QFrame *classifierFrame;
  QFrame *configurationFrame;

  QGroupBox	*classifierGroupBox;
  QGroupBox	*configurationGroupBox;
  QGridLayout *classifierGridLayout;
  QGridLayout *configurationGridLayout;

  QLineEdit	*cvValue;
  QLineEdit	*ttValue;
  QLineEdit	*mSplitModelDirectory;

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

    inputMaskName = new QLineEdit("");
    inputMaskName->setObjectName(QString::fromUtf8("inputMaskName"));
    sizePolicy.setHeightForWidth(inputMaskName->sizePolicy().hasHeightForWidth());
    inputMaskName->setSizePolicy(sizePolicy);
    inputMaskName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputMaskButton = new QPushButton(inputGroupBox);
    inputMaskButton->setObjectName(QString::fromUtf8("inputMaskButton"));
    inputMaskButton->setText(QString("Browse"));

    //inputGridLayout->addWidget(classifierGroup(), 0, 0);
    //inputGridLayout->addWidget(configurationGroup(), 1, 0);

    // output
    classifierGroupBox = new QGroupBox(fTrainingSimulator);
    classifierGroupBox->setTitle(QString::fromStdString("SVM Classification kernel"));
    classifierGridLayout = new QGridLayout(classifierGroupBox);
    classifierGridLayout->setObjectName(QString::fromUtf8("classifierGridLayout"));
    mLinearKernel = new QRadioButton("SVM: Linear");
    mLinearKernel->setEnabled(true);
    mRBFKernel = new QRadioButton("SVM: RBF");
    mRBFKernel->setEnabled(true);
    classifierGridLayout->addWidget(mLinearKernel, 0, 0, 1, 1);
    classifierGridLayout->addWidget(mRBFKernel, 0, 1, 1, 1);


    configurationGroupBox = new QGroupBox(fTrainingSimulator);
    configurationGroupBox->setTitle(QString::fromStdString("Configuration"));
    configurationGridLayout = new QGridLayout(configurationGroupBox);
    configurationGridLayout->setObjectName(QString::fromUtf8("configurationGridLayout"));
    mCrossValidation = new QRadioButton("CrossValidation");
    mCrossValidation->setEnabled(true);
    mSplitTrainTest = new QRadioButton("Split TrainTest");
    mSplitTrainTest->setEnabled(true);
    mSplitTrain = new QRadioButton("Split Train");
    mSplitTrain->setEnabled(true);
    mSplitTest = new QRadioButton("Split Test");
    mSplitTest->setEnabled(true);


    cvValue = new QLineEdit("");
    cvValue->setObjectName(QString::fromUtf8("cvValue"));
    sizePolicy.setHeightForWidth(cvValue->sizePolicy().hasHeightForWidth());
    cvValue->setSizePolicy(sizePolicy);
    cvValue->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    ttValue = new QLineEdit("");
    ttValue->setObjectName(QString::fromUtf8("ttValue"));
    sizePolicy.setHeightForWidth(ttValue->sizePolicy().hasHeightForWidth());
    ttValue->setSizePolicy(sizePolicy);
    ttValue->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    mSplitModelDirectory = new QLineEdit("");
    mSplitModelDirectory->setObjectName(QString::fromUtf8("mSplitModelDirectory"));
    sizePolicy.setHeightForWidth(mSplitModelDirectory->sizePolicy().hasHeightForWidth());
    mSplitModelDirectory->setSizePolicy(sizePolicy);
    mSplitModelDirectory->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    mSplitModelDirectoryButton = new QPushButton(configurationGroupBox);
    mSplitModelDirectoryButton->setObjectName(QString::fromUtf8("mSplitModelDirectoryButton"));
    mSplitModelDirectoryButton->setText(QString("Browse"));



    configurationGridLayout->addWidget(mCrossValidation, 0, 0, 1, 1);
    configurationGridLayout->addWidget(cvValue, 0, 1, 1, 1);
    configurationGridLayout->addWidget(mSplitTrainTest, 0, 2, 1, 1);
    configurationGridLayout->addWidget(ttValue, 0, 3, 1, 1);
    configurationGridLayout->addWidget(mSplitTrain, 0, 4, 1, 1);
    configurationGridLayout->addWidget(mSplitTest, 0, 5, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectory, 0, 6, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectoryButton, 0, 7, 1, 1);



    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 0, 1, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 0, 2, 1, 1);

    inputGridLayout->addWidget(inputMaskLabel, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputMaskName, 1, 1, 1, 1);
    inputGridLayout->addWidget(inputMaskButton, 1, 2, 1, 1);


    // output
    outputGroupBox = new QGroupBox(fTrainingSimulator);
    outputGroupBox->setTitle(QString::fromStdString("Output"));
    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputImageLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    //outputImageLabel->setSizePolicy(sizePolicy);

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

    //outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);
    outputGridLayout->addWidget(longRunningWarning, 1, 0, 1, 2);

    confirmGroupBox = new QGroupBox(fTrainingSimulator);
    confirmGroupBox->setTitle(QString::fromStdString(""));
    confirmGridLayout = new QGridLayout(confirmGroupBox);
    confirmGridLayout->setObjectName(QString::fromUtf8("confirmGridLayout"));

    confirmButton = new QPushButton(confirmGroupBox);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    cancelButton = new QPushButton(confirmGroupBox);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    
    confirmGridLayout->addWidget(confirmButton, 0, 0, 1, 1);
    confirmGridLayout->addWidget(cancelButton, 0, 1, 1, 1);


    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(classifierGroupBox, 2, 0, 1, 2);
    gridLayout->addWidget(configurationGroupBox, 3, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 4, 0, 1, 2);
    gridLayout->addWidget(confirmGroupBox, 5, 0, 1, 2);

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
  } // retranslateUi
};

namespace Ui {
  class fTrainingSimulator : public ui_fTrainingSimulator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fTrainingSimulator_H
