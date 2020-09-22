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

  QLabel		*inputFeaturesLabel;
  QLineEdit	*inputFeaturesName;
  QPushButton *inputFeaturesButton;

  QLabel		*inputTargetLabel;
  QLineEdit	*inputTargetName;
  QPushButton *inputTargetButton;


  QGroupBox	*outputGroupBox;
  QGridLayout *outputGridLayout;

  QGroupBox	*confirmGroupBox;
  QGridLayout *confirmGridLayout;


  QLineEdit	*outputDirectoryName;
  QLabel		*outputDirectoryLabel;
  QPushButton *outputDirectoryButton;
  QPushButton *mSplitModelDirectoryButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QRadioButton *mLinearKernel;
  QRadioButton *mRBFKernel;
  QRadioButton *mSVMFFS;
  QRadioButton *mESFS;
  QRadioButton *mResubstitution;
  QRadioButton *mFiveFold;

  QRadioButton *mCrossValidation;
  QRadioButton *mSplitTrain;
  QRadioButton *mSplitTest;

  
  QHBoxLayout * horizontalLayout;

  QLabel *longRunningWarning;
  QFrame *classifierFrame;
  QFrame *configurationFrame;

  QGroupBox	*classifierGroupBox;
  QGroupBox	*paramsGroupBox;
  QGroupBox	*fsGroupBox;

  QGroupBox	*configurationGroupBox;
  QGridLayout *classifierGridLayout;
  QGridLayout *paramsGridLayout;
  QGridLayout *fsGridLayout;
  QGridLayout *configurationGridLayout;

  QLabel		*cvLabel;
  QLineEdit	*cvValue;

  QLabel    *mSplitModelDirectoryLabel;
  QLineEdit	*mSplitModelDirectory;

  QLabel *parametersLabel;
  QCheckBox *mOptimization;
  QLabel *crossvalidationLabel;

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

    inputFeaturesLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputFeaturesLabel->sizePolicy().hasHeightForWidth());
    inputFeaturesLabel->setSizePolicy(sizePolicy);

    inputFeaturesName = new QLineEdit("");
    inputFeaturesName->setObjectName(QString::fromUtf8("inputFeaturesName"));
    sizePolicy.setHeightForWidth(inputFeaturesName->sizePolicy().hasHeightForWidth());
    inputFeaturesName->setSizePolicy(sizePolicy);
    inputFeaturesName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputFeaturesButton = new QPushButton(inputGroupBox);
    inputFeaturesButton->setObjectName(QString::fromUtf8("inputFeaturesButton"));
    inputFeaturesButton->setText(QString("Browse"));


    inputTargetLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputTargetLabel->sizePolicy().hasHeightForWidth());
    inputTargetLabel->setSizePolicy(sizePolicy);

    inputTargetName = new QLineEdit("");
    inputTargetName->setObjectName(QString::fromUtf8("inputTargetName"));
    sizePolicy.setHeightForWidth(inputTargetName->sizePolicy().hasHeightForWidth());
    inputTargetName->setSizePolicy(sizePolicy);
    inputTargetName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputTargetButton = new QPushButton(inputGroupBox);
    inputTargetButton->setObjectName(QString::fromUtf8("inputTargetButton"));
    inputTargetButton->setText(QString("Browse"));

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

    paramsGroupBox = new QGroupBox(fTrainingSimulator);
    paramsGroupBox->setTitle(QString::fromStdString("Parameter Optimization/Cross-validation"));
    paramsGridLayout = new QGridLayout(paramsGroupBox);
    paramsGridLayout->setObjectName(QString::fromUtf8("paramsGridLayout"));

    parametersLabel = new QLabel(paramsGroupBox);
    sizePolicy.setHeightForWidth(parametersLabel->sizePolicy().hasHeightForWidth());
    parametersLabel->setSizePolicy(sizePolicy);
    mOptimization = new QCheckBox("Yes/No");
    mOptimization->setEnabled(true);


    crossvalidationLabel = new QLabel(paramsGroupBox);
    sizePolicy.setHeightForWidth(crossvalidationLabel->sizePolicy().hasHeightForWidth());
    crossvalidationLabel->setSizePolicy(sizePolicy);
    mResubstitution = new QRadioButton("Resubstitution");
    mResubstitution->setEnabled(true);
    mFiveFold = new QRadioButton("5-fold");
    mFiveFold->setEnabled(true);
    paramsGridLayout->addWidget(parametersLabel, 0, 0, 1, 1); 
    paramsGridLayout->addWidget(crossvalidationLabel, 0, 2, 1, 1); 
    paramsGridLayout->addWidget(mOptimization, 1, 0, 1, 1);
    paramsGridLayout->addWidget(mResubstitution, 1, 2, 1, 1);
    paramsGridLayout->addWidget(mFiveFold, 1,3, 1, 1);

    fsGroupBox = new QGroupBox(fTrainingSimulator);
    fsGroupBox->setTitle(QString::fromStdString("Feature selection"));
    fsGridLayout = new QGridLayout(fsGroupBox);
    fsGridLayout->setObjectName(QString::fromUtf8("fsGridLayout"));
    mSVMFFS = new QRadioButton("SVM: Forward feature selection");
    mSVMFFS->setEnabled(true);
    mESFS = new QRadioButton("Effect Size feature selection");
    mESFS->setEnabled(true);
    fsGridLayout->addWidget(mSVMFFS, 0, 0, 1, 1);
    fsGridLayout->addWidget(mESFS, 0, 1, 1, 1);

    configurationGroupBox = new QGroupBox(fTrainingSimulator);
    configurationGroupBox->setTitle(QString::fromStdString("Configuration"));
    configurationGridLayout = new QGridLayout(configurationGroupBox);
    configurationGridLayout->setObjectName(QString::fromUtf8("configurationGridLayout"));
    cvLabel = new QLabel(configurationGroupBox);
    sizePolicy.setHeightForWidth(cvLabel->sizePolicy().hasHeightForWidth());
    cvLabel->setSizePolicy(sizePolicy);
    mSplitModelDirectoryLabel = new QLabel(configurationGroupBox);
    sizePolicy.setHeightForWidth(mSplitModelDirectoryLabel->sizePolicy().hasHeightForWidth());
    mSplitModelDirectoryLabel->setSizePolicy(sizePolicy);

    mCrossValidation = new QRadioButton("CrossValidation");
    mCrossValidation->setEnabled(true);
    mSplitTrain = new QRadioButton("Split Train");
    mSplitTrain->setEnabled(true);
    mSplitTest = new QRadioButton("Split Test");
    mSplitTest->setEnabled(true);

    cvValue = new QLineEdit("");
    cvValue->setObjectName(QString::fromUtf8("cvValue"));
    sizePolicy.setHeightForWidth(cvValue->sizePolicy().hasHeightForWidth());
    cvValue->setSizePolicy(sizePolicy);
    cvValue->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    mSplitModelDirectory = new QLineEdit("");
    mSplitModelDirectory->setObjectName(QString::fromUtf8("mSplitModelDirectory"));
    sizePolicy.setHeightForWidth(mSplitModelDirectory->sizePolicy().hasHeightForWidth());
    mSplitModelDirectory->setSizePolicy(sizePolicy);
    mSplitModelDirectory->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    mSplitModelDirectoryButton = new QPushButton(configurationGroupBox);
    mSplitModelDirectoryButton->setObjectName(QString::fromUtf8("mSplitModelDirectoryButton"));
    mSplitModelDirectoryButton->setText(QString("Browse"));

    configurationGridLayout->addWidget(mCrossValidation, 0, 0, 1, 1);
    configurationGridLayout->addWidget(cvLabel, 0, 1, 1, 1);
    configurationGridLayout->addWidget(cvValue, 0, 2, 1, 1);

    configurationGridLayout->addWidget(mSplitTrain, 2,0, 1, 1);
    configurationGridLayout->addWidget(mSplitTest, 3,0, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectoryLabel, 3, 1, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectory, 3,2, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectoryButton, 3,3, 1, 1);



    inputGridLayout->addWidget(inputFeaturesLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputFeaturesName, 0, 1, 1, 1);
    inputGridLayout->addWidget(inputFeaturesButton, 0, 2, 1, 1);

    inputGridLayout->addWidget(inputTargetLabel, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputTargetName, 1, 1, 1, 1);
    inputGridLayout->addWidget(inputTargetButton, 1, 2, 1, 1);


    // output
    outputGroupBox = new QGroupBox(fTrainingSimulator);
    outputGroupBox->setTitle(QString::fromStdString("Output"));
    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputDirectoryLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    //outputDirectoryLabel->setSizePolicy(sizePolicy);

    outputDirectoryName = new QLineEdit("");
    outputDirectoryName->setObjectName(QString::fromUtf8("outputDirectoryName"));
    sizePolicy.setHeightForWidth(outputDirectoryName->sizePolicy().hasHeightForWidth());
    outputDirectoryName->setSizePolicy(sizePolicy);
    outputDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputDirectoryButton = new QPushButton(outputGroupBox);
    outputDirectoryButton->setObjectName(QString::fromUtf8("outputDirectoryButton"));
    outputDirectoryButton->setText(QString("Browse"));

    longRunningWarning = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    longRunningWarning->setSizePolicy(sizePolicy);
    longRunningWarning->setAlignment(Qt::AlignRight);
    longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

    //outputGridLayout->addWidget(outputDirectoryLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirectoryName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirectoryButton, 0, 1, 1, 1);
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
    gridLayout->addWidget(paramsGroupBox, 3, 0, 1, 2);
    gridLayout->addWidget(fsGroupBox, 4, 0, 1, 2);
    gridLayout->addWidget(configurationGroupBox, 5, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 6, 0, 1, 2);
    gridLayout->addWidget(confirmGroupBox, 7, 0, 1, 2);

    retranslateUi(fTrainingSimulator);

    QMetaObject::connectSlotsByName(fTrainingSimulator);
  } // setupUi

  void retranslateUi(QDialog *fTrainingSimulator)
  {
    // NEW CHANGES
    fTrainingSimulator->setWindowTitle(QApplication::translate("fTrainingSimulator", "Training Module", 0));
    confirmButton->setText(QApplication::translate("fTrainingSimulator", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fTrainingSimulator", "Cancel", 0));
    inputFeaturesLabel->setText(QApplication::translate("fTrainingSimulator", "Features File:", 0));
    inputTargetLabel->setText(QApplication::translate("fTrainingSimulator", "Target File:", 0));
    cvLabel->setText(QApplication::translate("fTrainingSimulator", "No. of folds:", 0));
    mSplitModelDirectoryLabel->setText(QApplication::translate("fTrainingSimulator", "Model directory:", 0));
    parametersLabel->setText(QApplication::translate("fTrainingSimulator", "Optimize classifier's parameters", 0));
    crossvalidationLabel->setText(QApplication::translate("fTrainingSimulator", "Cross-validation", 0));

  } // retranslateUi
};

namespace Ui {
  class fTrainingSimulator : public ui_fTrainingSimulator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fTrainingSimulator_H
