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
#include <QtWidgets/QProgressBar>

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
  QRadioButton* mPolynomialKernel;
  QRadioButton* mIntersectionKernel;
  QRadioButton* mSigmoidKernel;
  QRadioButton* mChiSquaredKernel;
  QRadioButton* mRandomForest;
  QRadioButton* mSGDSVM;
  QRadioButton* mBoostedTrees;
  QRadioButton* mSelectOptimalML;

  QRadioButton *mSVMFFS;
  QRadioButton *mESFS;
  QRadioButton *mBackwardsFS;
  QRadioButton *mRandomForestFS;
  QRadioButton *mReliefFFS;
  QRadioButton* mSelectOptimalFS;

  QRadioButton *mResubstitution;
  QRadioButton *mFiveFold;
  QLabel *mNumberOfFeaturesLabel;

  QGroupBox* forwardfsTerminationGroupBox;
  QGridLayout* forwardfsTerminationGridLayout;
  QSpinBox *mNumberOfFeaturesSpinbox;
  QRadioButton* mUseFeatureCount;
  QRadioButton* mUseConvergenceThreshold;

  QLabel* validatePredictionsLabel;
  QCheckBox *mTestPredictionsAgainstProvidedLabels;

  QRadioButton *mCrossValidation;
  QRadioButton *mSplitTrain;
  QRadioButton *mSplitTest;

  
  QHBoxLayout * horizontalLayout;

  QLabel *longRunningWarning;
  QFrame *classifierFrame;
  QFrame *configurationFrame;

  QGroupBox	*classifierGroupBox;
  QGroupBox	*paramsGroupBox;
  QGroupBox *randomForestGroupBox;
  QGroupBox	*fsGroupBox;

  QGroupBox	*configurationGroupBox;
  QGridLayout *classifierGridLayout;
  QGridLayout *paramsGridLayout;
  QGridLayout* randomForestGridLayout;
  QGridLayout *fsGridLayout;
  QGridLayout *configurationGridLayout;

  QLabel		*cvLabel;
  QLineEdit	*cvValue;

  QLabel    *mSplitModelDirectoryLabel;
  QLineEdit	*mSplitModelDirectory;

  QLabel *parametersLabel;
  QCheckBox *mOptimization;
  QLabel *crossvalidationLabel;

  QLabel* gMinimumLabel;
  QLabel* gMaximumLabel;
  QLabel* cMinimumLabel;
  QLabel* cMaximumLabel;
  QDoubleSpinBox* gMinimumSpinbox;
  QDoubleSpinBox* gMaximumSpinbox;
  QDoubleSpinBox* cMinimumSpinbox;
  QDoubleSpinBox* cMaximumSpinbox;

  QDoubleSpinBox* mConvergenceThresholdSpinbox;

  QLabel* randomForestEpsilonLabel;
  QLabel* randomForestMaxIterationsLabel;
  QDoubleSpinBox* randomForestEpsilonSpinbox;
  QSpinBox* randomForestMaxIterationsSpinbox;

  QProgressBar* progressBar;

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
    classifierGroupBox->setTitle(QString::fromStdString("Classification"));
    classifierGridLayout = new QGridLayout(classifierGroupBox);
    classifierGridLayout->setObjectName(QString::fromUtf8("classifierGridLayout"));
    mLinearKernel = new QRadioButton("SVM: Linear");
    mLinearKernel->setEnabled(true);
    mRBFKernel = new QRadioButton("SVM: RBF");
    mRBFKernel->setEnabled(true);
    mPolynomialKernel = new QRadioButton("SVM: Polynomial");
    mPolynomialKernel->setEnabled(true);
    mSigmoidKernel = new QRadioButton("SVM: Sigmoid");
    mSigmoidKernel->setEnabled(true);
    mChiSquaredKernel = new QRadioButton("SVM: Chi-squared");
    mChiSquaredKernel->setEnabled(true);
    mIntersectionKernel = new QRadioButton("SVM: Histogram Intersection");
    mIntersectionKernel->setEnabled(true);
    mRandomForest = new QRadioButton("Random Forest");
    mRandomForest->setEnabled(true);
    mSGDSVM = new QRadioButton("Stochastic Gradient Descent SVM");
    mSGDSVM->setEnabled(true);
    mBoostedTrees = new QRadioButton("Boosted Trees");
    mBoostedTrees->setEnabled(true);

    mSelectOptimalML = new QRadioButton("Optimized Machine Learning");
    mSelectOptimalML->setToolTip("Iterates over each machine learning strategy and selects the best-performing. WARNING: Multiplies the time required to run!");
    mSelectOptimalML->setEnabled(true);

    classifierGridLayout->addWidget(mLinearKernel, 0, 0, 1, 1);
    classifierGridLayout->addWidget(mRBFKernel, 0, 1, 1, 1);
    classifierGridLayout->addWidget(mPolynomialKernel, 1, 0, 1, 1);
    classifierGridLayout->addWidget(mSigmoidKernel, 1, 1, 1, 1);
    classifierGridLayout->addWidget(mChiSquaredKernel, 2, 0, 1, 1);
    classifierGridLayout->addWidget(mIntersectionKernel, 2, 1, 1, 1);
    classifierGridLayout->addWidget(mRandomForest, 3, 0, 1, 1);
    classifierGridLayout->addWidget(mSGDSVM, 3, 1, 1, 1);
    classifierGridLayout->addWidget(mBoostedTrees, 4, 0, 1, 1);
    classifierGridLayout->addWidget(mSelectOptimalML, 4, 1, 1, 1);

    paramsGroupBox = new QGroupBox(fTrainingSimulator);
    paramsGroupBox->setTitle(QString::fromStdString("Parameter Optimization/Cross-validation"));
    paramsGridLayout = new QGridLayout(paramsGroupBox);
    paramsGridLayout->setObjectName(QString::fromUtf8("paramsGridLayout"));

    randomForestGroupBox = new QGroupBox(fTrainingSimulator);
    randomForestGroupBox->setVisible(false);
    randomForestGroupBox->setTitle(QString::fromStdString("Random Forest Parameters"));
    randomForestGridLayout = new QGridLayout(randomForestGroupBox);
    randomForestGridLayout->setObjectName(QString::fromUtf8("randomForestGridLayout"));

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

    cMinimumSpinbox = new QDoubleSpinBox(paramsGroupBox);
    cMinimumSpinbox->setMaximum(99.9);
    cMinimumSpinbox->setMinimum(-99.9);
    cMinimumSpinbox->setValue(-5.0);
    cMinimumSpinbox->setVisible(false);
    cMaximumSpinbox = new QDoubleSpinBox(paramsGroupBox);
    cMaximumSpinbox->setMaximum(99.9);
    cMaximumSpinbox->setMinimum(-99.9);
    cMaximumSpinbox->setValue(5.0);
    cMaximumSpinbox->setVisible(false);
    gMinimumSpinbox = new QDoubleSpinBox(paramsGroupBox);
    gMinimumSpinbox->setMaximum(99.9);
    gMinimumSpinbox->setMinimum(-99.9);
    gMinimumSpinbox->setValue(-5.0);
    gMinimumSpinbox->setVisible(false);
    gMaximumSpinbox = new QDoubleSpinBox(paramsGroupBox);
    gMaximumSpinbox->setMaximum(99.9);
    gMaximumSpinbox->setMinimum(-99.9);
    gMaximumSpinbox->setValue(5.0);
    gMaximumSpinbox->setVisible(false);

    cMinimumLabel = new QLabel("C Minimum power");
    cMinimumLabel->setVisible(false);
    cMaximumLabel = new QLabel("C Maximum power");
    cMaximumLabel->setVisible(false);
    gMinimumLabel = new QLabel("G Minimum power");
    gMinimumLabel->setVisible(false);
    gMaximumLabel = new QLabel("G Maximum power");
    gMaximumLabel->setVisible(false);

    paramsGridLayout->addWidget(parametersLabel, 0, 0, 1, 1); 
    paramsGridLayout->addWidget(crossvalidationLabel, 0, 2, 1, 1); 
    paramsGridLayout->addWidget(mOptimization, 1, 0, 1, 1);
    paramsGridLayout->addWidget(mResubstitution, 1, 2, 1, 1);
    paramsGridLayout->addWidget(mFiveFold, 1,3, 1, 1);
    
    paramsGridLayout->addWidget(cMinimumLabel, 3, 0, 1, 1);
    paramsGridLayout->addWidget(cMaximumLabel, 4, 0, 1, 1);
    paramsGridLayout->addWidget(cMinimumSpinbox, 3, 1, 1, 1);
    paramsGridLayout->addWidget(cMaximumSpinbox, 4, 1, 1, 1);
    paramsGridLayout->addWidget(gMinimumLabel, 3, 2, 1, 1);
    paramsGridLayout->addWidget(gMaximumLabel, 4, 2, 1, 1);
    paramsGridLayout->addWidget(gMinimumSpinbox, 3, 3, 1, 1);
    paramsGridLayout->addWidget(gMaximumSpinbox, 4, 3, 1, 1);

    randomForestEpsilonLabel = new QLabel("Termination accuracy threshold");
    randomForestEpsilonLabel->setParent(randomForestGroupBox);
    randomForestEpsilonLabel->setToolTip("Train until out-of-bag error is less than this value, unless the maximum iteration count is hit.");
    randomForestEpsilonLabel->setVisible(false);
    randomForestEpsilonSpinbox = new QDoubleSpinBox(randomForestGroupBox);
    randomForestEpsilonSpinbox->setMaximum(0.9999);
    randomForestEpsilonSpinbox->setMinimum(0.0001);
    randomForestEpsilonSpinbox->setSingleStep(0.01);
    randomForestEpsilonSpinbox->setValue(0.1);
    randomForestEpsilonSpinbox->setVisible(false);
    randomForestMaxIterationsLabel = new QLabel("Maximum iterations until termination");
    randomForestMaxIterationsLabel->setToolTip("Increases accuracy but falls off asymptotically. Computation time scales linearly.");
    randomForestMaxIterationsLabel->setVisible(false);
    randomForestMaxIterationsLabel->setParent(randomForestGroupBox);
    randomForestMaxIterationsSpinbox = new QSpinBox(randomForestGroupBox);
    randomForestMaxIterationsSpinbox->setMaximum(500);
    randomForestMaxIterationsSpinbox->setMinimum(1);
    randomForestMaxIterationsSpinbox->setValue(50);
    randomForestMaxIterationsSpinbox->setVisible(false);

    randomForestGridLayout->setObjectName(QString::fromStdString("randomForestGridLayout"));
    randomForestGridLayout->addWidget(randomForestEpsilonLabel, 0, 0, 1, 1);
    randomForestGridLayout->addWidget(randomForestEpsilonSpinbox, 0, 1, 1, 1);
    randomForestGridLayout->addWidget(randomForestMaxIterationsLabel, 1, 0, 1, 1);
    randomForestGridLayout->addWidget(randomForestMaxIterationsSpinbox, 1, 1, 1, 1);


    fsGroupBox = new QGroupBox(fTrainingSimulator);
    fsGroupBox->setTitle(QString::fromStdString("Feature selection"));
    fsGridLayout = new QGridLayout(fsGroupBox);
    fsGridLayout->setObjectName(QString::fromUtf8("fsGridLayout"));
    mSVMFFS = new QRadioButton("Forward feature selection");
    mSVMFFS->setEnabled(true);
    mESFS = new QRadioButton("Effect Size feature selection");
    mESFS->setEnabled(true);
    mBackwardsFS = new QRadioButton("Recursive feature elimination");
    mBackwardsFS->setEnabled(true);
    mRandomForestFS = new QRadioButton("Random Forest importance-based feature selection");
    mReliefFFS = new QRadioButton("RELIEF-F feature selection");
    mReliefFFS->setEnabled(true);
    mSelectOptimalFS = new QRadioButton("Optimized Feature Selection");
    mSelectOptimalFS->setEnabled(true);
    mSelectOptimalFS->setToolTip("Iterates over each feature selection method and selects the best performing. WARNING: multiplies the required run time!");
    fsGridLayout->addWidget(mSVMFFS, 0, 0, 1, 1);
    fsGridLayout->addWidget(mESFS, 0, 1, 1, 1);
    fsGridLayout->addWidget(mBackwardsFS, 0, 2, 1, 1);
    fsGridLayout->addWidget(mRandomForestFS, 0, 3, 1, 1);
    fsGridLayout->addWidget(mReliefFFS, 1, 0, 1, 1);
    fsGridLayout->addWidget(mSelectOptimalFS, 1, 1, 1, 1);

    forwardfsTerminationGroupBox = new QGroupBox(fsGroupBox);
    forwardfsTerminationGroupBox->setTitle("Termination Criteria");
    forwardfsTerminationGridLayout = new QGridLayout(forwardfsTerminationGroupBox);
    mNumberOfFeaturesLabel = new QLabel("Number of features");
    mNumberOfFeaturesLabel->setToolTip("If 0, include the best performing set of features. Otherwise, only include up to this many features. Must be greater than 0 for random forest");

    fsGridLayout->addWidget(forwardfsTerminationGroupBox, 2, 0, -1, -1);
    mNumberOfFeaturesSpinbox = new QSpinBox(fsGroupBox);
    mNumberOfFeaturesSpinbox->setMaximum(99999);
    mNumberOfFeaturesSpinbox->setMinimum(0);
    mNumberOfFeaturesSpinbox->setValue(0);
    mNumberOfFeaturesSpinbox->setToolTip("If 0, include the best performing set of features. Otherwise, only include up to this many features. Must be greater than 0 for random forest.");
    mNumberOfFeaturesSpinbox->setVisible(true);
    mUseFeatureCount = new QRadioButton("Terminate using feature count");
    mUseConvergenceThreshold = new QRadioButton("Terminate using convergence threshold");
    mConvergenceThresholdSpinbox = new QDoubleSpinBox(fsGroupBox);
    mConvergenceThresholdSpinbox->setMaximum(1.0);
    mConvergenceThresholdSpinbox->setMinimum(0.000001);
    mConvergenceThresholdSpinbox->setValue(0.01);
    mConvergenceThresholdSpinbox->setEnabled(false);
    mConvergenceThresholdSpinbox->setVisible(false);

    //forwardfsTerminationGridLayout->addWidget(mUseFeatureCount, 2, 0, 1, 2);
    //forwardfsTerminationGridLayout->addWidget(mUseConvergenceThreshold, 2, 2, 1, 2);
    forwardfsTerminationGridLayout->addWidget(mNumberOfFeaturesLabel, 3, 0, 1, 2);
    forwardfsTerminationGridLayout->addWidget(mNumberOfFeaturesSpinbox, 3, 1, 1, 2);
    //forwardfsTerminationGridLayout->addWidget(mConvergenceThresholdSpinbox, 3, 2, 1, 2);



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
    mSplitTrain = new QRadioButton("Training");
    mSplitTrain->setEnabled(true);
    mSplitTest = new QRadioButton("Testing");
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
    validatePredictionsLabel = new QLabel(configurationGroupBox);
    mTestPredictionsAgainstProvidedLabels = new QCheckBox("Yes/No");
    configurationGridLayout->addWidget(validatePredictionsLabel, 3, 1, 1, 1);
    configurationGridLayout->addWidget(mTestPredictionsAgainstProvidedLabels, 3, 2, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectoryLabel, 4, 0, 1, 1);
    configurationGridLayout->addWidget(mSplitModelDirectory, 4, 1, 1, 2);
    configurationGridLayout->addWidget(mSplitModelDirectoryButton, 4, 3, 1, 1);





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

    progressBar = new QProgressBar(confirmGroupBox);
    confirmGridLayout->addWidget(progressBar, 1, 0, 1, -1);

    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(classifierGroupBox, 2, 0, 1, 2);
    gridLayout->addWidget(paramsGroupBox, 3, 0, 1, 2);
    gridLayout->addWidget(randomForestGroupBox, 4, 0, 1, 2);
    gridLayout->addWidget(fsGroupBox, 5, 0, 1, 2);
    gridLayout->addWidget(configurationGroupBox, 6, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 7, 0, 1, 2);
    gridLayout->addWidget(confirmGroupBox, 8, 0, 1, 2);

    

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
    validatePredictionsLabel->setText(QApplication::translate("fTrainingSimulator", "Validate predictions against given labels", 0));
  } // retranslateUi
};

namespace Ui {
  class fTrainingSimulator : public ui_fTrainingSimulator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fTrainingSimulator_H
