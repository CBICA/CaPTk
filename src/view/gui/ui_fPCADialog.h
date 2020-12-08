#ifndef ui_fPCADialog_H
#define ui_fPCADialog_H

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

class ui_fPCADialog
{
public:
  //fDisclaimerDialog disclaimer_dialog;
  QGridLayout * gridLayout_3;
  QVBoxLayout * verticalLayout;
  QGroupBox *classificationGroupBox;
  QGroupBox *trainingGroupBox;
  QGridLayout * classificationGridLayout;
  QGroupBox *outputGroupBox;
  //QRadioButton * rdNewClassification;
  QRadioButton * rdLoadedClassification;
  QRadioButton * rdExistingClassification;
  QRadioButton *rdCreateModel;
  QLineEdit * svmModelFileName1;
  QLineEdit * svmModelFileName2;
  QPushButton * svmModelButton1;
  QPushButton * svmModelButton2;
  QLineEdit * testSubjectsDirectoryName;
  QPushButton * testSubjectsDirectoryButton;

  QLabel *numberPCALabelCreate; 
  QLineEdit * numberPCANameCreate;
  QLabel *numberPCALabelApply;
  QLineEdit * numberPCANameApply;


  QLabel *longRunningWarning;


  QLineEdit   * existingMaskDirectoryName;
  QPushButton * existingMasksButton;
  QPushButton * disclaimerButton;

  QGridLayout *outputGridLayout;
  QLabel	*outputDirectoryLabel;
  //QLabel     *disclaimerLabel;
  QLabel     *linkLabel;
  QLabel	*trainingDirectoryLabel;
  QLabel	*testDirectoryLabel;
  QLabel	*modelDirectoryLabel1;
  QLabel	*modelDirectoryLabel2;


  QLineEdit *outputDirectoryName;
  QPushButton *outputDirectoryButton;

  QPushButton * confirmButton;
  QPushButton * cancelButton;

  QCheckBox * cbT1Data;
  QCheckBox * cbDTIData;
  QCheckBox * cbPerfData;
  QCheckBox * cbDistanceData;

  QHBoxLayout * horizontalLayout;

  QGridLayout *overallGridLayout, *optionalgroupBoxLayout;
  QRadioButton *extractPCA, *applyPCA;
  QLabel *inputDirLabel, *outputDirLabel, *nPCsLabel, *varLabel, *pcaParamsLabel;
  QLineEdit *inputDirLE, *outputDirLE, *nPCsLE, *varLE, *pcaParamsLE;
  QPushButton *pbInput, *pbOutput, *pbParams;
  QGroupBox *optionalGroupBox;

  void setupUi(QDialog *fPCADialog)
  {

    if (fPCADialog->objectName().isEmpty())
      fPCADialog->setObjectName(QString::fromUtf8("fPCADialog"));
    fPCADialog->setWindowModality(Qt::ApplicationModal);
    fPCADialog->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPCADialog->sizePolicy().hasHeightForWidth());
    fPCADialog->setSizePolicy(sizePolicy);
    fPCADialog->setMinimumSize(QSize(0, 0));

	overallGridLayout = new QGridLayout(fPCADialog);
	extractPCA = new QRadioButton("Extract new parameters",fPCADialog);
	applyPCA = new QRadioButton("Apply extracted Parameters",fPCADialog);
	inputDirLabel = new QLabel("Input Directory",fPCADialog);
	outputDirLabel = new QLabel("Output Directory",fPCADialog);
	optionalGroupBox = new QGroupBox("Optional Parameters", fPCADialog);
	optionalgroupBoxLayout = new QGridLayout(optionalGroupBox);
	//TBD: add line here
	//TBD:: add text 'optional parameters'
	nPCsLabel = new QLabel("Number of PCA images to produce",fPCADialog);
	varLabel = new QLabel("Variance threshold",fPCADialog);
	pcaParamsLabel = new QLabel("Extracted parameters directory", fPCADialog);
	inputDirLE = new QLineEdit(fPCADialog);
	outputDirLE = new QLineEdit(fPCADialog);
	nPCsLE = new QLineEdit(fPCADialog);
	varLE = new QLineEdit(fPCADialog);
	pcaParamsLE = new QLineEdit(fPCADialog);
	pbInput = new QPushButton("Browse",fPCADialog);
	pbOutput = new QPushButton("Browse",fPCADialog);
	pbParams = new QPushButton("Browse",fPCADialog);

	overallGridLayout->addWidget(extractPCA, 0, 0/*, 1, 2*/);
	overallGridLayout->addWidget(applyPCA, 0, 1/*, 1, 2*/);
	overallGridLayout->addWidget(inputDirLabel, 1, 0);
	overallGridLayout->addWidget(inputDirLE, 1, 1);
	overallGridLayout->addWidget(pbInput, 1, 2);
	overallGridLayout->addWidget(outputDirLabel, 2, 0);
	overallGridLayout->addWidget(outputDirLE, 2, 1);
	overallGridLayout->addWidget(pbOutput, 2, 2);
	//TBD: in case of apply, this(pcaParamsLabel) should appear on the top ( do it if easy enough)
	overallGridLayout->addWidget(pcaParamsLabel, 3, 0); 
	overallGridLayout->addWidget(pcaParamsLE, 3, 1);
	overallGridLayout->addWidget(pbParams, 3, 2);
	
	overallGridLayout->addWidget(optionalGroupBox, 4, 0, 1, 4);
	optionalgroupBoxLayout->addWidget(nPCsLabel, 0, 0, 1, 2);
	optionalgroupBoxLayout->addWidget(nPCsLE, 0, 2, 1, 1);
	optionalgroupBoxLayout->addWidget(varLabel, 1, 0, 1, 2);
	optionalgroupBoxLayout->addWidget(varLE, 1, 2, 1, 1);

    //fPCADialog->setModal(true);
    gridLayout_3 = new QGridLayout(/*fPCADialog*/);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //--------------------------------------------------------------------

    classificationGroupBox = new QGroupBox(/*fPCADialog*/);
    classificationGroupBox->setTitle(QString::fromStdString("PCA Estimation"));
    classificationGridLayout = new QGridLayout(classificationGroupBox);
    classificationGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));

    //rdNewClassification = new QRadioButton(classificationGroupBox);
    //rdNewClassification->setObjectName(QString::fromUtf8("newClassification"));

    rdLoadedClassification = new QRadioButton(classificationGroupBox);
    rdLoadedClassification->setObjectName(QString::fromUtf8("loadedClassification"));
    rdExistingClassification = new QRadioButton(classificationGroupBox);
    rdExistingClassification->setObjectName(QString::fromUtf8("existingClassification"));
    rdCreateModel = new QRadioButton(classificationGroupBox);
    rdCreateModel->setObjectName(QString::fromUtf8("createModel"));



    svmModelFileName1 = new QLineEdit(classificationGroupBox);
    svmModelFileName1->setObjectName(QString::fromUtf8("svmModeFileName1"));
    QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy13.setHorizontalStretch(0);
    sizePolicy13.setVerticalStretch(0);
    sizePolicy13.setHeightForWidth(svmModelFileName1->sizePolicy().hasHeightForWidth());
    svmModelFileName1->setSizePolicy(sizePolicy13);
    svmModelFileName1->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    svmModelFileName1->setText("");


    svmModelFileName2 = new QLineEdit(classificationGroupBox);
    svmModelFileName2->setObjectName(QString::fromUtf8("svmModeFileName2"));
    sizePolicy13.setHeightForWidth(svmModelFileName2->sizePolicy().hasHeightForWidth());
    svmModelFileName2->setSizePolicy(sizePolicy13);
    svmModelFileName2->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    svmModelFileName2->setText("");






    testSubjectsDirectoryName = new QLineEdit(classificationGroupBox);
    testSubjectsDirectoryName->setObjectName(QString::fromUtf8("testSubjectsDirectoryName"));
    sizePolicy13.setHeightForWidth(testSubjectsDirectoryName->sizePolicy().hasHeightForWidth());
    testSubjectsDirectoryName->setSizePolicy(sizePolicy13);
    testSubjectsDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

	//this is the blank line edit field underneath the selected subjects
    existingMaskDirectoryName = new QLineEdit(classificationGroupBox);
    existingMaskDirectoryName->setObjectName(QString::fromUtf8("existingMaskDirectoryName"));
    sizePolicy13.setHeightForWidth(existingMaskDirectoryName->sizePolicy().hasHeightForWidth());
    existingMaskDirectoryName->setSizePolicy(sizePolicy13);
    existingMaskDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    svmModelButton1 = new QPushButton(classificationGroupBox);
    svmModelButton1->setObjectName(QString::fromUtf8("svmModelButton1"));
    svmModelButton1->setText(QString("Browse"));

    svmModelButton2 = new QPushButton(classificationGroupBox);
    svmModelButton2->setObjectName(QString::fromUtf8("svmModelButton2"));
    svmModelButton2->setText(QString("Browse"));

    testSubjectsDirectoryButton = new QPushButton(classificationGroupBox);
    testSubjectsDirectoryButton->setObjectName(QString::fromUtf8("testSubjectsDirectoryButton"));
    //testSubjectsDirectoryButton->setIcon(ButtonIcon);
    //testSubjectsDirectoryButton->setIconSize(QSize(20, 20));
    testSubjectsDirectoryButton->setText(QString("Browse"));
    testSubjectsDirectoryButton->setToolTip(QString("Directory containing Test subjects"));

	//where is this used?
    existingMasksButton = new QPushButton(classificationGroupBox);
    existingMasksButton->setObjectName(QString::fromUtf8("existingMasksButton"));
    //existingMasksButton->setIcon(ButtonIcon);
    //existingMasksButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    existingMasksButton->setText(QString("Browse"));
    // existingMasksButton->setToolTip(QString("Directory containing Training subjects"));
    existingMasksButton->setWhatsThis(QString("&nbsp;The meaning of the Source field depends "
      "on the Type field:"
      "<ul>"
      "<li><b>Books</b> have a Publisher"
      "<li><b>Articles</b> have a Journal name with "
      "volume and issue number"
      "<li><b>Theses</b> have an Institution name "
      "and a Department name"
      "</ul>"));

    trainingDirectoryLabel = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(trainingDirectoryLabel->sizePolicy().hasHeightForWidth());
    trainingDirectoryLabel->setSizePolicy(sizePolicy13);
    trainingDirectoryLabel->setAlignment(Qt::AlignRight);

    testDirectoryLabel = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(testDirectoryLabel->sizePolicy().hasHeightForWidth());
    testDirectoryLabel->setSizePolicy(sizePolicy13);
    testDirectoryLabel->setAlignment(Qt::AlignRight);

    modelDirectoryLabel1 = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(modelDirectoryLabel1->sizePolicy().hasHeightForWidth());
    modelDirectoryLabel1->setSizePolicy(sizePolicy13);
    modelDirectoryLabel1->setAlignment(Qt::AlignRight);

    modelDirectoryLabel2 = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(modelDirectoryLabel2->sizePolicy().hasHeightForWidth());
    modelDirectoryLabel2->setSizePolicy(sizePolicy13);
    modelDirectoryLabel2->setAlignment(Qt::AlignRight);

    //disclaimerLabel = new QLabel(classificationGroupBox);
    //sizePolicy13.setHeightForWidth(disclaimerLabel->sizePolicy().hasHeightForWidth());
    //disclaimerLabel->setSizePolicy(sizePolicy13);
    //disclaimerLabel->setAlignment(Qt::AlignRight);

    //disclaimerButton = new QPushButton(classificationGroupBox);
    //disclaimerButton->setObjectName(QString::fromUtf8("disclaimerButton"));
    ////testSubjectsDirectoryButton->setIcon(ButtonIcon);
    ////testSubjectsDirectoryButton->setIconSize(QSize(20, 20));
    //disclaimerButton->setText(QString("here."));
    //disclaimerButton->setToolTip(QString("disclaimerButton"));
    //disclaimerButton->setFlat(true);
    //disclaimerButton->setStyleSheet("Text-align:left");


    //QPalette* palette1 = new QPalette();
    //palette1->setColor(QPalette::ButtonText, Qt::blue);
    //disclaimerButton->setPalette(*palette1);

    //QFont font("Bavaria");
    //font.setPointSize(8);
    //font.setWeight(QFont::Bold);
    //font.setUnderline(true);
    //disclaimerButton->setFont(font);



    //linkLabel = new QLabel(classificationGroupBox);
    //sizePolicy13.setHeightForWidth(disclaimerLabel->sizePolicy().hasHeightForWidth());
    //linkLabel->setSizePolicy(sizePolicy13);
    //linkLabel->setAlignment(Qt::AlignRight);
    //linkLabel->setText("<a href=\"ftp://www.nitrc.org/home/groups/captk/downloads/SampleData_1.6.0/RecurrenceEstimator.zip\">Click Here!</a>");
    //linkLabel->setTextFormat(Qt::RichText);
    //linkLabel->setTextInteractionFlags(Qt:::TextBrowserInteraction);
    //linkLabel->setOpenExternalLinks(true);

    numberPCALabelApply = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(numberPCALabelApply->sizePolicy().hasHeightForWidth());
    numberPCALabelApply->setSizePolicy(sizePolicy13);
    numberPCALabelApply->setAlignment(Qt::AlignRight);


    numberPCANameApply = new QLineEdit(classificationGroupBox);
    numberPCANameApply->setObjectName(QString::fromUtf8("numberPCANameApply"));
    sizePolicy13.setHeightForWidth(numberPCANameApply->sizePolicy().hasHeightForWidth());
    numberPCANameApply->setSizePolicy(sizePolicy13);
    numberPCANameApply->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    numberPCALabelCreate = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(numberPCALabelCreate->sizePolicy().hasHeightForWidth());
    numberPCALabelCreate->setSizePolicy(sizePolicy13);
    numberPCALabelCreate->setAlignment(Qt::AlignRight);

    numberPCANameCreate = new QLineEdit(classificationGroupBox);
    numberPCANameCreate->setObjectName(QString::fromUtf8("numberPCANameCreate"));
    sizePolicy13.setHeightForWidth(numberPCANameCreate->sizePolicy().hasHeightForWidth());
    numberPCANameCreate->setSizePolicy(sizePolicy13);
    numberPCANameCreate->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    //    classificationGridLayout->addWidget(rdNewClassification, 0, 0, 1, 6);
    //classificationGridLayout->addWidget(rdLoadedClassification, 1, 0, 1, 6);
    //classificationGridLayout->addWidget(modelDirectoryLabel1, 2, 0, 1, 1);
    //classificationGridLayout->addWidget(svmModelFileName1, 2, 1, 1, 4);
    //classificationGridLayout->addWidget(svmModelButton1, 2, 5, 1, 1);

    classificationGridLayout->addWidget(rdExistingClassification, 3, 0, 1, 6);
    classificationGridLayout->addWidget(modelDirectoryLabel2, 4, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName2, 4, 1, 1, 4);
    classificationGridLayout->addWidget(svmModelButton2, 4, 5, 1, 1);


    classificationGridLayout->addWidget(testDirectoryLabel, 5, 0, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryName, 5, 1, 1, 4);
    classificationGridLayout->addWidget(testSubjectsDirectoryButton, 5, 5, 1, 1);
    classificationGridLayout->addWidget(numberPCALabelApply, 6, 0, 1, 1);
    classificationGridLayout->addWidget(numberPCANameApply, 6, 1, 1, 4);

 /*   classificationGridLayout->addWidget(disclaimerLabel, 7, 1, 1, 4);
    classificationGridLayout->addWidget(disclaimerButton, 7, 5, 1, 1);*/


    QFrame* line = new QFrame();
    line->setGeometry(QRect(1, 1, 10, 10));
    line->setFrameShape(QFrame::HLine); // Replace by VLine for vertical line
    line->setFrameShadow(QFrame::Sunken);
    classificationGridLayout->addWidget(line, 9, 0, 1, 6);

    classificationGridLayout->addWidget(rdCreateModel, 10, 0, 1, 6);
    classificationGridLayout->addWidget(trainingDirectoryLabel, 11, 0, 1, 1);
    classificationGridLayout->addWidget(existingMaskDirectoryName, 11, 1, 1, 4);
    classificationGridLayout->addWidget(existingMasksButton, 11, 5, 1, 1);

    classificationGridLayout->addWidget(numberPCALabelCreate, 12, 0, 1, 1);
    classificationGridLayout->addWidget(numberPCANameCreate, 12, 1, 1, 4);


    // classificationGridLayout->addWidget(line_3, 4, 1, 1, 2);

    cbT1Data = new QCheckBox(classificationGroupBox);
    cbT1Data->setObjectName(QString::fromUtf8("cbT1Data"));
    cbDTIData = new QCheckBox(classificationGroupBox);
    cbDTIData->setObjectName(QString::fromUtf8("cbDTIData"));
    cbPerfData = new QCheckBox(classificationGroupBox);
    cbPerfData->setObjectName(QString::fromUtf8("cbPerfData"));
    cbDistanceData = new QCheckBox(classificationGroupBox);
    cbDistanceData->setObjectName(QString::fromUtf8("cbDistanceData"));


    horizontalLayout = new QHBoxLayout();
    horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));

    horizontalLayout->addWidget(cbT1Data);
    //horizontalLayout->addWidget(cbT2Data);
    //horizontalLayout->addWidget(cbT1ceData);
    //horizontalLayout->addWidget(cbT2FlairData);
    horizontalLayout->addWidget(cbDTIData);
    horizontalLayout->addWidget(cbPerfData);
    horizontalLayout->addWidget(cbDistanceData);

    classificationGridLayout->addLayout(horizontalLayout, 8, 0, 1, 3);
    //--------------------------output-------------------------------------------------------


    outputGroupBox = new QGroupBox(/*fPCADialog*/);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory"));

    longRunningWarning = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    longRunningWarning->setSizePolicy(sizePolicy);
    //longRunningWarning->setAlignment(Qt::AlignRight);
    longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputDirectoryLabel = new QLabel(outputGroupBox);
    sizePolicy13.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    outputDirectoryLabel->setSizePolicy(sizePolicy13);
    outputDirectoryLabel->setAlignment(Qt::AlignRight);




    outputDirectoryName = new QLineEdit(outputGroupBox);
    outputDirectoryName->setObjectName(QString::fromUtf8("outputDirectoryName"));
    sizePolicy13.setHeightForWidth(outputDirectoryName->sizePolicy().hasHeightForWidth());
    outputDirectoryName->setSizePolicy(sizePolicy13);
    outputDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputDirectoryButton = new QPushButton(outputGroupBox);
    outputDirectoryButton->setObjectName(QString::fromUtf8("outputDirectoryButton"));
    //outputDirectoryButton->setIcon(ButtonIcon);
    //outputDirectoryButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    outputDirectoryButton->setText(QString("Browse"));

    outputGridLayout->addWidget(outputDirectoryLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirectoryName, 0, 1, 1, 4);
    outputGridLayout->addWidget(outputDirectoryButton, 0, 5, 1, 1);
    outputGridLayout->addWidget(longRunningWarning, 1, 0, 1, 2);

    gridLayout_3->addWidget(classificationGroupBox, 0, 0, 1, 2);
    gridLayout_3->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fPCADialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fPCADialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout_3->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout_3->addWidget(cancelButton, 2, 1, 1, 1);

	overallGridLayout->addWidget(longRunningWarning, 5, 0, 1, 4);
	overallGridLayout->addWidget(confirmButton, 6, 0);
	overallGridLayout->addWidget(cancelButton, 6, 1);

    retranslateUi(fPCADialog);

    QMetaObject::connectSlotsByName(fPCADialog);
  } // setupUi

  void retranslateUi(QDialog *fPCADialog)
  {
    fPCADialog->setWindowTitle(QApplication::translate("fPCADialog", "PCA Parameter Extractor", 0));
    //rdNewClassification->setText(QApplication::translate("fPCADialog", "Loaded subject: Near&Far", 0));
    rdExistingClassification->setText(QApplication::translate("fPCADialog", "Apply PCA", 0));
    rdCreateModel->setText(QApplication::translate("fPCADialog", "Extract PCA Parameters", 0));
    confirmButton->setText(QApplication::translate("fPCADialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fPCADialog", "Cancel", 0));

    modelDirectoryLabel1->setText(QApplication::translate("fPCADialog", "Model Directory:", 0));
    modelDirectoryLabel2->setText(QApplication::translate("fPCADialog", "Model Directory:", 0));
    trainingDirectoryLabel->setText(QApplication::translate("fPCADialog", "Selected Subjects:", 0));
    testDirectoryLabel->setText(QApplication::translate("fPCADialog", "Test Directory:", 0));
    outputDirectoryLabel->setText(QApplication::translate("fPCADialog", "Output Directory:", 0));


  } // retranslateUi

};

namespace Ui {
  class fPCADialog : public ui_fPCADialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPCADialog_H






