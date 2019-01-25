#ifndef ui_fPseudoProgressionDialog_H
#define ui_fPseudoProgressionDialog_H

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

class ui_fPseudoProgressionDialog
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

  QLabel *longRunningWarning;


  QLineEdit   * existingMaskDirectoryName;
  QPushButton * existingMasksButton;
  QPushButton * disclaimerButton;

  QGridLayout *outputGridLayout;
  QLabel	*outputDirectoryLabel;
  QLabel     *disclaimerLabel;
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

  void setupUi(QDialog *fPseudoProgressionDialog)
  {

    if (fPseudoProgressionDialog->objectName().isEmpty())
      fPseudoProgressionDialog->setObjectName(QString::fromUtf8("fPseudoProgressionDialog"));
    fPseudoProgressionDialog->setWindowModality(Qt::ApplicationModal);
    fPseudoProgressionDialog->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPseudoProgressionDialog->sizePolicy().hasHeightForWidth());
    fPseudoProgressionDialog->setSizePolicy(sizePolicy);
    fPseudoProgressionDialog->setMinimumSize(QSize(0, 0));

    //fPseudoProgressionDialog->setModal(true);
    gridLayout_3 = new QGridLayout(fPseudoProgressionDialog);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //--------------------------------------------------------------------

    classificationGroupBox = new QGroupBox(fPseudoProgressionDialog);
    classificationGroupBox->setTitle(QString::fromStdString("Pseudo-Progression Estimation Modeling"));
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

    disclaimerLabel = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(disclaimerLabel->sizePolicy().hasHeightForWidth());
    disclaimerLabel->setSizePolicy(sizePolicy13);
    disclaimerLabel->setAlignment(Qt::AlignRight);

    disclaimerButton = new QPushButton(classificationGroupBox);
    disclaimerButton->setObjectName(QString::fromUtf8("disclaimerButton"));
    //testSubjectsDirectoryButton->setIcon(ButtonIcon);
    //testSubjectsDirectoryButton->setIconSize(QSize(20, 20));
    disclaimerButton->setText(QString("here."));
    disclaimerButton->setToolTip(QString("disclaimerButton"));
    disclaimerButton->setFlat(true);
    disclaimerButton->setStyleSheet("Text-align:left");


    QPalette* palette1 = new QPalette();
    palette1->setColor(QPalette::ButtonText, Qt::blue);
    disclaimerButton->setPalette(*palette1);

    QFont font("Bavaria");
    font.setPointSize(8);
    font.setWeight(QFont::Bold);
    font.setUnderline(true);
    disclaimerButton->setFont(font);



    //linkLabel = new QLabel(classificationGroupBox);
    //sizePolicy13.setHeightForWidth(disclaimerLabel->sizePolicy().hasHeightForWidth());
    //linkLabel->setSizePolicy(sizePolicy13);
    //linkLabel->setAlignment(Qt::AlignRight);
    //linkLabel->setText("<a href=\"ftp://www.nitrc.org/home/groups/captk/downloads/SampleData_1.6.0/RecurrenceEstimator.zip\">Click Here!</a>");
    //linkLabel->setTextFormat(Qt::RichText);
    //linkLabel->setTextInteractionFlags(Qt:::TextBrowserInteraction);
    //linkLabel->setOpenExternalLinks(true);




    //    classificationGridLayout->addWidget(rdNewClassification, 0, 0, 1, 6);
    classificationGridLayout->addWidget(rdLoadedClassification, 1, 0, 1, 6);
    classificationGridLayout->addWidget(modelDirectoryLabel1, 2, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName1, 2, 1, 1, 4);
    classificationGridLayout->addWidget(svmModelButton1, 2, 5, 1, 1);

    classificationGridLayout->addWidget(rdExistingClassification, 3, 0, 1, 6);
    classificationGridLayout->addWidget(modelDirectoryLabel2, 4, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName2, 4, 1, 1, 4);
    classificationGridLayout->addWidget(svmModelButton2, 4, 5, 1, 1);


    classificationGridLayout->addWidget(testDirectoryLabel, 5, 0, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryName, 5, 1, 1, 4);
    classificationGridLayout->addWidget(testSubjectsDirectoryButton, 5, 5, 1, 1);
    classificationGridLayout->addWidget(disclaimerLabel, 6, 1, 1, 4);
    classificationGridLayout->addWidget(disclaimerButton, 6, 5, 1, 1);


    QFrame* line = new QFrame();
    line->setGeometry(QRect(1, 1, 10, 10));
    line->setFrameShape(QFrame::HLine); // Replace by VLine for vertical line
    line->setFrameShadow(QFrame::Sunken);
    classificationGridLayout->addWidget(line, 8, 0, 1, 6);

    classificationGridLayout->addWidget(rdCreateModel, 9, 0, 1, 6);
    classificationGridLayout->addWidget(trainingDirectoryLabel, 10, 0, 1, 1);
    classificationGridLayout->addWidget(existingMaskDirectoryName, 10, 1, 1, 4);
    classificationGridLayout->addWidget(existingMasksButton, 10, 5, 1, 1);

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


    outputGroupBox = new QGroupBox(fPseudoProgressionDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory"));

    longRunningWarning = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    longRunningWarning->setSizePolicy(sizePolicy);
    longRunningWarning->setAlignment(Qt::AlignRight);
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


    confirmButton = new QPushButton(fPseudoProgressionDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fPseudoProgressionDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout_3->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout_3->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fPseudoProgressionDialog);

    QMetaObject::connectSlotsByName(fPseudoProgressionDialog);
  } // setupUi

  void retranslateUi(QDialog *fPseudoProgressionDialog)
  {
    fPseudoProgressionDialog->setWindowTitle(QApplication::translate("fPseudoProgressionDialog", "Glioblastoma Pseudo-Progression Estimator", 0));
    //rdNewClassification->setText(QApplication::translate("fPseudoProgressionDialog", "Loaded subject: Near&Far", 0));
    rdLoadedClassification->setText(QApplication::translate("fPseudoProgressionDialog", "Pseudo-Progression estimation on loaded subject", 0));
    rdExistingClassification->setText(QApplication::translate("fPseudoProgressionDialog", "Pseudo-Progression estimation on a batch of subjects", 0));
    rdCreateModel->setText(QApplication::translate("fPseudoProgressionDialog", "Train new model", 0));
    confirmButton->setText(QApplication::translate("fPseudoProgressionDialog", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fPseudoProgressionDialog", "Cancel", 0));

    cbT1Data->setText(QApplication::translate("fPseudoProgressionDialog", "Conventional (T1, T2, Flair, T1Gd)", 0));
    cbDTIData->setText(QApplication::translate("fPseudoProgressionDialog", "DTI", 0));
    cbPerfData->setText(QApplication::translate("fPseudoProgressionDialog", "Perfusion", 0));
    cbDistanceData->setText(QApplication::translate("fPseudoProgressionDialog", "Distance", 0));

    modelDirectoryLabel1->setText(QApplication::translate("fPseudoProgressionDialog", "Model Directory:", 0));
    modelDirectoryLabel2->setText(QApplication::translate("fPseudoProgressionDialog", "Model Directory:", 0));
    trainingDirectoryLabel->setText(QApplication::translate("fPseudoProgressionDialog", "Selected Subjects:", 0));
    testDirectoryLabel->setText(QApplication::translate("fPseudoProgressionDialog", "Test Directory:", 0));
    outputDirectoryLabel->setText(QApplication::translate("fPseudoProgressionDialog", "Output Directory:", 0));


  } // retranslateUi

};

namespace Ui {
  class fPseudoProgressionDialog : public ui_fPseudoProgressionDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPseudoProgressionDialog_H






