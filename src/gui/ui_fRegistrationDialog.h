#ifndef ui_fRegistrationDialog_H
#define ui_fRegistrationDialog_H

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

QT_BEGIN_NAMESPACE

class ui_fRegistrationDialog
{
public:
  QGridLayout *gridLayout_3;
  QPushButton * confirmButton;
  QPushButton * cancelButton;

  QLabel *fixedFileLabel;
  QLabel *movingFileLabel1;
  QLabel *movingFileLabel2;
  QLabel *movingFileLabel3;
  QLabel *movingFileLabel4;
  QLabel *movingFileLabel5;

  QLineEdit *fixedFileName;
  QLineEdit *movingFileName1;
  QLineEdit *movingFileName2;
  QLineEdit *movingFileName3;
  QLineEdit *movingFileName4;
  QLineEdit *movingFileName5;

  QPushButton	* fixedFileButton;
  QPushButton *movingFileButton1;
  QPushButton *movingFileButton2;
  QPushButton *movingFileButton3;
  QPushButton *movingFileButton4;
  QPushButton *movingFileButton5;


  QLabel *movingFileOutputLabel1;
  QLabel *movingFileOutputLabel2;
  QLabel *movingFileOutputLabel3;
  QLabel *movingFileOutputLabel4;
  QLabel *movingFileOutputLabel5;
  QLineEdit *movingFileOutputName1;
  QLineEdit *movingFileOutputName2;
  QLineEdit *movingFileOutputName3;
  QLineEdit *movingFileOutputName4;
  QLineEdit *movingFileOutputName5;
  QPushButton *movingFileOutputButton1;
  QPushButton *movingFileOutputButton2;
  QPushButton *movingFileOutputButton3;
  QPushButton *movingFileOutputButton4;
  QPushButton *movingFileOutputButton5;

  //QLabel *longRunningWarning;

  void setupUi(QDialog *fRegistrationDialog)
  {

    if (fRegistrationDialog->objectName().isEmpty())
      fRegistrationDialog->setObjectName(QString::fromUtf8("fRegistrationDialog"));

    fRegistrationDialog->setWindowModality(Qt::ApplicationModal);
    fRegistrationDialog->resize(300, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fRegistrationDialog->sizePolicy().hasHeightForWidth());
    fRegistrationDialog->setSizePolicy(sizePolicy);
    fRegistrationDialog->setMinimumSize(QSize(0, 0));
    QFont font;
    font.setFamily(QString::fromUtf8("Calibri"));
    fRegistrationDialog->setFont(font);
    //fRegistrationDialog->setModal(true);
    gridLayout_3 = new QGridLayout(fRegistrationDialog);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //------------------------------fixed file data-------------------------------
    fixedFileLabel = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(fixedFileLabel->sizePolicy().hasHeightForWidth());
    fixedFileLabel->setSizePolicy(sizePolicy);

    fixedFileName = new QLineEdit(fRegistrationDialog);
    fixedFileName->setObjectName(QString::fromUtf8("existingMaskDirectoryName"));
    sizePolicy.setHeightForWidth(fixedFileName->sizePolicy().hasHeightForWidth());
    fixedFileName->setSizePolicy(sizePolicy);
    fixedFileName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    fixedFileButton = new QPushButton(fRegistrationDialog);
    fixedFileButton->setObjectName(QString::fromUtf8("fixedFileButton"));
    fixedFileButton->setText("Browse");

    //-----------------------------moving file input data---------------------------------
    movingFileLabel1 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileLabel1->sizePolicy().hasHeightForWidth());
    movingFileLabel1->setSizePolicy(sizePolicy);

    movingFileLabel2 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileLabel2->sizePolicy().hasHeightForWidth());
    movingFileLabel2->setSizePolicy(sizePolicy);

    movingFileLabel3 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileLabel3->sizePolicy().hasHeightForWidth());
    movingFileLabel3->setSizePolicy(sizePolicy);

    movingFileLabel4 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileLabel4->sizePolicy().hasHeightForWidth());
    movingFileLabel4->setSizePolicy(sizePolicy);

    movingFileLabel5 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileLabel5->sizePolicy().hasHeightForWidth());
    movingFileLabel5->setSizePolicy(sizePolicy);


    movingFileName1 = new QLineEdit(fRegistrationDialog);
    movingFileName1->setObjectName(QString::fromUtf8("movingFileName1"));
    sizePolicy.setHeightForWidth(movingFileName1->sizePolicy().hasHeightForWidth());
    movingFileName1->setSizePolicy(sizePolicy);
    movingFileName1->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileName2 = new QLineEdit(fRegistrationDialog);
    movingFileName2->setObjectName(QString::fromUtf8("movingFileName2"));
    sizePolicy.setHeightForWidth(movingFileName1->sizePolicy().hasHeightForWidth());
    movingFileName2->setSizePolicy(sizePolicy);
    movingFileName2->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileName3 = new QLineEdit(fRegistrationDialog);
    movingFileName3->setObjectName(QString::fromUtf8("movingFileName3"));
    sizePolicy.setHeightForWidth(movingFileName3->sizePolicy().hasHeightForWidth());
    movingFileName3->setSizePolicy(sizePolicy);
    movingFileName3->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileName4 = new QLineEdit(fRegistrationDialog);
    movingFileName4->setObjectName(QString::fromUtf8("movingFileName1"));
    sizePolicy.setHeightForWidth(movingFileName4->sizePolicy().hasHeightForWidth());
    movingFileName4->setSizePolicy(sizePolicy);
    movingFileName4->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileName5 = new QLineEdit(fRegistrationDialog);
    movingFileName5->setObjectName(QString::fromUtf8("movingFileName1"));
    sizePolicy.setHeightForWidth(movingFileName5->sizePolicy().hasHeightForWidth());
    movingFileName5->setSizePolicy(sizePolicy);
    movingFileName5->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    movingFileButton1 = new QPushButton(fRegistrationDialog);
    movingFileButton1->setObjectName(QString::fromUtf8("movingFileButton1"));
    movingFileButton1->setText("Browse");

    movingFileButton2 = new QPushButton(fRegistrationDialog);
    movingFileButton2->setObjectName(QString::fromUtf8("movingFileButton2"));
    movingFileButton2->setText("Browse");

    movingFileButton3 = new QPushButton(fRegistrationDialog);
    movingFileButton3->setObjectName(QString::fromUtf8("movingFileButton3"));
    movingFileButton3->setText("Browse");

    movingFileButton4 = new QPushButton(fRegistrationDialog);
    movingFileButton4->setObjectName(QString::fromUtf8("movingFileButton4"));
    movingFileButton4->setText("Browse");

    movingFileButton5 = new QPushButton(fRegistrationDialog);
    movingFileButton5->setObjectName(QString::fromUtf8("movingFileButton5"));
    movingFileButton5->setText("Browse");


    //-----------------------------moving file output data---------------------------------
    movingFileOutputLabel1 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileOutputLabel1->sizePolicy().hasHeightForWidth());
    movingFileOutputLabel1->setSizePolicy(sizePolicy);

    movingFileOutputLabel2 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileOutputLabel2->sizePolicy().hasHeightForWidth());
    movingFileOutputLabel2->setSizePolicy(sizePolicy);

    movingFileOutputLabel3 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileOutputLabel3->sizePolicy().hasHeightForWidth());
    movingFileOutputLabel3->setSizePolicy(sizePolicy);

    movingFileOutputLabel4 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileOutputLabel4->sizePolicy().hasHeightForWidth());
    movingFileOutputLabel4->setSizePolicy(sizePolicy);

    movingFileOutputLabel5 = new QLabel(fRegistrationDialog);
    sizePolicy.setHeightForWidth(movingFileOutputLabel5->sizePolicy().hasHeightForWidth());
    movingFileLabel5->setSizePolicy(sizePolicy);


    movingFileOutputName1 = new QLineEdit(fRegistrationDialog);
    movingFileOutputName1->setObjectName(QString::fromUtf8("movingFileOutputName1"));
    sizePolicy.setHeightForWidth(movingFileOutputName1->sizePolicy().hasHeightForWidth());
    movingFileOutputName1->setSizePolicy(sizePolicy);
    movingFileOutputName1->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileOutputName2 = new QLineEdit(fRegistrationDialog);
    movingFileOutputName2->setObjectName(QString::fromUtf8("movingFileOutputName2"));
    sizePolicy.setHeightForWidth(movingFileOutputName1->sizePolicy().hasHeightForWidth());
    movingFileOutputName2->setSizePolicy(sizePolicy);
    movingFileOutputName2->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileOutputName3 = new QLineEdit(fRegistrationDialog);
    movingFileOutputName3->setObjectName(QString::fromUtf8("movingFileOutputName3"));
    sizePolicy.setHeightForWidth(movingFileOutputName3->sizePolicy().hasHeightForWidth());
    movingFileOutputName3->setSizePolicy(sizePolicy);
    movingFileOutputName3->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileOutputName4 = new QLineEdit(fRegistrationDialog);
    movingFileOutputName4->setObjectName(QString::fromUtf8("movingFileOutputName1"));
    sizePolicy.setHeightForWidth(movingFileOutputName4->sizePolicy().hasHeightForWidth());
    movingFileOutputName4->setSizePolicy(sizePolicy);
    movingFileOutputName4->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    movingFileOutputName5 = new QLineEdit(fRegistrationDialog);
    movingFileOutputName5->setObjectName(QString::fromUtf8("movingFileOutputName1"));
    sizePolicy.setHeightForWidth(movingFileOutputName5->sizePolicy().hasHeightForWidth());
    movingFileOutputName5->setSizePolicy(sizePolicy);
    movingFileOutputName5->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    movingFileOutputButton1 = new QPushButton(fRegistrationDialog);
    movingFileOutputButton1->setObjectName(QString::fromUtf8("movingFileOutputButton1"));
    movingFileOutputButton1->setText("Browse");

    movingFileOutputButton2 = new QPushButton(fRegistrationDialog);
    movingFileOutputButton2->setObjectName(QString::fromUtf8("movingFileOutputButton2"));
    movingFileOutputButton2->setText("Browse");

    movingFileOutputButton3 = new QPushButton(fRegistrationDialog);
    movingFileOutputButton3->setObjectName(QString::fromUtf8("movingFileOutputButton3"));
    movingFileOutputButton3->setText("Browse");

    movingFileOutputButton4 = new QPushButton(fRegistrationDialog);
    movingFileOutputButton4->setObjectName(QString::fromUtf8("movingFileOutputButton4"));
    movingFileOutputButton4->setText("Browse");

    movingFileOutputButton5 = new QPushButton(fRegistrationDialog);
    movingFileOutputButton5->setObjectName(QString::fromUtf8("movingFileOutputButton5"));
    movingFileOutputButton5->setText("Browse");

    gridLayout_3->addWidget(fixedFileLabel, 1, 1, 1, 1);
    gridLayout_3->addWidget(fixedFileName, 1, 2, 1, 5);
    gridLayout_3->addWidget(fixedFileButton, 1, 7, 1, 1);

    gridLayout_3->addWidget(movingFileLabel1, 2, 1, 1, 1);
    gridLayout_3->addWidget(movingFileName1, 2, 2, 1, 5);
    gridLayout_3->addWidget(movingFileButton1, 2, 7, 1, 1);

    gridLayout_3->addWidget(movingFileLabel2, 3, 1, 1, 1);
    gridLayout_3->addWidget(movingFileName2, 3, 2, 1, 5);
    gridLayout_3->addWidget(movingFileButton2, 3, 7, 1, 1);

    gridLayout_3->addWidget(movingFileLabel3, 4, 1, 1, 1);
    gridLayout_3->addWidget(movingFileName3, 4, 2, 1, 5);
    gridLayout_3->addWidget(movingFileButton3, 4, 7, 1, 1);

    gridLayout_3->addWidget(movingFileLabel4, 5, 1, 1, 1);
    gridLayout_3->addWidget(movingFileName4, 5, 2, 1, 5);
    gridLayout_3->addWidget(movingFileButton4, 5, 7, 1, 1);

    gridLayout_3->addWidget(movingFileLabel5, 6, 1, 1, 1);
    gridLayout_3->addWidget(movingFileName5, 6, 2, 1, 5);
    gridLayout_3->addWidget(movingFileButton5, 6, 7, 1, 1);






    gridLayout_3->addWidget(movingFileOutputLabel1, 2, 8, 1, 1);
    gridLayout_3->addWidget(movingFileOutputName1, 2, 9, 1, 5);
    gridLayout_3->addWidget(movingFileOutputButton1, 2, 14, 1, 1);

    gridLayout_3->addWidget(movingFileOutputLabel2, 3, 8, 1, 1);
    gridLayout_3->addWidget(movingFileOutputName2, 3, 9, 1, 5);
    gridLayout_3->addWidget(movingFileOutputButton2, 3, 14, 1, 1);

    gridLayout_3->addWidget(movingFileOutputLabel3, 4, 8, 1, 1);
    gridLayout_3->addWidget(movingFileOutputName3, 4, 9, 1, 5);
    gridLayout_3->addWidget(movingFileOutputButton3, 4, 14, 1, 1);

    gridLayout_3->addWidget(movingFileOutputLabel4, 5, 8, 1, 1);
    gridLayout_3->addWidget(movingFileOutputName4, 5, 9, 1, 5);
    gridLayout_3->addWidget(movingFileOutputButton4, 5, 14, 1, 1);

    gridLayout_3->addWidget(movingFileOutputLabel5, 6, 8, 1, 1);
    gridLayout_3->addWidget(movingFileOutputName5, 6, 9, 1, 5);
    gridLayout_3->addWidget(movingFileOutputButton5, 6, 14, 1, 1);

    //longRunningWarning = new QLabel(fRegistrationDialog);
    //sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    //longRunningWarning->setSizePolicy(sizePolicy);
    //longRunningWarning->setAlignment(Qt::AlignRight);
    //longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

    //gridLayout_3->addWidget(longRunningWarning, 7, 0, 1, 15);

    //classificationGroupBox = new QGroupBox(fRegistrationDialog);
    //classificationGroupBox->setTitle(QString::fromStdString("Classification"));
    //classificationGridLayout = new QGridLayout(classificationGroupBox);
    //classificationGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));

    //rdNewClassification = new QRadioButton(classificationGroupBox);
    //rdNewClassification->setObjectName(QString::fromUtf8("newClassification"));
    //rdExistingClassification = new QRadioButton(classificationGroupBox);
    //rdExistingClassification->setObjectName(QString::fromUtf8("existingClassification"));
    //rdCreateModel = new QRadioButton(classificationGroupBox);
    //rdCreateModel->setObjectName(QString::fromUtf8("createModel"));


    //svmModelFileName = new QLineEdit(classificationGroupBox);
    //svmModelFileName->setObjectName(QString::fromUtf8("svmModeFileName"));
    //QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);
    //sizePolicy13.setHorizontalStretch(0);
    //sizePolicy13.setVerticalStretch(0);
    //sizePolicy13.setHeightForWidth(svmModelFileName->sizePolicy().hasHeightForWidth());
    //svmModelFileName->setSizePolicy(sizePolicy13);
    //svmModelFileName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    //testSubjectsDirectoryName = new QLineEdit(classificationGroupBox);
    //testSubjectsDirectoryName->setObjectName(QString::fromUtf8("testSubjectsDirectoryName"));
    //sizePolicy13.setHeightForWidth(testSubjectsDirectoryName->sizePolicy().hasHeightForWidth());
    //testSubjectsDirectoryName->setSizePolicy(sizePolicy13);
    //testSubjectsDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);





    //svmModelButton = new QPushButton(classificationGroupBox);
    //svmModelButton->setObjectName(QString::fromUtf8("svmModelButton"));
    //QPixmap pixmap("../data/images/images/OpenIcon.png");
    //QIcon ButtonIcon(pixmap);
    //svmModelButton->setIcon(ButtonIcon);
    //svmModelButton->setIconSize(QSize(20, 20));




    //existingMasksButton = new QPushButton(classificationGroupBox);
    //existingMasksButton->setObjectName(QString::fromUtf8("existingMasksButton"));
    //existingMasksButton->setIcon(ButtonIcon);
    //existingMasksButton->setIconSize(QSize(20, 20));


    //classificationGridLayout->addWidget(rdNewClassification, 0, 0, 1, 1);
    //classificationGridLayout->addWidget(rdExistingClassification, 1, 0, 1, 1);
    //classificationGridLayout->addWidget(svmModelFileName, 1, 1, 1, 1);
    //classificationGridLayout->addWidget(svmModelButton, 1, 2, 1, 1);
    //classificationGridLayout->addWidget(testSubjectsDirectoryName, 2, 1, 1, 1);
    //classificationGridLayout->addWidget(testSubjectsDirectoryButton, 2, 2, 1, 1);

    //classificationGridLayout->addWidget(rdCreateModel, 3, 0, 1, 1);
    //classificationGridLayout->addWidget(existingMaskDirectoryName, 3, 1, 1, 1);
    //classificationGridLayout->addWidget(existingMasksButton, 3, 2, 1, 1);

    //line_3 = new QFrame(classificationGroupBox);
    //line_3->setFrameStyle(QFrame::HLine);
    //line_3->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    //line_3->setFrameShadow(QFrame::Sunken);

    //classificationGridLayout->addWidget(line_3, 4,1, 1, 2);

    //cbT1Data = new QCheckBox(classificationGroupBox);
    //cbT1Data->setObjectName(QString::fromUtf8("cbT1Data"));
    //cbT2Data = new QCheckBox(classificationGroupBox);
    //cbT2Data->setObjectName(QString::fromUtf8("cbT2Data"));
    //cbT1ceData = new QCheckBox(classificationGroupBox);
    //cbT1ceData->setObjectName(QString::fromUtf8("cbT1ceData"));
    //cbT2FlairData = new QCheckBox(classificationGroupBox);
    //cbT2FlairData->setObjectName(QString::fromUtf8("cbT2FlairData"));
    //cbDTIData = new QCheckBox(classificationGroupBox);
    //cbDTIData->setObjectName(QString::fromUtf8("cbDTIData"));
    //cbPerfData = new QCheckBox(classificationGroupBox);
    //cbPerfData->setObjectName(QString::fromUtf8("cbPerfData"));
    //cbDistanceData = new QCheckBox(classificationGroupBox);
    //cbDistanceData->setObjectName(QString::fromUtf8("cbDistanceData"));


    //horizontalLayout = new QHBoxLayout();
    //horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
    //
    //horizontalLayout->addWidget(cbT1Data);
    //horizontalLayout->addWidget(cbT2Data);
    //horizontalLayout->addWidget(cbT1ceData);
    //horizontalLayout->addWidget(cbT2FlairData);
    //horizontalLayout->addWidget(cbDTIData);
    //horizontalLayout->addWidget(cbPerfData);
    //horizontalLayout->addWidget(cbDistanceData);

    //classificationGridLayout->addLayout(horizontalLayout, 5,0, 1, 3);
    ////--------------------------output-------------------------------------------------------


    //outputGroupBox = new QGroupBox(fRegistrationDialog);
    //outputGroupBox->setTitle(QString::fromStdString("Output"));

    //outputGridLayout = new QGridLayout(outputGroupBox);
    //outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputDirectoryLabel = new QLabel(outputGroupBox);
    //sizePolicy13.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    //outputDirectoryLabel->setSizePolicy(sizePolicy13);

    //outputDirectoryName = new QLineEdit(outputGroupBox);
    //outputDirectoryName->setObjectName(QString::fromUtf8("outputDirectoryName"));
    //sizePolicy13.setHeightForWidth(outputDirectoryName->sizePolicy().hasHeightForWidth());
    //outputDirectoryName->setSizePolicy(sizePolicy13);
    //outputDirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    //outputDirectoryButton = new QPushButton(outputGroupBox);
    //outputDirectoryButton->setObjectName(QString::fromUtf8("outputDirectoryButton"));
    //outputDirectoryButton->setIcon(ButtonIcon);
    //outputDirectoryButton->setIconSize(QSize(20, 20));

    //outputGridLayout->addWidget(outputDirectoryLabel, 0, 0, 1, 1);
    //outputGridLayout->addWidget(outputDirectoryName, 1, 0, 1, 1);
    //outputGridLayout->addWidget(outputDirectoryButton, 1, 1, 1, 1);

    //gridLayout_3->addWidget(classificationGroupBox, 0, 0, 1, 2);
    //gridLayout_3->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fRegistrationDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fRegistrationDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout_3->addWidget(confirmButton, 7, 2, 1, 1);
    gridLayout_3->addWidget(cancelButton, 7, 3, 1, 1);

    //gridLayout_3->addWidget(confirmButton,2,0,1,1);
    //gridLayout_3->addWidget(cancelButton,2,1,1,1);

    retranslateUi(fRegistrationDialog);

    QMetaObject::connectSlotsByName(fRegistrationDialog);
  } // setupUi

  void retranslateUi(QDialog *fRegistrationDialog)
  {

    fixedFileLabel->setText(QApplication::translate("fRegistrationDialog", "Target file:", 0, QApplication::UnicodeUTF8));
    movingFileLabel1->setText(QApplication::translate("fRegistrationDialog", "Source file 1:", 0, QApplication::UnicodeUTF8));
    movingFileLabel2->setText(QApplication::translate("fRegistrationDialog", "Source file 2:", 0, QApplication::UnicodeUTF8));
    movingFileLabel3->setText(QApplication::translate("fRegistrationDialog", "Source file 3:", 0, QApplication::UnicodeUTF8));
    movingFileLabel4->setText(QApplication::translate("fRegistrationDialog", "Source file 4:", 0, QApplication::UnicodeUTF8));
    movingFileLabel5->setText(QApplication::translate("fRegistrationDialog", "Source file 5:", 0, QApplication::UnicodeUTF8));

    movingFileOutputLabel1->setText(QApplication::translate("fRegistrationDialog", "Output file 1:", 0, QApplication::UnicodeUTF8));
    movingFileOutputLabel2->setText(QApplication::translate("fRegistrationDialog", "Output file 2:", 0, QApplication::UnicodeUTF8));
    movingFileOutputLabel3->setText(QApplication::translate("fRegistrationDialog", "Output file 3:", 0, QApplication::UnicodeUTF8));
    movingFileOutputLabel4->setText(QApplication::translate("fRegistrationDialog", "Output file 4:", 0, QApplication::UnicodeUTF8));
    movingFileOutputLabel5->setText(QApplication::translate("fRegistrationDialog", "Output file 5:", 0, QApplication::UnicodeUTF8));

    cancelButton->setText(QApplication::translate("fRegistrationDialog", "Cancel", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fRegistrationDialog", "Confirm", 0, QApplication::UnicodeUTF8));

  } // retranslateUi

};

namespace Ui {
  class fRegistrationDialog : public ui_fRegistrationDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fRegistrationDialog_H