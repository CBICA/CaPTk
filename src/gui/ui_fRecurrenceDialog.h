#ifndef ui_fRecurrenceDialog_H
#define ui_fRecurrenceDialog_H

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

class ui_fRecurrenceDialog
{
public:
  QGridLayout *gridLayout_3;
  QVBoxLayout * verticalLayout;
  QGroupBox *classificationGroupBox;
  QGridLayout * classificationGridLayout;
  QGroupBox *outputGroupBox;
  QRadioButton * rdNewClassification;
  QRadioButton * rdLoadedClassification;
  QRadioButton * rdExistingClassification;
  QRadioButton *rdCreateModel;
  QLineEdit * svmModelFileName;
  QPushButton * svmModelButton;
  QLineEdit * testSubjectsDirectoryName;
  QPushButton * testSubjectsDirectoryButton;


  QLineEdit * existingMaskDirectoryName;
  QPushButton * existingMasksButton;

  QGridLayout *outputGridLayout;
  QLabel	*outputDirectoryLabel;
  QLabel	*trainingDirectoryLabel;
  QLabel	*testDirectoryLabel;
  QLabel	*modelDirectoryLabel;


  QLineEdit *outputDirectoryName;
  QPushButton *outputDirectoryButton;

  QPushButton * confirmButton;
  QPushButton * cancelButton;

  QCheckBox * cbT1Data;
  QCheckBox * cbDTIData;
  QCheckBox * cbPerfData;
  QCheckBox * cbDistanceData;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fRecurrenceDialog)
  {

    if (fRecurrenceDialog->objectName().isEmpty())
      fRecurrenceDialog->setObjectName(QString::fromUtf8("fRecurrenceDialog"));
    fRecurrenceDialog->setWindowModality(Qt::ApplicationModal);
    fRecurrenceDialog->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fRecurrenceDialog->sizePolicy().hasHeightForWidth());
    fRecurrenceDialog->setSizePolicy(sizePolicy);
    fRecurrenceDialog->setMinimumSize(QSize(0, 0));

    fRecurrenceDialog->setModal(true);
    gridLayout_3 = new QGridLayout(fRecurrenceDialog);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //--------------------------------------------------------------------
    classificationGroupBox = new QGroupBox(fRecurrenceDialog);
    classificationGroupBox->setTitle(QString::fromStdString("Classification"));
    classificationGridLayout = new QGridLayout(classificationGroupBox);
    classificationGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));

    rdNewClassification = new QRadioButton(classificationGroupBox);
    rdNewClassification->setObjectName(QString::fromUtf8("newClassification"));

	  rdLoadedClassification = new QRadioButton(classificationGroupBox);
	  rdLoadedClassification->setObjectName(QString::fromUtf8("loadedClassification"));
    rdExistingClassification = new QRadioButton(classificationGroupBox);
    rdExistingClassification->setObjectName(QString::fromUtf8("existingClassification"));
    rdCreateModel = new QRadioButton(classificationGroupBox);
    rdCreateModel->setObjectName(QString::fromUtf8("createModel"));



    svmModelFileName = new QLineEdit(classificationGroupBox);
    svmModelFileName->setObjectName(QString::fromUtf8("svmModeFileName"));
    QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy13.setHorizontalStretch(0);
    sizePolicy13.setVerticalStretch(0);
    sizePolicy13.setHeightForWidth(svmModelFileName->sizePolicy().hasHeightForWidth());
    svmModelFileName->setSizePolicy(sizePolicy13);
    svmModelFileName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);
    svmModelFileName->setText("../data/recurrence/");

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

    svmModelButton = new QPushButton(classificationGroupBox);
    svmModelButton->setObjectName(QString::fromUtf8("svmModelButton"));
    //QPixmap pixmap("../data/images/images/OpenIcon.png");
    //QIcon ButtonIcon(pixmap);
    //svmModelButton->setIcon(ButtonIcon);
    //svmModelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    svmModelButton->setText(QString("Browse"));

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

    modelDirectoryLabel = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(modelDirectoryLabel->sizePolicy().hasHeightForWidth());
    modelDirectoryLabel->setSizePolicy(sizePolicy13);
    modelDirectoryLabel->setAlignment(Qt::AlignRight);

    classificationGridLayout->addWidget(rdNewClassification, 0, 0, 1, 6);
	  classificationGridLayout->addWidget(rdLoadedClassification, 1, 0, 1, 6);
    classificationGridLayout->addWidget(rdExistingClassification, 2, 0, 1, 6);
    classificationGridLayout->addWidget(modelDirectoryLabel, 3, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName, 3, 1, 1, 4);
    classificationGridLayout->addWidget(svmModelButton, 3, 5, 1, 1);

    classificationGridLayout->addWidget(testDirectoryLabel, 4, 0, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryName, 4, 1, 1, 4);
    classificationGridLayout->addWidget(testSubjectsDirectoryButton, 4, 5, 1, 1);

    classificationGridLayout->addWidget(rdCreateModel, 5, 0, 1, 6);
    classificationGridLayout->addWidget(trainingDirectoryLabel, 6, 0, 1, 1);
    classificationGridLayout->addWidget(existingMaskDirectoryName, 6, 1, 1, 4);
    classificationGridLayout->addWidget(existingMasksButton, 6, 5, 1, 1);

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

    classificationGridLayout->addLayout(horizontalLayout, 5, 0, 1, 3);
    //--------------------------output-------------------------------------------------------


    outputGroupBox = new QGroupBox(fRecurrenceDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output Directory"));

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

    gridLayout_3->addWidget(classificationGroupBox, 0, 0, 1, 2);
    gridLayout_3->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fRecurrenceDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fRecurrenceDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout_3->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout_3->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fRecurrenceDialog);

    QMetaObject::connectSlotsByName(fRecurrenceDialog);
  } // setupUi

  void retranslateUi(QDialog *fRecurrenceDialog)
  {
    fRecurrenceDialog->setWindowTitle(QApplication::translate("fRecurrenceDialog", "Glioblastoma Infiltration Index", 0, QApplication::UnicodeUTF8));
    rdNewClassification->setText(QApplication::translate("fRecurrenceDialog", "Loaded subject: Near&Far", 0, QApplication::UnicodeUTF8));
	  rdLoadedClassification->setText(QApplication::translate("fRecurrenceDialog", "Loaded subject: Existing model", 0, QApplication::UnicodeUTF8));
    rdExistingClassification->setText(QApplication::translate("fRecurrenceDialog", "Use existing model", 0, QApplication::UnicodeUTF8));
    rdCreateModel->setText(QApplication::translate("fRecurrenceDialog", "Train new model", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fRecurrenceDialog", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fRecurrenceDialog", "Cancel", 0, QApplication::UnicodeUTF8));

    cbT1Data->setText(QApplication::translate("fRecurrenceDialog", "Conventional (T1, T2, Flair, T1Gd)", 0, QApplication::UnicodeUTF8));
    cbDTIData->setText(QApplication::translate("fRecurrenceDialog", "DTI", 0, QApplication::UnicodeUTF8));
    cbPerfData->setText(QApplication::translate("fRecurrenceDialog", "Perfusion", 0, QApplication::UnicodeUTF8));
    cbDistanceData->setText(QApplication::translate("fRecurrenceDialog", "Distance", 0, QApplication::UnicodeUTF8));

    modelDirectoryLabel->setText(QApplication::translate("fRecurrenceDialog", "Model Directory:", 0, QApplication::UnicodeUTF8));
    trainingDirectoryLabel->setText(QApplication::translate("fRecurrenceDialog", "Training Directory:", 0, QApplication::UnicodeUTF8));
    testDirectoryLabel->setText(QApplication::translate("fRecurrenceDialog", "Test Directory:", 0, QApplication::UnicodeUTF8));
    outputDirectoryLabel->setText(QApplication::translate("fRecurrenceDialog", "Output Directory:", 0, QApplication::UnicodeUTF8));


  } // retranslateUi

};

namespace Ui {
  class fRecurrenceDialog : public ui_fRecurrenceDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fRecurrenceDialog_H






