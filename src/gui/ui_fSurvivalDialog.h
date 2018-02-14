#ifndef ui_fSurvivalPredictor_H
#define ui_fSurvivalPredictor_H

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
#include <QtGui/QLabel>

QT_BEGIN_NAMESPACE

class ui_fSurvivalPredictor
{
public:
  QGridLayout *gridLayout_3;
  QVBoxLayout * verticalLayout;
  QGroupBox *classificationGroupBox;
  QGridLayout * classificationGridLayout;
  QGroupBox *outputGroupBox;
  QRadioButton * rdExistingClassification;
  QRadioButton *rdCreateModel;
  QLineEdit * svmModelFileName;
  QPushButton * svmModelButton;
  QLineEdit * testSubjectsDirectoryName;
  QPushButton * testSubjectsDirectoryButton;

  QLabel *lModelPath = new QLabel;
  QLabel *lTestSubjectsPath = new QLabel;
  QLabel *lTrainingSubjectsPath = new QLabel;
  QLabel *lTrainingNewSubject = new QLabel;


  QLineEdit * existingMaskDirectoryName;
  QPushButton * existingMasksButton;

  QGridLayout *outputGridLayout;
  QLabel	*outputDirectoryLabel;
  QLineEdit *outputDirectoryName;
  QPushButton *outputDirectoryButton;

  QPushButton * confirmButton;
  QPushButton * cancelButton;

  //QFrame *line_3;
  //QFrame *line_4;
  //QFrame *line_5;
  //QFrame *line_6;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fSurvivalPredictor)
  {
    if (fSurvivalPredictor->objectName().isEmpty())
      fSurvivalPredictor->setObjectName(QString::fromUtf8("fSurvivalPredictor"));
    fSurvivalPredictor->setWindowModality(Qt::ApplicationModal);
    fSurvivalPredictor->resize(200,200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fSurvivalPredictor->sizePolicy().hasHeightForWidth());
    fSurvivalPredictor->setSizePolicy(sizePolicy);
    fSurvivalPredictor->setMinimumSize(QSize(0, 0));

    fSurvivalPredictor->setModal(true);
    gridLayout_3 = new QGridLayout(fSurvivalPredictor);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //--------------------------------------------------------------------
    classificationGroupBox = new QGroupBox(fSurvivalPredictor);
    classificationGroupBox->setTitle(QString::fromStdString("Classification"));
    classificationGridLayout = new QGridLayout(classificationGroupBox);
    classificationGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));

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
    svmModelFileName->setText("../data/survival/");

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
    testSubjectsDirectoryButton->setToolTip(QString("Dir containing Test subjects"));


    existingMasksButton = new QPushButton(classificationGroupBox);
    existingMasksButton->setObjectName(QString::fromUtf8("existingMasksButton"));
    //existingMasksButton->setIcon(ButtonIcon);
    //existingMasksButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    existingMasksButton->setText(QString("Browse"));
    existingMasksButton->setToolTip(QString("Dir containing Training subjects"));


//    classificationGridLayout->addWidget(rdCurrentSubject, 1, 0, 1, 1);
    //line_3 = new QFrame(classificationGroupBox);
    //line_3->setFrameStyle(QFrame::HLine);
    //line_3->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    //line_3->setFrameShadow(QFrame::Sunken);

    //classificationGridLayout->addWidget(line_3, 3, 0, 1, 4);

    classificationGridLayout->addWidget(rdExistingClassification, 5, 0, 1, 1);
    classificationGridLayout->addWidget(lModelPath, 6, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName, 6, 1, 1, 1);
    classificationGridLayout->addWidget(svmModelButton, 6, 2, 1, 1);

    //line_4 = new QFrame(classificationGroupBox);
    //line_4->setFrameStyle(QFrame::HLine);
    //line_4->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    //line_4->setFrameShadow(QFrame::Sunken);

    //classificationGridLayout->addWidget(line_4, 7, 1, 2, 3);

    classificationGridLayout->addWidget(lTestSubjectsPath, 8, 0, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryName, 8, 1, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryButton, 8, 2, 1, 1);

    //line_5 = new QFrame(classificationGroupBox);
    //line_5->setFrameStyle(QFrame::HLine);
    //line_5->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    //line_5->setFrameShadow(QFrame::Sunken);

    //classificationGridLayout->addWidget(line_5, 9, 1, 2, 2);

    auto tempLabel = new QLabel;
    tempLabel->setText(QApplication::translate("fSurvivalPredictor", "Select Subjects:", 0, QApplication::UnicodeUTF8));

    classificationGridLayout->addWidget(rdCreateModel, 11, 0, 1, 1);
    classificationGridLayout->addWidget(tempLabel, 12, 0, 1, 1);
    classificationGridLayout->addWidget(existingMaskDirectoryName, 12, 1, 1, 1);
    classificationGridLayout->addWidget(existingMasksButton, 12, 2, 1, 1);

    //line_6 = new QFrame(classificationGroupBox);
    //line_6->setFrameStyle(QFrame::HLine);
    //line_6->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    //line_6->setFrameShadow(QFrame::Sunken);

    //classificationGridLayout->addWidget(line_6, 14, 1, 2, 2);


    //--------------------------output-------------------------------------------------------


    outputGroupBox = new QGroupBox(fSurvivalPredictor);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputDirectoryLabel = new QLabel(outputGroupBox);
    sizePolicy13.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    outputDirectoryLabel->setSizePolicy(sizePolicy13);

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
    outputGridLayout->addWidget(outputDirectoryName, 0, 1, 1, 1);
    outputGridLayout->addWidget(outputDirectoryButton, 0, 2, 1, 1);

    gridLayout_3->addWidget(classificationGroupBox, 0, 0, 1, 2);
    gridLayout_3->addWidget(outputGroupBox, 1, 0, 1, 2);


    confirmButton = new QPushButton(fSurvivalPredictor);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fSurvivalPredictor);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout_3->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout_3->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fSurvivalPredictor);

    QMetaObject::connectSlotsByName(fSurvivalPredictor);
  } // setupUi

  void retranslateUi(QDialog *fSurvivalPredictor)
  {
    fSurvivalPredictor->setWindowTitle(QApplication::translate("fSurvivalPredictor", "Survival Prediction Index", 0, QApplication::UnicodeUTF8));
//    rdCurrentSubject->setText(QApplication::translate("fSurvivalPredictor", "Loaded subject", 0, QApplication::UnicodeUTF8));
    rdExistingClassification->setText(QApplication::translate("fSurvivalPredictor", "Use existing model", 0, QApplication::UnicodeUTF8));
    lModelPath->setText(QApplication::translate("fSurvivalPredictor", "Model Directory:", 0, QApplication::UnicodeUTF8));
    lTestSubjectsPath->setText(QApplication::translate("fSurvivalPredictor", "Test Subjects:", 0, QApplication::UnicodeUTF8));
    outputDirectoryLabel->setText(QApplication::translate("fSurvivalPredictor", "Folder to Save Results:", 0, QApplication::UnicodeUTF8));
    lTrainingNewSubject->setText(QApplication::translate("fSurvivalPredictor", "Training Subjects:", 0, QApplication::UnicodeUTF8));

    rdCreateModel->setText(QApplication::translate("fSurvivalPredictor", "Train New Model", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fSurvivalPredictor", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fSurvivalPredictor", "Cancel", 0, QApplication::UnicodeUTF8));


  } // retranslateUi

};

namespace Ui {
  class fSurvivalPredictor : public ui_fSurvivalPredictor {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fSurvivalPredictor_H






