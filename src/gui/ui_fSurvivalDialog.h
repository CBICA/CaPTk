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
  QPushButton * disclaimerButton;
  QLabel     *disclaimerLabel;

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
    fSurvivalPredictor->resize(200, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fSurvivalPredictor->sizePolicy().hasHeightForWidth());
    fSurvivalPredictor->setSizePolicy(sizePolicy);
    fSurvivalPredictor->setMinimumSize(QSize(0, 0));

    //fSurvivalPredictor->setModal(true);
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
    //svmModelFileName->setText("../data/survival/");

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

    QPalette* palette1 = new QPalette();
    palette1->setColor(QPalette::ButtonText, Qt::blue);
    disclaimerButton->setPalette(*palette1);

    QFont font("Bavaria");
    font.setPointSize(8);
    font.setWeight(QFont::Bold);
    font.setUnderline(TRUE);
    disclaimerButton->setFont(font);
    disclaimerButton->setStyleSheet("Text-align:left");

    //    classificationGridLayout->addWidget(rdCurrentSubject, 1, 0, 1, 1);
    //line_3 = new QFrame(classificationGroupBox);
    //line_3->setFrameStyle(QFrame::HLine);
    //line_3->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Expanding);
    //line_3->setFrameShadow(QFrame::Sunken);

    //classificationGridLayout->addWidget(line_3, 3, 0, 1, 4);

    classificationGridLayout->addWidget(rdExistingClassification, 1, 0, 1, 6);
    classificationGridLayout->addWidget(lModelPath, 2, 0, 1, 1);
    classificationGridLayout->addWidget(svmModelFileName, 2, 1, 1, 4);
    classificationGridLayout->addWidget(svmModelButton, 2, 5, 1, 1);

    classificationGridLayout->addWidget(lTestSubjectsPath, 3, 0, 1, 1);
    classificationGridLayout->addWidget(testSubjectsDirectoryName, 3, 1, 1, 4);
    classificationGridLayout->addWidget(testSubjectsDirectoryButton, 3, 5, 1, 1);


    classificationGridLayout->addWidget(disclaimerLabel, 4, 1, 1, 4);
    classificationGridLayout->addWidget(disclaimerButton, 4, 5, 1, 1);

    auto tempLabel = new QLabel;
    tempLabel->setText(QApplication::translate("fSurvivalPredictor", "Selected Subjects:", 0, QApplication::UnicodeUTF8));

    classificationGridLayout->addWidget(rdCreateModel, 5, 0, 1, 6);
    classificationGridLayout->addWidget(tempLabel, 6, 0, 1, 1);
    classificationGridLayout->addWidget(existingMaskDirectoryName, 6, 1, 1, 4);
    classificationGridLayout->addWidget(existingMasksButton, 6, 5, 1, 1);

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
    outputGridLayout->addWidget(outputDirectoryName, 0, 1, 1, 4);
    outputGridLayout->addWidget(outputDirectoryButton, 0, 5, 1, 1);

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
    rdExistingClassification->setText(QApplication::translate("fSurvivalPredictor", "Use an existing survival prediction model", 0, QApplication::UnicodeUTF8));
    lModelPath->setText(QApplication::translate("fSurvivalPredictor", "Model Directory:", 0, QApplication::UnicodeUTF8));
    lTestSubjectsPath->setText(QApplication::translate("fSurvivalPredictor", "Test Subjects:", 0, QApplication::UnicodeUTF8));
    outputDirectoryLabel->setText(QApplication::translate("fSurvivalPredictor", "Output Folder:", 0, QApplication::UnicodeUTF8));
    lTrainingNewSubject->setText(QApplication::translate("fSurvivalPredictor", "Training Subjects:", 0, QApplication::UnicodeUTF8));

    rdCreateModel->setText(QApplication::translate("fSurvivalPredictor", "Train a new survival prediction model", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fSurvivalPredictor", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fSurvivalPredictor", "Cancel", 0, QApplication::UnicodeUTF8));


  } // retranslateUi

};

namespace Ui {
  class fSurvivalPredictor : public ui_fSurvivalPredictor {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fSurvivalPredictor_H






