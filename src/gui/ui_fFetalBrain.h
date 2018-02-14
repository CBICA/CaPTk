#ifndef ui_fFetalBrain_H
#define ui_fFetalBrain_H

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


class ui_fFetalBrain
{
  public:
  
    QGridLayout *gridLayout;
    QVBoxLayout * verticalLayout;
  
    QGroupBox *inputGroupBox;
    QGridLayout * inputGridLayout;
    QPushButton *Segment;
    QPushButton *Predict;
    QGridLayout *outputgridlayout;

    QGroupBox *classificationGroupBox;
    QGridLayout * classificationGridLayout;
    QGroupBox *Confirmlabelbox;

    QRadioButton * rdExistingClassification;
    QLabel	* SubjLabel;
    QLineEdit * svmModelFileName;
    QPushButton * svmModelButton;

    QLineEdit * testSubjectsDirectoryName;
    QPushButton * testSubjectsButton;

    QRadioButton *rdCreateModel;
    QLineEdit * trainDirectoryName;
    QPushButton * TrainingSubDirButton;
    QLabel	* testSubjLabel;
    QLabel	* modelDirectoryLabel;
    QLabel	* trainLabel;
    QLabel	* trainDirectoryLabel;

    
 //   QGroupBox *outputGroupBox;
  //  QGridLayout *outputGridLayout;
    QLabel	*outputDirectoryLabel;
    QLineEdit *outputDirectoryName;
    QPushButton *outputDirectoryButton;
  
    QPushButton * confirmButton;
    QPushButton * cancelButton;

  void setupUi(QDialog *fFetalBrain)
  {

    if (fFetalBrain->objectName().isEmpty())
      fFetalBrain->setObjectName(QString::fromUtf8("fFetalBrain"));
 
    fFetalBrain->resize(200, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fFetalBrain->sizePolicy().hasHeightForWidth());
    fFetalBrain->setSizePolicy(sizePolicy);
    fFetalBrain->setMinimumSize(QSize(0, 0));

    fFetalBrain->setModal(false);



   //--------------------------input-------------------------------------------------------//

    gridLayout = new QGridLayout(fFetalBrain);
    gridLayout->setObjectName(QString::fromUtf8("inputgridLayout"));

    inputGroupBox = new QGroupBox(fFetalBrain);
    inputGridLayout = new QGridLayout(inputGroupBox);
    inputGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));
    rdExistingClassification = new QRadioButton(inputGroupBox);
    rdExistingClassification->setObjectName(QString::fromUtf8("existingClassification"));
    rdExistingClassification->setText(QApplication::translate("fFetalBrain", "Single Subject", 0, QApplication::UnicodeUTF8));

    Segment = new QPushButton(inputGroupBox);
    Segment->setObjectName(QString::fromUtf8("Segment"));
    Segment->setText(QString("Segment"));

    Predict = new QPushButton(inputGroupBox);
    Predict->setObjectName(QString::fromUtf8("Predict"));
    Predict->setText(QString("Predict"));

    inputGridLayout->addWidget(rdExistingClassification, 0, 0, 1, 1);    
    inputGridLayout->addWidget(Segment, 1, 0, 1, 1);
    inputGridLayout->addWidget(Predict, 1, 1, 1, 1);
    


     //--------------------------Classification-----------------------------------------------//
    classificationGroupBox = new QGroupBox(fFetalBrain);
    classificationGridLayout = new QGridLayout(classificationGroupBox);
    classificationGridLayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));
    QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);


  
    rdCreateModel = new QRadioButton(classificationGroupBox);
    rdCreateModel->setObjectName(QString::fromUtf8("createModel"));
    rdCreateModel->setText(QApplication::translate("fFetalBrain", "Train Subject", 0, QApplication::UnicodeUTF8));

    trainDirectoryLabel = new QLabel(classificationGroupBox);
    trainDirectoryLabel->setText(QString("Input Dir:"));


    trainDirectoryName = new QLineEdit(classificationGroupBox);
    trainDirectoryName->setObjectName(QString::fromUtf8("DirectoryName"));
    sizePolicy13.setHeightForWidth(trainDirectoryName->sizePolicy().hasHeightForWidth());
    trainDirectoryName->setSizePolicy(sizePolicy13);
  //  trainDirectoryName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);

    TrainingSubDirButton = new QPushButton(classificationGroupBox);
    TrainingSubDirButton->setObjectName(QString::fromUtf8("TrainingSubDirButton"));
    TrainingSubDirButton->setText(QString("Browse"));


    outputDirectoryLabel = new QLabel(classificationGroupBox);
    sizePolicy13.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    outputDirectoryLabel->setSizePolicy(sizePolicy13);
   // outputDirectoryLabel->setAlignment(Qt::AlignRight);
    outputDirectoryLabel->setText("Output Dir:");

    outputDirectoryName = new QLineEdit(".\\");
    outputDirectoryName->setObjectName(QString::fromUtf8("outputDirectoryName"));
    sizePolicy13.setHeightForWidth(outputDirectoryName->sizePolicy().hasHeightForWidth());
    outputDirectoryName->setSizePolicy(sizePolicy13);
  //  outputDirectoryName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);

    outputDirectoryButton = new QPushButton(classificationGroupBox);
    outputDirectoryButton->setObjectName(QString::fromUtf8("outputDirectoryButton"));
    outputDirectoryButton->setText(QString("Browse"));


    classificationGridLayout->addWidget(rdCreateModel, 0, 0, 1, 1);
    classificationGridLayout->addWidget(trainDirectoryLabel, 1, 0, 1, 1);
    classificationGridLayout->addWidget(trainDirectoryName, 1,1, 1, 1);
    classificationGridLayout->addWidget(TrainingSubDirButton,1, 2, 1, 1);

    classificationGridLayout->addWidget(outputDirectoryLabel, 2, 0, 1, 1);
    classificationGridLayout->addWidget(outputDirectoryName, 2, 1, 1, 1);
    classificationGridLayout->addWidget(outputDirectoryButton, 2, 2, 1, 1);


    //--------------------------Classification-----------------------------------------------//
    Confirmlabelbox = new QGroupBox(fFetalBrain);
    outputgridlayout = new QGridLayout(Confirmlabelbox);
    outputgridlayout->setObjectName(QString::fromUtf8("imagestabgridLayout3"));


    confirmButton = new QPushButton(Confirmlabelbox);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    confirmButton->setText(QString("Confirm"));

    cancelButton = new QPushButton(Confirmlabelbox);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    cancelButton->setText(QString("Cancel"));


    outputgridlayout->addWidget(confirmButton, 0,0,1,1);
    outputgridlayout->addWidget(cancelButton, 0, 1, 1, 1);

    gridLayout->addWidget(inputGroupBox, 0, 0, 1, 1);
    gridLayout->addWidget(classificationGroupBox, 1, 0, 1, 1);
    gridLayout->addWidget(Confirmlabelbox, 2, 0, 1, 1);

    QMetaObject::connectSlotsByName(fFetalBrain);
  } // setupUi

};

namespace Ui {
  class fFetalBrain : public ui_fFetalBrain{};
} // namespace Ui


QT_END_NAMESPACE

#endif // ui_fFetalBrain_H






