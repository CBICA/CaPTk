#ifndef ui_fPerfusionEstimator_H
#define ui_fPerfusionEstimator_H

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

class ui_fPerfusionEstimator
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;


  QGroupBox *inputGroupBox;
  QGridLayout *inputGridLayout;
  QLabel *inputImageLabel;
  QLineEdit *inputImageName;
  QPushButton *inputImageButton;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputImageName;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QCheckBox* m_psr;
  QCheckBox* m_rcbv;
  QCheckBox* m_ph;

  QHBoxLayout * horizontalLayout;

  void setupUi(QDialog *fPerfusionEstimator)
  {
    if (fPerfusionEstimator->objectName().isEmpty())
      fPerfusionEstimator->setObjectName(QString::fromUtf8("fPerfusionEstimator"));
    //fPerfusionEstimator->setWindowModality(Qt::NonModal);
    fPerfusionEstimator->resize(400, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPerfusionEstimator->sizePolicy().hasHeightForWidth());
    fPerfusionEstimator->setSizePolicy(sizePolicy);
    fPerfusionEstimator->setMinimumSize(QSize(0, 0));

    fPerfusionEstimator->setModal(true);
    gridLayout = new QGridLayout(fPerfusionEstimator);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    //input 
    inputGroupBox = new QGroupBox(fPerfusionEstimator);
    inputGroupBox->setTitle(QString::fromStdString("Input Data"));

    inputGridLayout = new QGridLayout(inputGroupBox);
    inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    //inputImageLabel = new QLabel(inputGroupBox);
    //sizePolicy.setHeightForWidth(inputImageLabel->sizePolicy().hasHeightForWidth());
    //inputImageLabel->setSizePolicy(sizePolicy);

    inputImageName = new QLineEdit(" ");
    inputImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputImageName->sizePolicy().hasHeightForWidth());
    inputImageName->setSizePolicy(sizePolicy);
    inputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputImageButton = new QPushButton(inputGroupBox);
    inputImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputImageButton->setText(QString("Browse"));

    m_psr = new QCheckBox("Percent Signal Recovery");
    m_psr->setEnabled(true);
    m_rcbv = new QCheckBox("Reletive Cerebal Blood Volume");
    m_rcbv->setEnabled(true);
    m_ph = new QCheckBox("Peak Height");
    m_ph->setEnabled(true);
    
    //inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 0, 1, 1, 1);

    inputGridLayout->addWidget(m_rcbv, 2, 0, 1, 1);
    inputGridLayout->addWidget(m_ph, 3, 0, 1, 1);
    inputGridLayout->addWidget(m_psr, 4, 0, 1, 1);

    // output 
    outputGroupBox = new QGroupBox(fPerfusionEstimator);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputImageLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    //outputImageLabel->setSizePolicy(sizePolicy);

    outputImageName = new QLineEdit(" ");
    outputImageName->setObjectName(QString::fromUtf8("outputImageName"));
    sizePolicy.setHeightForWidth(outputImageName->sizePolicy().hasHeightForWidth());
    outputImageName->setSizePolicy(sizePolicy);
    outputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputImageButton = new QPushButton(outputGroupBox);
    outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
    //outputImageButton->setIcon(ButtonIcon);
    //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    outputImageButton->setText(QString("Browse"));

    //outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageName, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);

    // put the layout in perspective
    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fPerfusionEstimator);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fPerfusionEstimator);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fPerfusionEstimator);

    QMetaObject::connectSlotsByName(fPerfusionEstimator);
  } // setupUi

  void retranslateUi(QDialog *fPerfusionEstimator)
  {
    fPerfusionEstimator->setWindowTitle(QApplication::translate("fPerfusionEstimator", "Perfusion Derivatives", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fPerfusionEstimator", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fPerfusionEstimator", "Cancel", 0, QApplication::UnicodeUTF8));
    m_rcbv->setText(QApplication::translate("fPerfusionEstimator", "Relative Cerebral Blood Volume", 0, QApplication::UnicodeUTF8));
    m_psr->setText(QApplication::translate("fPerfusionEstimator", "Percent Signal Recovery", 0, QApplication::UnicodeUTF8));
    m_ph->setText(QApplication::translate("fPerfusionEstimator", "Peak Height", 0, QApplication::UnicodeUTF8));
  } // retranslateUi
};

namespace Ui {
  class fPerfusionEstimator : public ui_fPerfusionEstimator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPerfusionEstimator_H