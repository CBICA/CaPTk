#ifndef ui_fDirectionalityDialog_H
#define ui_fDirectionalityDialog_H

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

class ui_fDirectionalityDialog
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *roi1GroupBox;
  QGridLayout *roi1GridLayout;
  QLabel *roi1ImageLabel;
  QLineEdit *roi1ImageName;
  QPushButton *roi1ImageButton;

  QGroupBox *roi2GroupBox;
  QGridLayout *roi2GridLayout;
  QLabel *roi2ImageLabel;
  QLineEdit *roi2ImageName;
  QPushButton *roi2ImageButton;

  QGroupBox *outputGroupBox;
  QGridLayout *outputGridLayout;
  QLineEdit *outputDirectory;
  QLabel *outputImageLabel;
  QPushButton *outputImageButton;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  void setupUi(QDialog *fDirectionalityDialog)
  {

    if (fDirectionalityDialog->objectName().isEmpty())
      fDirectionalityDialog->setObjectName(QString::fromUtf8("fDirectionalityDialog"));
    //fDirectionalityDialog->setWindowModality(Qt::NonModal);
    fDirectionalityDialog->resize(200, 200); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fDirectionalityDialog->sizePolicy().hasHeightForWidth());
    fDirectionalityDialog->setSizePolicy(sizePolicy);
    fDirectionalityDialog->setMinimumSize(QSize(0, 0));

    fDirectionalityDialog->setModal(true);
    gridLayout = new QGridLayout(fDirectionalityDialog);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // ROI-1 
    roi1GroupBox = new QGroupBox(fDirectionalityDialog);
    roi1GroupBox->setTitle(QString::fromStdString("ROI Pre-Injection"));

    roi1GridLayout = new QGridLayout(roi1GroupBox);
    roi1GridLayout->setObjectName(QString::fromUtf8("roi1GridLayout"));

    //roi1ImageLabel = new QLabel(roi1GroupBox);
    //sizePolicy.setHeightForWidth(roi1ImageLabel->sizePolicy().hasHeightForWidth());
    //roi1ImageLabel->setSizePolicy(sizePolicy);

    roi1ImageName = new QLineEdit(" ");
    roi1ImageName->setObjectName(QString::fromUtf8("roi1ImageName"));
    sizePolicy.setHeightForWidth(roi1ImageName->sizePolicy().hasHeightForWidth());
    roi1ImageName->setSizePolicy(sizePolicy);
    roi1ImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    roi1ImageButton = new QPushButton(roi1GroupBox);
    roi1ImageButton->setObjectName(QString::fromUtf8("roi1ImageButton"));
    roi1ImageButton->setText(QString("File Name"));

    //roi1GridLayout->addWidget(roi1ImageLabel, 0, 0, 1, 1);
    roi1GridLayout->addWidget(roi1ImageName, 0, 0, 1, 1);
    roi1GridLayout->addWidget(roi1ImageButton, 0, 1, 1, 1);

    //referenceMaskButton->setToolTip(QString("Atlas based on which stripping is done"));

    // ROI-2
    roi2GroupBox = new QGroupBox(fDirectionalityDialog);
    roi2GroupBox->setTitle(QString::fromStdString("ROI Post-Injection)"));

    roi2GridLayout = new QGridLayout(roi2GroupBox);
    roi2GridLayout->setObjectName(QString::fromUtf8("roi2GridLayout"));

    //roi2ImageLabel = new QLabel(roi2GroupBox);
    //sizePolicy.setHeightForWidth(roi2ImageLabel->sizePolicy().hasHeightForWidth());
    //roi2ImageLabel->setSizePolicy(sizePolicy);

    roi2ImageName = new QLineEdit(" ");
    roi2ImageName->setObjectName(QString::fromUtf8("roi2ImageName"));
    sizePolicy.setHeightForWidth(roi2ImageName->sizePolicy().hasHeightForWidth());
    roi2ImageName->setSizePolicy(sizePolicy);
    roi2ImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    roi2ImageButton = new QPushButton(roi2GroupBox);
    roi2ImageButton->setObjectName(QString::fromUtf8("roi2ImageButton"));
    roi2ImageButton->setText(QString("File Name"));
    //referenceMaskButton->setToolTip(QString("Atlas based on which stripping is done"));

    //roi2GridLayout->addWidget(roi2ImageLabel, 0, 0, 1, 1);
    roi2GridLayout->addWidget(roi2ImageName, 0, 0, 1, 1);
    roi2GridLayout->addWidget(roi2ImageButton, 0, 1, 1, 1);

    // output 
    outputGroupBox = new QGroupBox(fDirectionalityDialog);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    //outputImageLabel = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    //outputImageLabel->setSizePolicy(sizePolicy);
    //outputImageLabel->setText("Select Folder");

    outputDirectory = new QLineEdit(" ");
    outputDirectory->setObjectName(QString::fromUtf8("outputDirectory"));
    sizePolicy.setHeightForWidth(outputDirectory->sizePolicy().hasHeightForWidth());
    outputDirectory->setSizePolicy(sizePolicy);
    outputDirectory->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputImageButton = new QPushButton(outputGroupBox);
    outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
    //outputImageButton->setIcon(ButtonIcon);
    //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
    outputImageButton->setText(QString("Select Folder"));

    //outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputDirectory, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);

    // put the layout in perspective
    gridLayout->addWidget(roi1GroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(roi2GroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fDirectionalityDialog);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fDirectionalityDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fDirectionalityDialog);

    QMetaObject::connectSlotsByName(fDirectionalityDialog);
  } // setupUi

  void retranslateUi(QDialog *fDirectionalityDialog)
  {
    fDirectionalityDialog->setWindowTitle(QApplication::translate("fDirectionalityDialog", "Directionality Estimator", 0, QApplication::UnicodeUTF8));
    confirmButton->setText(QApplication::translate("fDirectionalityDialog", "Confirm", 0, QApplication::UnicodeUTF8));
    cancelButton->setText(QApplication::translate("fDirectionalityDialog", "Cancel", 0, QApplication::UnicodeUTF8));


  } // retranslateUi

};

namespace Ui {
  class fDirectionalityDialog : public ui_fDirectionalityDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDirectionalityDialog_H






