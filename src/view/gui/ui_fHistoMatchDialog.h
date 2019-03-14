#ifndef ui_fHistoMatch_H
#define ui_fHistoMatch_H

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

QT_BEGIN_NAMESPACE

class ui_fHistoMatcher
{
public:
  QGridLayout *gridLayout;
  QVBoxLayout * verticalLayout;

  QGroupBox *referenceGroupBox;
  QGridLayout *referenceGridLayout;
  QLabel *referenceImageLabel;
  QLineEdit *referenceImageName;
  QPushButton *referenceImageButton;

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

  QFrame *line_3;

  QHBoxLayout * horizontalLayout;

  //QLabel *longRunningWarning;

  void setupUi(QDialog *fHistoMatcher)
  {

    if (fHistoMatcher->objectName().isEmpty())
      fHistoMatcher->setObjectName(QString::fromUtf8("fHistoMatcher"));
    //fHistoMatcher->setWindowModality(Qt::NonModal);
    fHistoMatcher->resize(200, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fHistoMatcher->sizePolicy().hasHeightForWidth());
    fHistoMatcher->setSizePolicy(sizePolicy);
    fHistoMatcher->setMinimumSize(QSize(0, 0));

    //fHistoMatcher->setModal(true);
    gridLayout = new QGridLayout(fHistoMatcher);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    // reference
    referenceGroupBox = new QGroupBox(fHistoMatcher);
    referenceGroupBox->setTitle(QString::fromStdString("Reference Image"));

    referenceGridLayout = new QGridLayout(referenceGroupBox);
    referenceGridLayout->setObjectName(QString::fromUtf8("referenceGridLayout"));

    referenceImageLabel = new QLabel(referenceGroupBox);
    sizePolicy.setHeightForWidth(referenceImageLabel->sizePolicy().hasHeightForWidth());
    referenceImageLabel->setSizePolicy(sizePolicy);

    referenceImageName = new QLineEdit("");
    referenceImageName->setObjectName(QString::fromUtf8("referenceImageName"));
    sizePolicy.setHeightForWidth(referenceImageName->sizePolicy().hasHeightForWidth());
    referenceImageName->setSizePolicy(sizePolicy);
    referenceImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    referenceImageButton = new QPushButton(referenceGroupBox);
    referenceImageButton->setObjectName(QString::fromUtf8("referenceImageButton"));
    referenceImageButton->setText(QString("Browse"));
    //referenceMaskButton->setToolTip(QString("Atlas based on which stripping is done"));

    referenceGridLayout->addWidget(referenceImageLabel, 0, 0, 1, 1);
    referenceGridLayout->addWidget(referenceImageName, 1, 0, 1, 1);
    referenceGridLayout->addWidget(referenceImageButton, 1, 1, 1, 1);

    // input
    inputGroupBox = new QGroupBox(fHistoMatcher);
    inputGroupBox->setTitle(QString::fromStdString("Input Image"));

    inputGridLayout = new QGridLayout(inputGroupBox);
    inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

    inputImageLabel = new QLabel(inputGroupBox);
    sizePolicy.setHeightForWidth(inputImageLabel->sizePolicy().hasHeightForWidth());
    inputImageLabel->setSizePolicy(sizePolicy);

    inputImageName = new QLineEdit("");
    inputImageName->setObjectName(QString::fromUtf8("inputImageName"));
    sizePolicy.setHeightForWidth(inputImageName->sizePolicy().hasHeightForWidth());
    inputImageName->setSizePolicy(sizePolicy);
    inputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    inputImageButton = new QPushButton(inputGroupBox);
    inputImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
    inputImageButton->setText(QString("Browse"));

    inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
    inputGridLayout->addWidget(inputImageName, 1, 0, 1, 1);
    inputGridLayout->addWidget(inputImageButton, 1, 1, 1, 1);

    // output
    outputGroupBox = new QGroupBox(fHistoMatcher);
    outputGroupBox->setTitle(QString::fromStdString("Output"));

    outputGridLayout = new QGridLayout(outputGroupBox);
    outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

    outputImageLabel = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
    outputImageLabel->setSizePolicy(sizePolicy);

    outputImageName = new QLineEdit("");
    outputImageName->setObjectName(QString::fromUtf8("outputImageName"));
    sizePolicy.setHeightForWidth(outputImageName->sizePolicy().hasHeightForWidth());
    outputImageName->setSizePolicy(sizePolicy);
    outputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputImageButton = new QPushButton(outputGroupBox);
    outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
    //outputImageButton->setIcon(ButtonIcon);
    //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent
    outputImageButton->setText(QString("Browse"));

    outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
    outputGridLayout->addWidget(outputImageName, 1, 0, 1, 1);
    outputGridLayout->addWidget(outputImageButton, 1, 1, 1, 1);

    // put the layout in perspective
    gridLayout->addWidget(referenceGroupBox, 0, 0, 1, 2);
    gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
    gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


    confirmButton = new QPushButton(fHistoMatcher);
    confirmButton->setObjectName(QString::fromUtf8("confirm"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    cancelButton = new QPushButton(fHistoMatcher);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

    gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

    retranslateUi(fHistoMatcher);

    QMetaObject::connectSlotsByName(fHistoMatcher);
  } // setupUi

  void retranslateUi(QDialog *fHistoMatcher)
  {
    fHistoMatcher->setWindowTitle(QApplication::translate("fHistoMatcher", "Histogram Matching", 0));
    confirmButton->setText(QApplication::translate("fHistoMatcher", "Confirm", 0));
    cancelButton->setText(QApplication::translate("fHistoMatcher", "Cancel", 0));

  } // retranslateUi

};

namespace Ui {
  class fHistoMatcher : public ui_fHistoMatcher {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fHistoMatcher_H
