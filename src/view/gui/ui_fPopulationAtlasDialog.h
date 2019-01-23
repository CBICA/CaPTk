#ifndef ui_fPopulationAtlasDialog_H
#define ui_fPopulationAtlasDialog_H

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

class ui_fPopulationAtlasDialog
{
public:
  QGridLayout *gridLayout_3;
  QLineEdit * inputdirectoryName;
  QLineEdit * inputlabelName;
  QLineEdit * inputAtlasName;
  QLineEdit * outputdirectoryName;

  QPushButton * inputdirectoryButton;
  QPushButton * inputlabelButton;
  QPushButton * inputAtlasButton;
  QPushButton * outputdirectoryButton;


  QPushButton * confirmButton;
  QPushButton * cancelButton;

  QLabel	*outputDirectoryLabel;
  QLabel	*inputDirectoryLabel;
  QLabel	*inputlabelLabel;
  QLabel	*inputAtlasLabel;


  void setupUi(QDialog *fPopulationAtlasDialog)
  {

    if (fPopulationAtlasDialog->objectName().isEmpty())
      fPopulationAtlasDialog->setObjectName(QString::fromUtf8("fPopulationAtlasDialog"));
    fPopulationAtlasDialog->setWindowModality(Qt::ApplicationModal);
    fPopulationAtlasDialog->resize(400, 200); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fPopulationAtlasDialog->sizePolicy().hasHeightForWidth());
    fPopulationAtlasDialog->setSizePolicy(sizePolicy);
    fPopulationAtlasDialog->setMinimumSize(QSize(0, 0));

    //fPopulationAtlasDialog->setModal(true);
    gridLayout_3 = new QGridLayout(fPopulationAtlasDialog);
    gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));

    //--------------------------------------------------------------------
    inputdirectoryName = new QLineEdit(fPopulationAtlasDialog);
    inputdirectoryName->setObjectName(QString::fromUtf8("inputdirectoryName"));
    QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy13.setHorizontalStretch(0);
    sizePolicy13.setVerticalStretch(0);
    sizePolicy13.setHeightForWidth(inputdirectoryName->sizePolicy().hasHeightForWidth());
    inputdirectoryName->setSizePolicy(sizePolicy13);
    inputdirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

    outputdirectoryName = new QLineEdit(fPopulationAtlasDialog);
    outputdirectoryName->setObjectName(QString::fromUtf8("outputdirectoryName"));
    sizePolicy13.setHeightForWidth(outputdirectoryName->sizePolicy().hasHeightForWidth());
    outputdirectoryName->setSizePolicy(sizePolicy13);
    outputdirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    inputlabelName = new QLineEdit(fPopulationAtlasDialog);
    inputlabelName->setObjectName(QString::fromUtf8("inputlabelName"));
    sizePolicy13.setHeightForWidth(inputlabelName->sizePolicy().hasHeightForWidth());
    inputlabelName->setSizePolicy(sizePolicy13);
    inputlabelName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    inputAtlasName = new QLineEdit(fPopulationAtlasDialog);
    inputAtlasName->setObjectName(QString::fromUtf8("inputAtlasName"));
    sizePolicy13.setHeightForWidth(inputAtlasName->sizePolicy().hasHeightForWidth());
    inputAtlasName->setSizePolicy(sizePolicy13);
    inputAtlasName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);




    inputdirectoryButton = new QPushButton(fPopulationAtlasDialog);
    inputdirectoryButton->setObjectName(QString::fromUtf8("inputdirectoryButton"));
    inputdirectoryButton->setText(QString("Browse"));
    inputdirectoryButton->setToolTip(QString("Directory containing Input subjects"));

    outputdirectoryButton = new QPushButton(fPopulationAtlasDialog);
    outputdirectoryButton->setObjectName(QString::fromUtf8("outputdirectoryButton"));
    outputdirectoryButton->setText(QString("Browse"));
    outputdirectoryButton->setWhatsThis(QString("&nbsp;The meaning of the Source field depends "
      "on the Type field:"
      "<ul>"
      "<li><b>Books</b> have a Publisher"
      "<li><b>Articles</b> have a Journal name with "
      "volume and issue number"
      "<li><b>Theses</b> have an Institution name "
      "and a Department name"
      "</ul>"));

    outputDirectoryLabel = new QLabel(fPopulationAtlasDialog);
    sizePolicy13.setHeightForWidth(outputDirectoryLabel->sizePolicy().hasHeightForWidth());
    outputDirectoryLabel->setSizePolicy(sizePolicy13);

    inputDirectoryLabel = new QLabel(fPopulationAtlasDialog);
    sizePolicy13.setHeightForWidth(inputDirectoryLabel->sizePolicy().hasHeightForWidth());
    inputDirectoryLabel->setSizePolicy(sizePolicy13);

    inputlabelLabel = new QLabel(fPopulationAtlasDialog);
    sizePolicy13.setHeightForWidth(inputlabelLabel->sizePolicy().hasHeightForWidth());
    inputlabelLabel->setSizePolicy(sizePolicy13);

    inputAtlasLabel = new QLabel(fPopulationAtlasDialog);
    sizePolicy13.setHeightForWidth(inputAtlasLabel->sizePolicy().hasHeightForWidth());
    inputAtlasLabel->setSizePolicy(sizePolicy13);

    inputlabelButton = new QPushButton(fPopulationAtlasDialog);
    inputlabelButton->setObjectName(QString::fromUtf8("Confirm"));
    inputlabelButton->setText(QString("Browse"));

    inputAtlasButton = new QPushButton(fPopulationAtlasDialog);
    inputAtlasButton->setObjectName(QString::fromUtf8("Cancel"));
    inputAtlasButton->setText(QString("Browse"));

    confirmButton = new QPushButton(fPopulationAtlasDialog);
    confirmButton->setObjectName(QString::fromUtf8("Confirm"));
    confirmButton->setText(QString("Confirm"));

    cancelButton = new QPushButton(fPopulationAtlasDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    cancelButton->setText(QString("cancel"));



    gridLayout_3->addWidget(inputDirectoryLabel, 1, 0, 1, 1);
    gridLayout_3->addWidget(inputdirectoryName, 1, 1, 1, 4);
    gridLayout_3->addWidget(inputdirectoryButton, 1, 5, 1, 1);

    gridLayout_3->addWidget(inputlabelLabel, 2, 0, 1, 1);
    gridLayout_3->addWidget(inputlabelName, 2, 1, 1, 4);
    gridLayout_3->addWidget(inputlabelButton, 2, 5, 1, 1);

    gridLayout_3->addWidget(inputAtlasLabel, 3, 0, 1, 1);
    gridLayout_3->addWidget(inputAtlasName, 3, 1, 1, 4);
    gridLayout_3->addWidget(inputAtlasButton, 3, 5, 1, 1);

    //longRunningWarning = new QLabel(outputGroupBox);
    //sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    //longRunningWarning->setSizePolicy(sizePolicy);
    //longRunningWarning->setAlignment(Qt::AlignRight);
    //longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

    gridLayout_3->addWidget(outputDirectoryLabel, 4, 0, 1, 1);
    gridLayout_3->addWidget(outputdirectoryName, 4, 1, 1, 4);
    gridLayout_3->addWidget(outputdirectoryButton, 4, 5, 1, 1);


    gridLayout_3->addWidget(confirmButton, 5, 0, 1, 1);
    gridLayout_3->addWidget(cancelButton, 5, 1, 1, 1);

    retranslateUi(fPopulationAtlasDialog);

    QMetaObject::connectSlotsByName(fPopulationAtlasDialog);
  } // setupUi

  void retranslateUi(QDialog *fPopulationAtlasDialog)
  {
    //  outputDirectoryLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Output Directory:", 0, QApplication::UnicodeUTF8));
    //  inputDirectoryLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Input Directory:", 0, QApplication::UnicodeUTF8));
   // inputlabelLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Label file:", 0, QApplication::UnicodeUTF8));
   // inputAtlasLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Atlas file:", 0, QApplication::UnicodeUTF8));
    // NEW CHANGES
    outputDirectoryLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Output Directory:", 0));
    inputDirectoryLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Input Directory:", 0));
    inputlabelLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Label File:", 0));
    inputAtlasLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Atlas File:", 0));
  } // retranslateUi

};

namespace Ui {
  class fPopulationAtlasDialog : public ui_fPopulationAtlasDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPopulationAtlasDialog_H
