#ifndef ui_fPopulationAtlasDialog_H
#define ui_fPopulationAtlasDialog_H

#include <QtCore/QVariant>
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
  QLineEdit * inputfileName;
  QLineEdit * inputAtlasName;
  QLineEdit * outputdirectoryName;

  QPushButton * inputfileButton;
  QPushButton * inputAtlasButton;
  QPushButton * outputdirectoryButton;


  QPushButton * confirmButton;
  QPushButton * cancelButton;

  QLabel	*outputDirectoryLabel;
  QLabel	*inputfileLabel;
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
    QSizePolicy sizePolicy13(QSizePolicy::Preferred, QSizePolicy::Fixed);
    sizePolicy13.setHorizontalStretch(0);
    sizePolicy13.setVerticalStretch(0);

    outputdirectoryName = new QLineEdit(fPopulationAtlasDialog);
    outputdirectoryName->setObjectName(QString::fromUtf8("outputdirectoryName"));
    sizePolicy13.setHeightForWidth(outputdirectoryName->sizePolicy().hasHeightForWidth());
    outputdirectoryName->setSizePolicy(sizePolicy13);
    outputdirectoryName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    inputfileName = new QLineEdit(fPopulationAtlasDialog);
    inputfileName->setObjectName(QString::fromUtf8("inputfileName"));
    sizePolicy13.setHeightForWidth(inputfileName->sizePolicy().hasHeightForWidth());
    inputfileName->setSizePolicy(sizePolicy13);
    inputfileName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


    inputAtlasName = new QLineEdit(fPopulationAtlasDialog);
    inputAtlasName->setObjectName(QString::fromUtf8("inputAtlasName"));
    sizePolicy13.setHeightForWidth(inputAtlasName->sizePolicy().hasHeightForWidth());
    inputAtlasName->setSizePolicy(sizePolicy13);
    inputAtlasName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

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

    inputfileLabel = new QLabel(fPopulationAtlasDialog);
    sizePolicy13.setHeightForWidth(inputfileLabel->sizePolicy().hasHeightForWidth());
    inputfileLabel->setSizePolicy(sizePolicy13);

    inputAtlasLabel = new QLabel(fPopulationAtlasDialog);
    sizePolicy13.setHeightForWidth(inputAtlasLabel->sizePolicy().hasHeightForWidth());
    inputAtlasLabel->setSizePolicy(sizePolicy13);

    inputfileButton = new QPushButton(fPopulationAtlasDialog);
    inputfileButton->setObjectName(QString::fromUtf8("Confirm"));
    inputfileButton->setText(QString("Browse"));

    inputAtlasButton = new QPushButton(fPopulationAtlasDialog);
    inputAtlasButton->setObjectName(QString::fromUtf8("Cancel"));
    inputAtlasButton->setText(QString("Browse"));

    confirmButton = new QPushButton(fPopulationAtlasDialog);
    confirmButton->setObjectName(QString::fromUtf8("Confirm"));
    confirmButton->setText(QString("Confirm"));

    cancelButton = new QPushButton(fPopulationAtlasDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    cancelButton->setText(QString("cancel"));

    gridLayout_3->addWidget(inputfileLabel, 2, 0, 1, 1);
    gridLayout_3->addWidget(inputfileName, 2, 1, 1, 4);
    gridLayout_3->addWidget(inputfileButton, 2, 5, 1, 1);

    gridLayout_3->addWidget(inputAtlasLabel, 3, 0, 1, 1);
    gridLayout_3->addWidget(inputAtlasName, 3, 1, 1, 4);
    gridLayout_3->addWidget(inputAtlasButton, 3, 5, 1, 1);

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
    outputDirectoryLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Output Directory:", 0));
    inputfileLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Input File:", 0));
    inputAtlasLabel->setText(QApplication::translate("fPopulationAtlasDialog", "Atlas File:", 0));
  } // retranslateUi

};

namespace Ui {
  class fPopulationAtlasDialog : public ui_fPopulationAtlasDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fPopulationAtlasDialog_H
