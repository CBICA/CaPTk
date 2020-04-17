#ifndef ui_fAppDownloadDialog_H
#define ui_fAppDownloadDialog_H

#include <QtCore/QVariant>
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
#include <QProgressDialog>

#include "CaPTkGUIUtils.h"

QT_BEGIN_NAMESPACE

class ui_fAppDownloadDialog
{
public:
  QGridLayout *gridLayout;

  QPushButton *confirmButton;
  QPushButton *cancelButton;

  QProgressDialog *progressDialog;

  void setupUi(QDialog *fAppDownloadDialog)
  {

    if (fAppDownloadDialog->objectName().isEmpty())
      fAppDownloadDialog->setObjectName(QString::fromUtf8("fAppDownloadDialog"));

    fAppDownloadDialog->setWindowModality(Qt::NonModal);
    fAppDownloadDialog->resize(100, 100); // needs to be screenSize dependent 
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fAppDownloadDialog->sizePolicy().hasHeightForWidth());
    fAppDownloadDialog->setSizePolicy(sizePolicy);
    fAppDownloadDialog->setMinimumSize(QSize(0, 0));

    //--------------------------------------------------------------------
    gridLayout = new QGridLayout(fAppDownloadDialog);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    QLabel *message = new QLabel(fAppDownloadDialog);
    message->setObjectName(QString::fromUtf8("message"));
    message->setText("This application has not been installed. Do you want to install it?");

    gridLayout->addWidget(message, 1, 0, 1, 2);

    confirmButton = new QPushButton(fAppDownloadDialog);
    confirmButton->setObjectName(QString::fromUtf8("Ok"));
    //confirmButton->setIcon(ButtonIcon);
    //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    cancelButton = new QPushButton(fAppDownloadDialog);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
    //cancelButton->setIcon(ButtonIcon);
    //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 

    gridLayout->addWidget(confirmButton, 2, 0, 1, 1);
    gridLayout->addWidget(cancelButton, 2, 1, 1, 1);

    retranslateUi(fAppDownloadDialog);

    QMetaObject::connectSlotsByName(fAppDownloadDialog);
  } // setupUi

  void setupDownload(QDialog *fAppDownloadDialog) {
    progressDialog = new QProgressDialog(fAppDownloadDialog);
    progressDialog->setObjectName(QString::fromUtf8("ProgressDialog"));
  }

  void retranslateUi(QDialog *fAppDownloadDialog)
  {
    fAppDownloadDialog->setWindowTitle(QApplication::translate("fAppDownloadDialog", "Applications Manager", 0));
    confirmButton->setText(QApplication::translate("fAppDownloadDialog", "Ok", 0));
    cancelButton->setText(QApplication::translate("fAppDownloadDialog", "Cancel", 0));
  } // retranslateUi

};

namespace Ui {
  class fAppDownloadDialog : public ui_fAppDownloadDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fAppDownloadDialog_H