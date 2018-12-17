/********************************************************************************
** Form generated from reading UI file 'fProgressDialog.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FPROGRESSDIALOG_H
#define UI_FPROGRESSDIALOG_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QDialog>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QProgressBar>
// #include <QtGui/QVBoxLayout>
// NEW CHANGES
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QVBoxLayout>


QT_BEGIN_NAMESPACE

class Ui_fProgressDialog
{
public:
  QVBoxLayout *verticalLayout;
  QLabel *textLabel;
  QProgressBar *progressBar;

  void setupUi(QDialog *fProgressDialog)
  {

    if (fProgressDialog->objectName().isEmpty())
      fProgressDialog->setObjectName(QString::fromUtf8("fProgressDialog"));

    fProgressDialog->setWindowModality(Qt::ApplicationModal);
    fProgressDialog->resize(400, 70); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fProgressDialog->sizePolicy().hasHeightForWidth());
    fProgressDialog->setSizePolicy(sizePolicy);
    QFont font;
    font.setFamily(QString::fromUtf8("Calibri"));
    fProgressDialog->setFont(font);
    //fProgressDialog->setModal(true);
    verticalLayout = new QVBoxLayout(fProgressDialog);
    verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
    verticalLayout->setContentsMargins(6, 6, 6, 3);
    textLabel = new QLabel(fProgressDialog);
    textLabel->setObjectName(QString::fromUtf8("textLabel"));
    QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(0);
    sizePolicy1.setHeightForWidth(textLabel->sizePolicy().hasHeightForWidth());
    textLabel->setSizePolicy(sizePolicy1);
    QFont font1;
    font1.setBold(true);
    font1.setWeight(75);
    textLabel->setFont(font1);
    textLabel->setAlignment(Qt::AlignLeading | Qt::AlignLeft | Qt::AlignVCenter);

    verticalLayout->addWidget(textLabel);

    progressBar = new QProgressBar(fProgressDialog);
    progressBar->setObjectName(QString::fromUtf8("progressBar"));
    progressBar->setValue(24);

    verticalLayout->addWidget(progressBar);


    retranslateUi(fProgressDialog);

    QMetaObject::connectSlotsByName(fProgressDialog);
  } // setupUi

  void retranslateUi(QDialog *fProgressDialog)
  {
    // fProgressDialog->setWindowTitle(QApplication::translate("fProgressDialog", "Progress", 0, QApplication::UnicodeUTF8));
    // textLabel->setText(QApplication::translate("fProgressDialog", "Opening image...", 0, QApplication::UnicodeUTF8));
    // NEW CHANGES
    fProgressDialog->setWindowTitle(QApplication::translate("fProgressDialog", "Progress", 0));
    textLabel->setText(QApplication::translate("fProgressDialog", "Opening image...", 0));
  } // retranslateUi

};

namespace Ui {
  class fProgressDialog : public Ui_fProgressDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FPROGRESSDIALOG_H
