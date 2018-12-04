/********************************************************************************
** Form generated from reading UI file 'fJobDialogGLISTR.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FHELPDIALOG_H
#define UI_FHELPDIALOG_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QCheckBox>
// #include <QtGui/QComboBox>
// #include <QtGui/QDoubleSpinBox>
// #include <QtGui/QFrame>
// #include <QtGui/QGridLayout>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QMainWindow>
// #include <QtGui/QMenu>
// #include <QtGui/QMenuBar>
// #include <QtGui/QSlider>
// #include <QtGui/QSpacerItem>
// #include <QtGui/QStatusBar>
// #include <QtGui/QTabWidget>
// #include <QtGui/QTableWidget>
// #include <QtGui/QWidget>
// NEW CHANGES
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QWidget>
#include "QVTKWidget.h"
#include "QVTKOpenGLWidget.h"

//#include "CAPTk.h"

QT_BEGIN_NAMESPACE

class Ui_fHelpDialog
{
public:
  QGridLayout *gridLayout;
  QGridLayout *mainGridLayout;
  QLabel *label1;
  QLabel *label2;
  QLabel *label3;
  QLabel *label4;
  QLabel *label5;
  QLabel *label6;
  QLabel *label7;
  QLabel *label8;
  QLabel *label9;
  QLabel *label1text;
  QLabel *label2text;
  QLabel *label3text;
  QLabel *label4text;
  QLabel *label5text;

  QWidget *TumorTab;
  QGridLayout *gridLayoutT;
  QLabel *label1t;
  QLabel *label2t;
  QLabel *label3t;
  QLabel *label4t;
  QLabel *label5t;
  QLabel *label6t;
  QLabel *label7t;
  QLabel *label9t;
  QLabel *label10t;
  QLabel *label11t;
  QLabel *label12t;
  QLabel *label1ttext;
  QLabel *label2ttext;
  QLabel *label3ttext;
  QLabel *label4ttext;
  QLabel *label5ttext;
  QLabel *label6ttext;
  QLabel *label7ttext;



  QWidget *ShortcutsTab;
  QVBoxLayout * verticalLayoutS;
  QFrame * FrameS;
  QGridLayout *gridLayoutS;
  QLabel *descriptionm;
  QLabel *mouseheading;
  QLabel *keyboardheading;
  QLabel *label1s;
  QLabel *label2s;
  QLabel *label3s;
  QLabel *label4s;
  QLabel *label5s;
  QLabel *label6s;
  QLabel *label7s;
  QLabel *label1stext;
  QLabel *label2stext;
  QLabel *label3stext;
  QLabel *label4stext;
  QLabel *label5stext;
  QLabel *label6stext;
  QLabel *label7stext;

  QWidget *DrawingTab;
  QVBoxLayout * verticalLayoutD;
  QFrame * FrameD;
  QGridLayout *gridLayoutD;
  QLabel *drawing;
  QLabel *erasing;

  QLabel *label1d;
  QLabel *label2d;
  QLabel *label3d;
  QLabel *label4d;
  QLabel *label5d;
  QLabel *label6d;
  QLabel *label7d;
  QLabel *label1dtext;
  QLabel *label2dtext;
  QLabel *label3dtext;
  QLabel *label4dtext;
  QLabel *label5dtext;
  QLabel *label6dtext;
  QLabel *label7dtext;





  QTabWidget *tabWidget;
  QWidget *LoadingTab;
  QWidget *ImagesTab;
  QWidget *shortcutsTab;
  QWidget *aboutTab;

  QWidget *mainLayout;

  QGridLayout *gridLayout_1;

  void setupUi(QDialog *fHelpDialog)
  {
    QSize dialogSize = QSize(600, 800);
    if (fHelpDialog->objectName().isEmpty())
      fHelpDialog->setObjectName(QString::fromUtf8("fHelpDialog"));
    fHelpDialog->setWindowModality(Qt::ApplicationModal);
    fHelpDialog->resize(dialogSize); // needs to be screenSize dependent
    QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(fHelpDialog);
    fHelpDialog->setSizePolicy(sizePolicy);
    fHelpDialog->setMinimumSize(QSize(0, 0));


    mainLayout = new QWidget(fHelpDialog);
    mainLayout->setObjectName(QString::fromUtf8("mainLayout"));
    QSizePolicy sizePolicy1(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy1.setHorizontalStretch(0);
    sizePolicy1.setVerticalStretch(0);
    sizePolicy1.setHeightForWidth(mainLayout->sizePolicy().hasHeightForWidth());
    mainLayout->setSizePolicy(sizePolicy1);



    mainGridLayout = new QGridLayout(mainLayout);
    mainGridLayout->setSpacing(1);
    mainGridLayout->setContentsMargins(1, 1, 1, 1);
    mainGridLayout->setObjectName(QString::fromUtf8("mainGridLayout"));
    mainGridLayout->setSizeConstraint(QLayout::SetMinimumSize);




    tabWidget = new QTabWidget(mainLayout);
    tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
    tabWidget->setEnabled(true);
    //QSizePolicy sizePolicy4(QSizePolicy::Fixed, QSizePolicy::Preferred);
    //sizePolicy4.setHorizontalStretch(0);
    //sizePolicy4.setVerticalStretch(0);
    //sizePolicy4.setHeightForWidth(tabWidget->sizePolicy().hasHeightForWidth());
    //tabWidget->setSizePolicy(sizePolicy4);
    tabWidget->setMaximumSize(QSize(600, 600)); // needs to be screenSize dependent

    //-----------------------loading tab-----------------------------------
    LoadingTab = new QWidget();
    LoadingTab->setObjectName(QString::fromUtf8("LoadingTab"));
    gridLayout = new QGridLayout(LoadingTab);
    gridLayout->setSpacing(0);
    gridLayout->setContentsMargins(0, 0, 0, 0);
    gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

    label1 = new QLabel(LoadingTab);
    label1->setObjectName(QString::fromUtf8("Label1"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1->setSizePolicy(sizePolicy1);
    label2 = new QLabel(LoadingTab);
    label2->setObjectName(QString::fromUtf8("Label2"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2->setSizePolicy(sizePolicy1);
    label3 = new QLabel(LoadingTab);
    label3->setObjectName(QString::fromUtf8("Label3"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3->setSizePolicy(sizePolicy1);
    label4 = new QLabel(LoadingTab);
    label4->setObjectName(QString::fromUtf8("Label4"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4->setSizePolicy(sizePolicy1);
    label5 = new QLabel(LoadingTab);
    label5->setObjectName(QString::fromUtf8("Label5"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5->setSizePolicy(sizePolicy1);

    label6 = new QLabel(LoadingTab);
    label6->setObjectName(QString::fromUtf8("label6"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6->setSizePolicy(sizePolicy1);

    label7 = new QLabel(LoadingTab);
    label7->setObjectName(QString::fromUtf8("label7"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7->setSizePolicy(sizePolicy1);

    label8 = new QLabel(LoadingTab);
    label8->setObjectName(QString::fromUtf8("label8"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label8->setSizePolicy(sizePolicy1);

    label9 = new QLabel(LoadingTab);
    label9->setObjectName(QString::fromUtf8("label9"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label9->setSizePolicy(sizePolicy1);




    label1text = new QLabel(LoadingTab);
    label1text->setObjectName(QString::fromUtf8("label1text"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1text->setSizePolicy(sizePolicy1);


    label2text = new QLabel(LoadingTab);
    label2text->setObjectName(QString::fromUtf8("label2text"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2text->setSizePolicy(sizePolicy1);


    label3text = new QLabel(LoadingTab);
    label3text->setObjectName(QString::fromUtf8("label3text"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3text->setSizePolicy(sizePolicy1);


    label4text = new QLabel(LoadingTab);
    label4text->setObjectName(QString::fromUtf8("label4text"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4text->setSizePolicy(sizePolicy1);

    label5text = new QLabel(LoadingTab);
    label5text->setObjectName(QString::fromUtf8("label5text"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5text->setSizePolicy(sizePolicy1);


    gridLayout->addWidget(label1, 1, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label2, 2, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label3, 3, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label4, 4, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label5, 5, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label6, 6, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label7, 7, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label8, 8, 1, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label9, 9, 1, 1, 1, Qt::AlignTop);


    gridLayout->addWidget(label1text, 1, 2, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label2text, 2, 2, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label3text, 3, 2, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label4text, 4, 2, 1, 1, Qt::AlignTop);
    gridLayout->addWidget(label5text, 5, 2, 1, 1, Qt::AlignTop);


    //-----------------------Tumor tab-----------------------------------
    TumorTab = new QWidget();
    TumorTab->setObjectName(QString::fromUtf8("TumorTab"));
    gridLayoutT = new QGridLayout(TumorTab);
    gridLayoutT->setSpacing(0);
    gridLayoutT->setContentsMargins(0, 0, 0, 0);
    gridLayoutT->setObjectName(QString::fromUtf8("gridLayoutT"));


    label1t = new QLabel(TumorTab);
    label1t->setObjectName(QString::fromUtf8("label1t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1t->setSizePolicy(sizePolicy1);

    label2t = new QLabel(TumorTab);
    label2t->setObjectName(QString::fromUtf8("label2t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2t->setSizePolicy(sizePolicy1);


    label3t = new QLabel(TumorTab);
    label3t->setObjectName(QString::fromUtf8("label3t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3t->setSizePolicy(sizePolicy1);



    label4t = new QLabel(TumorTab);
    label4t->setObjectName(QString::fromUtf8("label3t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4t->setSizePolicy(sizePolicy1);



    label5t = new QLabel(TumorTab);
    label5t->setObjectName(QString::fromUtf8("label5t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5t->setSizePolicy(sizePolicy1);


    label6t = new QLabel(TumorTab);
    label6t->setObjectName(QString::fromUtf8("label6t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6t->setSizePolicy(sizePolicy1);



    label7t = new QLabel(TumorTab);
    label7t->setObjectName(QString::fromUtf8("label7t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7t->setSizePolicy(sizePolicy1);


    label9t = new QLabel(TumorTab);
    label9t->setObjectName(QString::fromUtf8("label9t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label9t->setSizePolicy(sizePolicy1);


    label10t = new QLabel(TumorTab);
    label10t->setObjectName(QString::fromUtf8("label10t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label10t->setSizePolicy(sizePolicy1);



    label11t = new QLabel(TumorTab);
    label11t->setObjectName(QString::fromUtf8("label11t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label11t->setSizePolicy(sizePolicy1);


    label12t = new QLabel(TumorTab);
    label12t->setObjectName(QString::fromUtf8("label12t"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label12t->setSizePolicy(sizePolicy1);


    label1ttext = new QLabel(TumorTab);
    label1ttext->setObjectName(QString::fromUtf8("label1ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1ttext->setSizePolicy(sizePolicy1);


    label2ttext = new QLabel(TumorTab);
    label2ttext->setObjectName(QString::fromUtf8("label2ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2ttext->setSizePolicy(sizePolicy1);


    label3ttext = new QLabel(TumorTab);
    label3ttext->setObjectName(QString::fromUtf8("label3ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3ttext->setSizePolicy(sizePolicy1);


    label4ttext = new QLabel(TumorTab);
    label4ttext->setObjectName(QString::fromUtf8("label4ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4ttext->setSizePolicy(sizePolicy1);


    label5ttext = new QLabel(TumorTab);
    label5ttext->setObjectName(QString::fromUtf8("label5ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5ttext->setSizePolicy(sizePolicy1);

    label6ttext = new QLabel(TumorTab);
    label6ttext->setObjectName(QString::fromUtf8("label6ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6ttext->setSizePolicy(sizePolicy1);

    label7ttext = new QLabel(TumorTab);
    label7ttext->setObjectName(QString::fromUtf8("label7ttext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7ttext->setSizePolicy(sizePolicy1);


    gridLayoutT->addWidget(label2t, 2, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label3t, 3, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label4t, 4, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label5t, 5, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label6t, 6, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label7t, 7, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label9t, 9, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label10t, 10, 1, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label11t, 11, 1, 1, 1, Qt::AlignTop);

    gridLayoutT->addWidget(label2ttext, 2, 2, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label3ttext, 3, 2, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label4ttext, 4, 2, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label5ttext, 5, 2, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label6ttext, 6, 2, 1, 1, Qt::AlignTop);
    gridLayoutT->addWidget(label7ttext, 7, 2, 1, 1, Qt::AlignTop);




    //-----------------------Shortcuts tab-----------------------------------
    ShortcutsTab = new QWidget();
    ShortcutsTab->setObjectName(QString::fromUtf8("ShortcutsTab"));

    verticalLayoutS = new QVBoxLayout(ShortcutsTab);
    verticalLayoutS->setObjectName(QString::fromUtf8("verticalLayout"));




    FrameS = new QFrame(ShortcutsTab);
    FrameS->setObjectName(QString::fromUtf8("FrameS"));

    gridLayoutS = new QGridLayout(FrameS);
    gridLayoutS->setSpacing(2);
    gridLayoutS->setContentsMargins(2, 2, 2, 2);
    gridLayoutS->setObjectName(QString::fromUtf8("gridLayoutS"));

    descriptionm = new QLabel(FrameS);
    descriptionm->setObjectName(QString::fromUtf8("descriptionm"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    descriptionm->setSizePolicy(sizePolicy1);



    mouseheading = new QLabel(FrameS);
    mouseheading->setObjectName(QString::fromUtf8("mouseheading"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    mouseheading->setSizePolicy(sizePolicy1);


    keyboardheading = new QLabel(FrameS);
    keyboardheading->setObjectName(QString::fromUtf8("keyboardheading"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    keyboardheading->setSizePolicy(sizePolicy1);


    label1s = new QLabel(FrameS);
    label1s->setObjectName(QString::fromUtf8("label1s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1s->setSizePolicy(sizePolicy1);


    label2s = new QLabel(FrameS);
    label2s->setObjectName(QString::fromUtf8("label2s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2s->setSizePolicy(sizePolicy1);

    label3s = new QLabel(FrameS);
    label3s->setObjectName(QString::fromUtf8("label3s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3s->setSizePolicy(sizePolicy1);

    label4s = new QLabel(FrameS);
    label4s->setObjectName(QString::fromUtf8("label4s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4s->setSizePolicy(sizePolicy1);

    label5s = new QLabel(FrameS);
    label5s->setObjectName(QString::fromUtf8("label5s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5s->setSizePolicy(sizePolicy1);

    label6s = new QLabel(FrameS);
    label6s->setObjectName(QString::fromUtf8("label6s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6s->setSizePolicy(sizePolicy1);

    label7s = new QLabel(FrameS);
    label7s->setObjectName(QString::fromUtf8("label7s"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7s->setSizePolicy(sizePolicy1);

    label1stext = new QLabel(FrameS);
    label1stext->setObjectName(QString::fromUtf8("label1stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1stext->setSizePolicy(sizePolicy1);

    label2stext = new QLabel(FrameS);
    label2stext->setObjectName(QString::fromUtf8("label2stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2stext->setSizePolicy(sizePolicy1);

    label3stext = new QLabel(ShortcutsTab);
    label3stext->setObjectName(QString::fromUtf8("label3stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3stext->setSizePolicy(sizePolicy1);

    label4stext = new QLabel(FrameS);
    label4stext->setObjectName(QString::fromUtf8("label4stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4stext->setSizePolicy(sizePolicy1);


    label5stext = new QLabel(FrameS);
    label5stext->setObjectName(QString::fromUtf8("label5stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5stext->setSizePolicy(sizePolicy1);

    label6stext = new QLabel(FrameS);
    label6stext->setObjectName(QString::fromUtf8("label6stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6stext->setSizePolicy(sizePolicy1);

    label7stext = new QLabel(FrameS);
    label7stext->setObjectName(QString::fromUtf8("label7stext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7stext->setSizePolicy(sizePolicy1);


    gridLayoutS->addWidget(mouseheading, 6, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label1s, 7, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label2s, 8, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label3s, 9, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(keyboardheading, 1, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label4s, 2, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label5s, 3, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label6s, 4, 1, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label7s, 5, 1, 1, 1, Qt::AlignTop);

    gridLayoutS->addWidget(label1stext, 7, 2, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label2stext, 8, 2, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label3stext, 9, 2, 1, 1, Qt::AlignTop);


    gridLayoutS->addWidget(label4stext, 2, 2, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label5stext, 3, 2, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label6stext, 4, 2, 1, 1, Qt::AlignTop);
    gridLayoutS->addWidget(label7stext, 5, 2, 1, 1, Qt::AlignTop);


    verticalLayoutS->addWidget(FrameS);
    verticalLayoutS->addWidget(descriptionm, Qt::AlignTop);


    //-----------------------Drawing tab-----------------------------------
    DrawingTab = new QWidget();
    DrawingTab->setObjectName(QString::fromUtf8("DrawingTab"));

    verticalLayoutD = new QVBoxLayout(DrawingTab);
    verticalLayoutD->setObjectName(QString::fromUtf8("verticalLayout"));

    FrameD = new QFrame(DrawingTab);
    FrameD->setObjectName(QString::fromUtf8("FrameD"));

    gridLayoutD = new QGridLayout(FrameD);
    gridLayoutD->setSpacing(2);
    gridLayoutD->setContentsMargins(2, 2, 2, 2);
    gridLayoutD->setObjectName(QString::fromUtf8("gridLayoutD"));

    drawing = new QLabel(FrameD);
    drawing->setObjectName(QString::fromUtf8("drawing"));
    sizePolicy1.setHeightForWidth(drawing->sizePolicy().hasHeightForWidth());
    drawing->setSizePolicy(sizePolicy1);

    erasing = new QLabel(FrameD);
    erasing->setObjectName(QString::fromUtf8("erasing"));
    sizePolicy1.setHeightForWidth(erasing->sizePolicy().hasHeightForWidth());
    erasing->setSizePolicy(sizePolicy1);


    label1d = new QLabel(FrameD);
    label1d->setObjectName(QString::fromUtf8("label1d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1d->setSizePolicy(sizePolicy1);



    label2d = new QLabel(FrameD);
    label2d->setObjectName(QString::fromUtf8("label2d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2d->setSizePolicy(sizePolicy1);


    label3d = new QLabel(FrameD);
    label3d->setObjectName(QString::fromUtf8("label3d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3d->setSizePolicy(sizePolicy1);

    label4d = new QLabel(FrameD);
    label4d->setObjectName(QString::fromUtf8("label4d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4d->setSizePolicy(sizePolicy1);


    label5d = new QLabel(FrameD);
    label5d->setObjectName(QString::fromUtf8("label5d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5d->setSizePolicy(sizePolicy1);

    label6d = new QLabel(FrameD);
    label6d->setObjectName(QString::fromUtf8("label6d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6d->setSizePolicy(sizePolicy1);


    label7d = new QLabel(FrameD);
    label7d->setObjectName(QString::fromUtf8("label7d"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7d->setSizePolicy(sizePolicy1);


    label1dtext = new QLabel(FrameD);
    label1dtext->setObjectName(QString::fromUtf8("label1dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label1dtext->setSizePolicy(sizePolicy1);


    label2dtext = new QLabel(FrameD);
    label2dtext->setObjectName(QString::fromUtf8("label2dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label2dtext->setSizePolicy(sizePolicy1);

    label3dtext = new QLabel(FrameD);
    label3dtext->setObjectName(QString::fromUtf8("label3dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label3dtext->setSizePolicy(sizePolicy1);


    label4dtext = new QLabel(FrameD);
    label4dtext->setObjectName(QString::fromUtf8("label4dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label4dtext->setSizePolicy(sizePolicy1);


    label5dtext = new QLabel(FrameD);
    label5dtext->setObjectName(QString::fromUtf8("label5dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label5dtext->setSizePolicy(sizePolicy1);


    label6dtext = new QLabel(FrameD);
    label6dtext->setObjectName(QString::fromUtf8("label6dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label6dtext->setSizePolicy(sizePolicy1);


    label7dtext = new QLabel(FrameD);
    label7dtext->setObjectName(QString::fromUtf8("label7dtext"));
    sizePolicy1.setHeightForWidth(label1->sizePolicy().hasHeightForWidth());
    label7dtext->setSizePolicy(sizePolicy1);


    gridLayoutD->addWidget(drawing, 1, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label1d, 2, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label2d, 3, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label3d, 4, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label4d, 5, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(erasing, 6, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label5d, 7, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label6d, 8, 1, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label7d, 9, 1, 1, 1, Qt::AlignTop);

    gridLayoutD->addWidget(label1dtext, 2, 2, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label2dtext, 3, 2, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label3dtext, 4, 2, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label4dtext, 5, 2, 1, 1, Qt::AlignTop);

    gridLayoutD->addWidget(label5dtext, 7, 2, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label6dtext, 8, 2, 1, 1, Qt::AlignTop);
    gridLayoutD->addWidget(label7dtext, 9, 2, 1, 1, Qt::AlignTop);


    verticalLayoutD->addWidget(FrameD);

    //---------------------------------------------------------------------------------------
    tabWidget->addTab(ShortcutsTab, QString());
    tabWidget->addTab(LoadingTab, QString());
    tabWidget->addTab(TumorTab, QString());
    tabWidget->addTab(DrawingTab, QString());


    mainGridLayout->addWidget(tabWidget, 1, 1, 1, 1);

    QWidget::setTabOrder(tabWidget, label1);
    QWidget::setTabOrder(label1, label2);


    retranslateUi(fHelpDialog);

    QMetaObject::connectSlotsByName(fHelpDialog);
  } // setupUi

  void retranslateUi(QDialog *fHelpDialogs)
  {
    // fHelpDialogs->setWindowTitle(QApplication::translate("fHelpDialog", "Help for Interactions", 0, QApplication::UnicodeUTF8));
    // tabWidget->setTabText(tabWidget->indexOf(LoadingTab), QApplication::translate("fMainWindow", "Images Panel", 0, QApplication::UnicodeUTF8));
    // tabWidget->setTabText(tabWidget->indexOf(TumorTab), QApplication::translate("fMainWindow", "Seed-points Panel", 0, QApplication::UnicodeUTF8));
    // tabWidget->setTabText(tabWidget->indexOf(DrawingTab), QApplication::translate("fMainWindow", "Drawing Panel", 0, QApplication::UnicodeUTF8));
    // tabWidget->setTabText(tabWidget->indexOf(ShortcutsTab), QApplication::translate("fMainWindow", "Controls", 0, QApplication::UnicodeUTF8));
    // NEW CHANGES
    fHelpDialogs->setWindowTitle(QApplication::translate("fHelpDialog", "Help for Interactions", 0));
    tabWidget->setTabText(tabWidget->indexOf(LoadingTab), QApplication::translate("fMainWindow", "Images Panel", 0));
    tabWidget->setTabText(tabWidget->indexOf(TumorTab), QApplication::translate("fMainWindow", "Seed-points Panel", 0));
    tabWidget->setTabText(tabWidget->indexOf(DrawingTab), QApplication::translate("fMainWindow", "Drawing Panel", 0));
    tabWidget->setTabText(tabWidget->indexOf(ShortcutsTab), QApplication::translate("fMainWindow", "Controls", 0));

  }

};

namespace Ui
{
  class fHelpDialog : public Ui_fHelpDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FJOBDIALOGGLISTR_H
