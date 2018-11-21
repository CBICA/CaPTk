/********************************************************************************
** Form generated from reading UI file 'fBottomImageInfoTip.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_fBOTTOMIMAGEINFOTIP_H
#define UI_fBOTTOMIMAGEINFOTIP_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QFrame>
// #include <QtGui/QGridLayout>
// #include <QtGui/QHBoxLayout>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QLineEdit>
// #include <QtGui/QPushButton>
// #include <QtGui/QScrollArea>
// #include <QtGui/QSpacerItem>
// #include <QtGui/QVBoxLayout>
// #include <QtGui/QWidget>
// NEW CHANGES
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_fBottomImageInfoTip
{
public:
  QWidget *centralwidget;

  QLabel *imageLabel;
  QLabel *sizePixelLabel;
  QLabel *spacingLabel;
  QLabel *originLabel;
  QLabel *valueLabel;
  QPushButton *pixelPosButton;
  QLineEdit *pixelPosX;
  QLineEdit *pixelPosZ;
  QLineEdit *pixelPosY;


  void setupUi(QWidget *parent)
  {
    QVBoxLayout *mainLayout = new QVBoxLayout(parent);
    QHBoxLayout *subLayout1 = new QHBoxLayout();
    QHBoxLayout *subLayout2 = new QHBoxLayout();
    //int buttonWidth = QLabel().fontMetrics().width("-----------------");
    int pixelButtonWidth = QLabel().fontMetrics().width("---------");
    imageLabel = new QLabel();
    imageLabel->setFixedWidth(0);//TBD fix this if want to show
    sizePixelLabel = new QLabel();
    //sizePixelLabel->setFixedWidth(pixelButtonWidth); // [TBD] I don't think these should be fixed
    spacingLabel = new QLabel();
    //spacingLabel->setFixedWidth(buttonWidth);
    originLabel = new QLabel();
    //originLabel->setFixedWidth(buttonWidth);
    valueLabel = new QLabel();
    //valueLabel->setFixedWidth(pixelButtonWidth);
    pixelPosX = new QLineEdit();
    pixelPosX->setFixedWidth(pixelButtonWidth);
    pixelPosX->setToolTip("Right Panel");
    pixelPosX->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    pixelPosY = new QLineEdit();
    pixelPosY->setFixedWidth(pixelButtonWidth);
    pixelPosY->setToolTip("Left Panel");
    pixelPosY->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    pixelPosZ = new QLineEdit();
    pixelPosZ->setFixedWidth(pixelButtonWidth);
    pixelPosZ->setToolTip("Center Panel");
    pixelPosZ->setAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    pixelPosButton = new QPushButton("Update");
    pixelPosButton->setFixedWidth(2 * pixelButtonWidth);

    subLayout1->addWidget(new QLabel("[X, Y, Z]   |   Image Size:"));
    subLayout1->addWidget(sizePixelLabel);
    subLayout1->addWidget(new QLabel("|   Origin:"));
    subLayout1->addWidget(originLabel);
    subLayout1->addWidget(new QLabel("|   Spacing:"));
    subLayout1->addWidget(spacingLabel);
    //subLayout1->addWidget(new QLabel("Name:"));
    subLayout1->addWidget(imageLabel);
    subLayout1->addStretch();


    subLayout2->addWidget(new QLabel("Current Slice:"));
    subLayout2->addWidget(new QLabel("Y:"));
    subLayout2->addWidget(pixelPosY);
    subLayout2->addWidget(new QLabel("Z:"));
    subLayout2->addWidget(pixelPosZ);
    subLayout2->addWidget(new QLabel("X:"));
    subLayout2->addWidget(pixelPosX);
    subLayout2->addWidget(pixelPosButton);
    subLayout2->addWidget(new QLabel("Value:"));
    subLayout2->addWidget(valueLabel);
    subLayout2->addStretch();

    mainLayout->addLayout(subLayout1);
    mainLayout->addLayout(subLayout2);
    QMetaObject::connectSlotsByName(parent);
  } // setupUi



};

namespace Ui {
  class fBottomImageInfoTip : public Ui_fBottomImageInfoTip {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_fBottomImageInfoTip_H
