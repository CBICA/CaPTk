/********************************************************************************
** Form generated from reading UI file 'fMainWindow.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FMAINWINDOW_H
#define UI_FMAINWINDOW_H

#include <QtCore/QVariant>
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
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include <QtWidgets/QDockWidget>
#include <QSpacerItem>
#include <QGroupBox>
#include <QRadioButton>
#include <QtWidgets/QSpacerItem>

QT_BEGIN_NAMESPACE

class Ui_fMainWindow
{
private:

public:
 
  QSlider *image4DSlider;
  QGroupBox *preferencesGroupBox;
  QWidget *centralwidget;

  QWidget *AxialWidget;
  QWidget *SaggitalWidget;
  QWidget *CoronalWidget;

  QGridLayout * SaggitalWidgetGridLayout;
  QGridLayout * AxialWidgetGridLayout;
  QGridLayout * CoronalWidgetGridLayout;

  QSlider *SaggitalViewSlider;
  QSlider *AxialViewSlider;
  QSlider *CoronalViewSlider;

  QWidget *infoWidget;
  QGridLayout *gridLayout_5;
  QTabWidget *m_tabWidget;
  QDockWidget *m_toolTabdock;

  //preference group box members
  QGridLayout *PrefGridLayout;
  QLabel *thresholdLabel;
  QDoubleSpinBox *levelSpinBox;
  QLabel *presetLabel;
  QDoubleSpinBox *windowSpinBox;
  QLabel *windowLabel;
  QDoubleSpinBox *thresholdSpinBox;
  QComboBox *presetComboBox;
  QLabel *levelLabel;

  QStatusBar *statusbar;
  QGridLayout *overallGridLayout;

  void setupUi(QMainWindow *fMainWindow)
  {

    centralwidget = new QWidget(fMainWindow);
    //--------------------main vertical layout--------------------------
    overallGridLayout = new QGridLayout(centralwidget);
    overallGridLayout->setSpacing(2);
    overallGridLayout->setContentsMargins(2, 2, 2, 2);
    overallGridLayout->setSizeConstraint(QLayout::SetMinimumSize);

    //TopLeveloverallGridLayout: item 1: view widget
    SaggitalWidget = new QWidget(centralwidget);

    QSizePolicy sizePolicy2(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy2.setHorizontalStretch(20);
    sizePolicy2.setVerticalStretch(20);
    sizePolicy2.setHeightForWidth(SaggitalWidget->sizePolicy().hasHeightForWidth());
    SaggitalWidget->setSizePolicy(sizePolicy2);

    SaggitalWidgetGridLayout = new QGridLayout(SaggitalWidget);
    SaggitalWidgetGridLayout->setContentsMargins(0, 0, 0, 0);
    SaggitalWidgetGridLayout->setSizeConstraint(QLayout::SetMinimumSize);
    SaggitalWidgetGridLayout->setHorizontalSpacing(0);
    SaggitalViewSlider = new QSlider(SaggitalWidget);
    SaggitalViewSlider->setOrientation(Qt::Vertical);
    SaggitalWidgetGridLayout->addWidget(SaggitalViewSlider, 0, 1, 1, 1);
    overallGridLayout->addWidget(SaggitalWidget, 2, 2, 1, 1);

    AxialWidget = new QWidget(centralwidget);
    QSizePolicy sizePolicy6(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy6.setHorizontalStretch(20);
    sizePolicy6.setVerticalStretch(20);
    sizePolicy6.setHeightForWidth(AxialWidget->sizePolicy().hasHeightForWidth());
    AxialWidget->setSizePolicy(sizePolicy6);
    AxialWidget->setMinimumSize(QSize(256, 256));//TBD not working

    AxialWidgetGridLayout = new QGridLayout(AxialWidget);
    AxialWidgetGridLayout->setContentsMargins(0, 0, 0, 0);
    AxialWidgetGridLayout->setSizeConstraint(QLayout::SetMinimumSize);
    AxialWidgetGridLayout->setHorizontalSpacing(0);

    AxialViewSlider = new QSlider(AxialWidget);
    AxialViewSlider->setOrientation(Qt::Vertical);
    AxialWidgetGridLayout->addWidget(AxialViewSlider, 0, 1, 1, 1);
    overallGridLayout->addWidget(AxialWidget, 2, 1, 1, 1);

    CoronalWidget = new QWidget(centralwidget);
    QSizePolicy sizePolicy7(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sizePolicy7.setHorizontalStretch(20);
    sizePolicy7.setVerticalStretch(20);
    sizePolicy7.setHeightForWidth(CoronalWidget->sizePolicy().hasHeightForWidth());
    CoronalWidget->setSizePolicy(sizePolicy7);
    CoronalWidgetGridLayout = new QGridLayout(CoronalWidget);
    CoronalWidgetGridLayout->setContentsMargins(0, 0, 0, 0);
    CoronalWidgetGridLayout->setSizeConstraint(QLayout::SetMinimumSize);
    CoronalWidgetGridLayout->setHorizontalSpacing(0);

    CoronalViewSlider = new QSlider(CoronalWidget);
    CoronalViewSlider->setOrientation(Qt::Vertical);
    CoronalWidgetGridLayout->addWidget(CoronalViewSlider, 0, 1, 1, 1);

    overallGridLayout->addWidget(CoronalWidget, 2, 0, 1, 1);
    //-------------------------------second line of fMainWindow---------------------------
    image4DSlider = new QSlider(centralwidget);
    image4DSlider->setMaximum(45);
    image4DSlider->setPageStep(1);
    image4DSlider->setSliderPosition(5);
    image4DSlider->setOrientation(Qt::Horizontal);
    overallGridLayout->addWidget(image4DSlider, 1, 0, 1, 3);
    //---------------------------------------------------------------------------------------
    preferencesGroupBox = new QGroupBox(centralwidget);
    //-------------------preferences related objects----------------------

    PrefGridLayout = new QGridLayout(preferencesGroupBox);

    thresholdLabel = new QLabel(preferencesGroupBox);
    levelSpinBox = new QDoubleSpinBox(preferencesGroupBox);
    levelSpinBox->setDecimals(3);
    levelSpinBox->setMinimum(-66000);
    levelSpinBox->setMaximum(66000);
    levelSpinBox->setSingleStep(10);

    presetLabel = new QLabel(preferencesGroupBox);
    windowSpinBox = new QDoubleSpinBox(preferencesGroupBox);
    windowSpinBox->setDecimals(3);
    windowSpinBox->setMinimum(-66000);
    windowSpinBox->setMaximum(66000);
    windowSpinBox->setSingleStep(10);

    windowLabel = new QLabel(preferencesGroupBox);
    thresholdSpinBox = new QDoubleSpinBox(preferencesGroupBox);
    thresholdSpinBox->setDecimals(3);
    thresholdSpinBox->setMinimum(-99999);
    thresholdSpinBox->setMaximum(99999);
    thresholdSpinBox->setSingleStep(1);

    presetComboBox = new QComboBox(preferencesGroupBox);
    QSizePolicy sizePolicy8(QSizePolicy::Minimum, QSizePolicy::Fixed);
    sizePolicy8.setHorizontalStretch(0);
    sizePolicy8.setVerticalStretch(0);
    sizePolicy8.setHeightForWidth(presetComboBox->sizePolicy().hasHeightForWidth());
    presetComboBox->setSizePolicy(sizePolicy8);
    presetComboBox->setMaximumSize(QSize(100, 16777215));
    levelLabel = new QLabel(preferencesGroupBox);
    levelLabel->setAlignment(Qt::AlignCenter);
    thresholdLabel->setAlignment(Qt::AlignCenter);

    PrefGridLayout->addWidget(windowLabel, 0, 0, 1, 1);
    PrefGridLayout->addWidget(windowSpinBox, 0, 1, 1, 1);

    PrefGridLayout->addWidget(levelLabel, 0, 2, 1, 1);
    PrefGridLayout->addWidget(levelSpinBox, 0, 3, 1, 1);

    PrefGridLayout->addWidget(presetLabel, 1, 0, 1, 1);
    PrefGridLayout->addWidget(presetComboBox, 1, 1, 1, 1);


    PrefGridLayout->addWidget(thresholdLabel, 1, 2, 1, 1);
    PrefGridLayout->addWidget(thresholdSpinBox, 1, 3, 1, 1);

    m_tabWidget = new QTabWidget();
    m_toolTabdock = new QDockWidget();
    statusbar = new QStatusBar();

    retranslateUi(fMainWindow);

    QMetaObject::connectSlotsByName(fMainWindow);
  }

  void retranslateUi(QMainWindow *fMainWindow)
  {
 
    fMainWindow->setWindowTitle(QApplication::translate("fMainWindow", "fMainWindow", 0));
    thresholdLabel->setText(QApplication::translate("fMainWindow", "Threshold", 0));
    presetLabel->setText(QApplication::translate("fMainWindow", "Preset", 0));
    windowLabel->setText(QApplication::translate("fMainWindow", "Window", 0));
    levelLabel->setText(QApplication::translate("fMainWindow", "Level", 0));
  }

};

namespace Ui {
  class fMainWindow : public Ui_fMainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FMAINWINDOW_H
