/********************************************************************************
** Form generated from reading UI file 'fTumorPanel.ui'
**
** Created by: Qt User Interface Compiler version 4.8.6
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FTUMORPANEL_H
#define UI_FTUMORPANEL_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QTableWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <QtWidgets/QGroupBox>
#include <QDir>
#include "CaPTkGUIUtils.h"

struct ButtonStructure
{
  QRadioButton *radioButton;
  std::string text_toolTip, text_label;
  int enum_label;
  size_t counter;
};

QT_BEGIN_NAMESPACE

class Ui_fTumorPanel
{
public:
  // type selector buttons
  QRadioButton *	m_typeRadBtnTumor;
  QRadioButton *	m_typeRadBtnAllTissues;
  QRadioButton *	m_typeRadBtnGlister;
  QRadioButton *	m_typeRadBtnPorterPre;
  QRadioButton *	m_typeRadBtnPorterPost;
  QPushButton * HelpButton;

  std::vector< ButtonStructure > m_vectorOfSeedPointLabels; // genericLabels

  QTableWidget *sTableWidget; // tumor point table
  QPushButton *sRemoveButton; // tumor point remove selected
  QPushButton *sRemoveAllButton; // tumor point remove all
  QPushButton *sLoadButton; // tumor point load from previous
  QPushButton *sSaveButton; // tumor point save to file
  QTableWidget *tTableWidget; // tissue point table
  QVBoxLayout *finalVerticalLayout;
  QPushButton *tRemoveButton; // tissue point remove selected
  QPushButton *tRemoveAllButton; // tissue point remove all
  QPushButton *tLoadButton; // tissue point load from previous
  QPushButton *tSaveButton; // tissue point save to file

  void setupUi(QWidget *parent)
  {
    m_typeRadBtnTumor = new QRadioButton();
    m_typeRadBtnTumor->setObjectName(QString::fromUtf8("m_typeRadBtnTumor"));
    m_typeRadBtnTumor->setFocusPolicy(Qt::TabFocus);
    m_typeRadBtnTumor->setToolTip("Define 'Tumor Points' using coordinate and radius information");

    m_typeRadBtnAllTissues = new QRadioButton();
    m_typeRadBtnAllTissues->setObjectName(QString::fromUtf8("m_typeRadBtnAllTissues"));
    m_typeRadBtnAllTissues->setFocusPolicy(Qt::TabFocus);
    m_typeRadBtnAllTissues->setToolTip("Define 'Tissue Points' using only coordinate");

    m_typeRadBtnGlister = new QRadioButton();
    m_typeRadBtnGlister->setObjectName(QString::fromUtf8("m_typeRadBtnGlister"));
    m_typeRadBtnGlister->setFocusPolicy(Qt::TabFocus);

    m_typeRadBtnPorterPre = new QRadioButton();
    m_typeRadBtnPorterPre->setObjectName(QString::fromUtf8("m_typeRadBtnPorterPre"));
    m_typeRadBtnPorterPre->setFocusPolicy(Qt::TabFocus);

    m_typeRadBtnPorterPost = new QRadioButton();
    m_typeRadBtnPorterPost->setObjectName(QString::fromUtf8("m_typeRadBtnPorterPost"));
    m_typeRadBtnPorterPost->setFocusPolicy(Qt::TabFocus);

    sTableWidget = new QTableWidget();
    if (sTableWidget->columnCount() < 4)
      sTableWidget->setColumnCount(4);
    QTableWidgetItem *__qtablewidgetitem = new QTableWidgetItem();
    __qtablewidgetitem->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    sTableWidget->setHorizontalHeaderItem(0, __qtablewidgetitem);
    QTableWidgetItem *__qtablewidgetitem1 = new QTableWidgetItem();
    __qtablewidgetitem1->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    sTableWidget->setHorizontalHeaderItem(1, __qtablewidgetitem1);
    QTableWidgetItem *__qtablewidgetitem2 = new QTableWidgetItem();
    __qtablewidgetitem2->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    sTableWidget->setHorizontalHeaderItem(2, __qtablewidgetitem2);
    QTableWidgetItem *__qtablewidgetitem3 = new QTableWidgetItem();
    __qtablewidgetitem3->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    sTableWidget->setHorizontalHeaderItem(3, __qtablewidgetitem3);
    sTableWidget->setObjectName(QString::fromUtf8("sTableWidget"));
    sTableWidget->setAutoScroll(true);
    sTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
    sTableWidget->setTabKeyNavigation(false);
    sTableWidget->setAlternatingRowColors(false);
    sTableWidget->setTextElideMode(Qt::ElideNone);
    sTableWidget->setWordWrap(false);
    sTableWidget->horizontalHeader()->setHighlightSections(false);
    sTableWidget->verticalHeader()->setHighlightSections(false);
    //http://stackoverflow.com/questions/13025579/qt-how-can-i-set-qtablewidget-row-index-starting-from-64-instead-of-starting-w

    sRemoveButton = new QPushButton();
    sRemoveButton->setObjectName(QString::fromUtf8("sRemoveButton"));
    sRemoveButton->setFocusPolicy(Qt::TabFocus);
    sRemoveAllButton = new QPushButton();
    sRemoveAllButton->setObjectName(QString::fromUtf8("sRemoveAllButton"));
    sRemoveAllButton->setFocusPolicy(Qt::TabFocus);
    sLoadButton = new QPushButton();
    sLoadButton->setObjectName(QString::fromUtf8("sLoadButton"));
    sLoadButton->setFocusPolicy(Qt::TabFocus);
    sSaveButton = new QPushButton();
    sSaveButton->setObjectName(QString::fromUtf8("sSaveButton"));
    sSaveButton->setFocusPolicy(Qt::TabFocus);

    QHBoxLayout * seedLayout = new QHBoxLayout();
    QVBoxLayout * seedLayoutSub = new QVBoxLayout();
    seedLayoutSub->addWidget(sRemoveButton);
    seedLayoutSub->addWidget(sRemoveAllButton);
    seedLayoutSub->addWidget(sLoadButton);
    seedLayoutSub->addWidget(sSaveButton);
    seedLayoutSub->addStretch();

    seedLayout->addWidget(sTableWidget);
    seedLayout->addLayout(seedLayoutSub);

    QGroupBox*  seedGroupBox = new QGroupBox();
    seedGroupBox->setTitle(QString::fromStdString("Tumor Points [World Coordinates]"));
    seedGroupBox->setLayout(seedLayout);
    //----------------------------------------------------------------------------

    tTableWidget = new QTableWidget();
    if (tTableWidget->columnCount() < 4)
      tTableWidget->setColumnCount(4);
    QTableWidgetItem *__qtablewidgetitem15 = new QTableWidgetItem();
    __qtablewidgetitem15->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    tTableWidget->setHorizontalHeaderItem(0, __qtablewidgetitem15);
    QTableWidgetItem *__qtablewidgetitem16 = new QTableWidgetItem();
    __qtablewidgetitem16->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    tTableWidget->setHorizontalHeaderItem(1, __qtablewidgetitem16);
    QTableWidgetItem *__qtablewidgetitem17 = new QTableWidgetItem();
    __qtablewidgetitem17->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    tTableWidget->setHorizontalHeaderItem(2, __qtablewidgetitem17);
    QTableWidgetItem *__qtablewidgetitem57 = new QTableWidgetItem();
    __qtablewidgetitem57->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter | Qt::AlignCenter);
    tTableWidget->setHorizontalHeaderItem(3, __qtablewidgetitem57);

    tTableWidget->setObjectName(QString::fromUtf8("tTableWidget"));
    tTableWidget->setAutoScroll(true);
    tTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);
    tTableWidget->setTabKeyNavigation(false);
    tTableWidget->setAlternatingRowColors(false);
    tTableWidget->setTextElideMode(Qt::ElideNone);
    tTableWidget->setWordWrap(false);
    tTableWidget->horizontalHeader()->setHighlightSections(false);
    tTableWidget->verticalHeader()->setHighlightSections(false);

    //----------------------------------------------------
    QGridLayout*  tissueTypeLayout = new QGridLayout();;
    m_vectorOfSeedPointLabels.resize(NumberOfTissueTypes); // genericLabels
    int tissueTypeLayoutRowCounter = 1;

    QLabel * tissueTypeLabel = new QLabel("Types :");
    tissueTypeLayout->addWidget(tissueTypeLabel, tissueTypeLayoutRowCounter, 1, 1, 4);
    tissueTypeLayoutRowCounter++;

    // genericLabels
    for (size_t i = 0; i < NumberOfTissueTypes; tissueTypeLayoutRowCounter++)
    {
      for (size_t colCounter = 1; colCounter < 6; colCounter++, i++)
      {
        if (i < NumberOfTissueTypes)
        {
          m_vectorOfSeedPointLabels[i].radioButton = new QRadioButton();
          tissueTypeLayout->addWidget(m_vectorOfSeedPointLabels[i].radioButton, tissueTypeLayoutRowCounter, colCounter, 1, 1);
        }
      }
    }

    tRemoveButton = new QPushButton();
    tRemoveButton->setObjectName(QString::fromUtf8("tRemoveButton"));
    tRemoveButton->setFocusPolicy(Qt::TabFocus);
    tRemoveAllButton = new QPushButton();
    tRemoveAllButton->setObjectName(QString::fromUtf8("tRemoveAllButton"));
    tRemoveAllButton->setFocusPolicy(Qt::TabFocus);
    tLoadButton = new QPushButton();
    tLoadButton->setObjectName(QString::fromUtf8("tLoadButton"));
    tLoadButton->setFocusPolicy(Qt::TabFocus);
    tSaveButton = new QPushButton();
    tSaveButton->setObjectName(QString::fromUtf8("tSaveButton"));
    tSaveButton->setFocusPolicy(Qt::TabFocus);

    QHBoxLayout * tissueLayou = new QHBoxLayout();
    QHBoxLayout * tissueLayoutSub1 = new QHBoxLayout();
    QVBoxLayout * tissueLayoutSubSub1 = new QVBoxLayout();
    tissueLayoutSubSub1->addWidget(tRemoveButton);
    tissueLayoutSubSub1->addWidget(tRemoveAllButton);
    tissueLayoutSubSub1->addWidget(tLoadButton);
    tissueLayoutSubSub1->addWidget(tSaveButton);
    tissueLayoutSubSub1->addStretch();

    tissueLayoutSub1->addWidget(tTableWidget);
    tissueLayoutSub1->addLayout(tissueLayoutSubSub1);

    QVBoxLayout * tissueTypeLayoutV = new QVBoxLayout();
    tissueTypeLayoutV->addLayout(tissueTypeLayout);
    tissueTypeLayoutV->addStretch();
    tissueLayou->addLayout(tissueTypeLayoutV);
    tissueLayou->addLayout(tissueLayoutSub1);


    QGroupBox* tissueGroupBox = new QGroupBox();
    tissueGroupBox->setTitle(QString::fromStdString("Tissue Points [World Coordinates]"));
    tissueGroupBox->setLayout(tissueLayou);

    HelpButton = new QPushButton();
    std::string iconDir;

    iconDir = getCaPTkDataDir() + "/icons/";
    if (!QDir(cbica::normPath(iconDir).c_str()).exists()) // packaged binary
    {
      if (QDir((std::string(PROJECT_SOURCE_DIR) + "/data/icons/").c_str()).exists()) // developer_mode
      {
        iconDir = std::string(PROJECT_SOURCE_DIR) + "/data/icons/";
      }
    }

    HelpButton->setIcon(QIcon((iconDir + "help.png").c_str()));
    HelpButton->setToolTip("Get Help");
    //HelpButton->setIconSize(QSize(30, 30));
    QHBoxLayout * radBtnLayout = new QHBoxLayout();
    radBtnLayout->addWidget(new QLabel("Point Types: "));
    radBtnLayout->addWidget(m_typeRadBtnTumor);
    radBtnLayout->addStretch();
    radBtnLayout->addStretch();
    radBtnLayout->addWidget(m_typeRadBtnAllTissues);
    radBtnLayout->addWidget(m_typeRadBtnGlister);
    radBtnLayout->addWidget(m_typeRadBtnPorterPre);
    radBtnLayout->addWidget(m_typeRadBtnPorterPost);
    radBtnLayout->addStretch();
    radBtnLayout->addWidget(HelpButton);

    QSizePolicy spLeft(QSizePolicy::Preferred, QSizePolicy::Preferred);
    spLeft.setHorizontalStretch(4);
    QSizePolicy spRight(QSizePolicy::Preferred, QSizePolicy::Preferred);
    spRight.setHorizontalStretch(5);

    seedGroupBox->setSizePolicy(spLeft);
    tissueGroupBox->setSizePolicy(spRight);

    QHBoxLayout* mainLayoutSub = new QHBoxLayout();
    mainLayoutSub->addWidget(seedGroupBox);
    mainLayoutSub->addWidget(tissueGroupBox);

    QVBoxLayout* mainLayout = new QVBoxLayout(parent);
    mainLayout->addLayout(radBtnLayout);
    mainLayout->addLayout(mainLayoutSub);
    retranslateUi(parent);

  }

  void retranslateUi(QWidget *fTumorPanel)
  {
    fTumorPanel->setWindowTitle(QApplication::translate("fTumorPanel", "Form", 0));
    m_typeRadBtnTumor->setText(QApplication::translate("fTumorPanel", "Tumor Points", 0));
    m_typeRadBtnAllTissues->setText(QApplication::translate("fTumorPanel", "All Tissue Points", 0));

    m_typeRadBtnGlister->setText(QApplication::translate("fTumorPanel", "GLISTR/GLISTRboost", 0));
    m_typeRadBtnPorterPre->setText(QApplication::translate("fTumorPanel", "PORTR-PRE", 0));
    m_typeRadBtnPorterPost->setText(QApplication::translate("fTumorPanel", "PORTR-POST", 0));


    auto sHeader = sTableWidget->horizontalHeader();
    sHeader->setSectionResizeMode(QHeaderView::Stretch);
    QTableWidgetItem *___qtablewidgetitem = sTableWidget->horizontalHeaderItem(0);
    ___qtablewidgetitem->setText(QApplication::translate("fTumorPanel", "      x      ", 0));
    QTableWidgetItem *___qtablewidgetitem1 = sTableWidget->horizontalHeaderItem(1);
    ___qtablewidgetitem1->setText(QApplication::translate("fTumorPanel", "      y      ", 0));
    QTableWidgetItem *___qtablewidgetitem2 = sTableWidget->horizontalHeaderItem(2);
    ___qtablewidgetitem2->setText(QApplication::translate("fTumorPanel", "      z      ", 0));
    QTableWidgetItem *___qtablewidgetitem3 = sTableWidget->horizontalHeaderItem(3);
    ___qtablewidgetitem3->setText(QApplication::translate("fTumorPanel", "      r      ", 0));
    const bool __sortingEnabled = sTableWidget->isSortingEnabled();
    sTableWidget->setSortingEnabled(false);
    sTableWidget->setSortingEnabled(__sortingEnabled);
    sTableWidget->verticalHeader()->setVisible(false);
    sRemoveButton->setText(QApplication::translate("fTumorPanel", "Remove", 0));
    sRemoveAllButton->setText(QApplication::translate("fTumorPanel", "Clear all", 0));
    sLoadButton->setText(QApplication::translate("fTumorPanel", "Load", 0));
    sSaveButton->setText(QApplication::translate("fTumorPanel", "Save", 0));

    auto tHeader = tTableWidget->horizontalHeader();
    tHeader->setSectionResizeMode(QHeaderView::Stretch);
    QTableWidgetItem *___qtablewidgetitem14 = tTableWidget->horizontalHeaderItem(1);
    ___qtablewidgetitem14->setText(QApplication::translate("fTumorPanel", "      x      ", 0));
    QTableWidgetItem *___qtablewidgetitem15 = tTableWidget->horizontalHeaderItem(2);
    ___qtablewidgetitem15->setText(QApplication::translate("fTumorPanel", "      y      ", 0));
    QTableWidgetItem *___qtablewidgetitem16 = tTableWidget->horizontalHeaderItem(3);
    ___qtablewidgetitem16->setText(QApplication::translate("fTumorPanel", "      z      ", 0));

    tRemoveButton->setText(QApplication::translate("fTumorPanel", "Remove", 0));
    tRemoveAllButton->setText(QApplication::translate("fTumorPanel", "Clear all", 0));
    tLoadButton->setText(QApplication::translate("fTumorPanel", "Load", 0));
    tSaveButton->setText(QApplication::translate("fTumorPanel", "Save", 0));

    for (size_t i = 0; i < NumberOfTissueTypes; i++)
    {
      m_vectorOfSeedPointLabels[i].radioButton->setText(m_vectorOfSeedPointLabels[i].text_label.c_str());
    }

  }

};

namespace Ui {
  class fTumorPanel : public Ui_fTumorPanel {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FTUMORPANEL_H