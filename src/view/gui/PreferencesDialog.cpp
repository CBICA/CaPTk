/**
\file  PreferencesDialog.cpp

\brief Implementation of PreferencesDialog class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "PreferencesDialog.h"
#include "AppearancePage.h"

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QStackedWidget>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

#include <QFontDialog>
#include <QColorDialog>
#include <QSizePolicy>
#include <QDebug>

PreferencesDialog::PreferencesDialog(QWidget *parent) :
    QDialog(parent)
{
    this->SetupUi();

	this->m_AppearancePage = new AppearancePage(this->stackedWidget);
    this->m_FontDialog = new QFontDialog(this->stackedWidget);
    this->m_FontDialog->setWindowFlags(Qt::SubWindow);
    this->m_FontDialog->setOptions(
                /* do not use native dialog */
                QFontDialog::DontUseNativeDialog
                /* you don't need to set it, but if you don't set this
                                        the "OK" and "Cancel" buttons will show up, I don't
                                        think you'd want that. */
                | QFontDialog::NoButtons
                );
    this->m_FontDialog->setSizeGripEnabled(false);

	this->stackedWidget->insertWidget(0, this->m_AppearancePage);
    //this->stackedWidget->insertWidget(1,this->m_FontDialog);

	//! signals and slots
    connect(this->listWidget,SIGNAL(itemSelectionChanged()),this,SLOT(OnItemSelectionChanged()));

}

PreferencesDialog::~PreferencesDialog()
{
}

QFontDialog *PreferencesDialog::GetFontDialog()
{
    return this->m_FontDialog;
}

void PreferencesDialog::OnItemSelectionChanged()
{
    int currentItemIndex = this->listWidget->currentRow();
    this->stackedWidget->setCurrentIndex(currentItemIndex);
}

void PreferencesDialog::SetupUi()
{
	if (this->objectName().isEmpty())
		this->setObjectName(QStringLiteral("PreferencesDialog"));
	this->resize(800, 400);
	verticalLayout = new QVBoxLayout(this);
	verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
	widget = new QWidget(this);
	widget->setObjectName(QStringLiteral("widget"));
	widget->setMinimumSize(QSize(0, 0));
	horizontalLayout = new QHBoxLayout(widget);
	horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
	listWidget = new QListWidget(widget);
	new QListWidgetItem(listWidget);
	new QListWidgetItem(listWidget);

	listWidget->setObjectName(QStringLiteral("listWidget"));
	listWidget->setFixedWidth(175);
	listWidget->setMaximumSize(QSize(16777215, 16777215));

	horizontalLayout->addWidget(listWidget);

	stackedWidget = new QStackedWidget(widget);
	stackedWidget->setObjectName(QStringLiteral("stackedWidget"));
	stackedWidget->setStyleSheet(QStringLiteral(""));

	horizontalLayout->addWidget(stackedWidget);
	
	verticalLayout->addWidget(widget);

	buttonBox = new QDialogButtonBox(this);
	buttonBox->setObjectName(QStringLiteral("buttonBox"));
	buttonBox->setOrientation(Qt::Horizontal);
	buttonBox->setStandardButtons(QDialogButtonBox::Cancel | QDialogButtonBox::Ok);

	verticalLayout->addWidget(buttonBox);
	
	retranslateUi(this);
	QObject::connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	QObject::connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

	stackedWidget->setCurrentIndex(-1);
	
	QMetaObject::connectSlotsByName(this);
}

void PreferencesDialog::retranslateUi(QDialog *PreferencesDialog)
{
	this->setWindowTitle(QApplication::translate("PreferencesDialog", "Preferences", nullptr));

	const bool __sortingEnabled = listWidget->isSortingEnabled();
	listWidget->setSortingEnabled(false);
	QListWidgetItem *___qlistwidgetitem = listWidget->item(0);
	___qlistwidgetitem->setText(QApplication::translate("PreferencesDialog", "Appearance", nullptr));
	QListWidgetItem *___qlistwidgetitem1 = listWidget->item(1);
	___qlistwidgetitem1->setText(QApplication::translate("PreferencesDialog", "Fonts", nullptr));
	listWidget->setSortingEnabled(__sortingEnabled);
}