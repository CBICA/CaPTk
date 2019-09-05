/**
\file  AppearancePage.cpp

\brief Definition of AppearancePage class

https://www.med.upenn.edu/sbia/software/
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved.
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#include "AppearancePage.h"
#include "ui_AppearancePage.h"
#include <QFontDialog>
#include <QDebug>
#include <QMetaEnum>
#include "CaPTkGUIUtils.h"
#include "cbicaLogging.h"
#include "ApplicationPreferences.h"

AppearancePage::AppearancePage(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::AppearancePage)
{
    ui->setupUi(this);
    ui->currentFontLabel->setText(qApp->font().family());

	//! default
	this->m_PreviousFont = this->m_SelectedFont = qApp->font();
	this->m_PreviousStyleSheet = this->m_SelectedStyleSheet = qApp->styleSheet();
	this->m_PreviousTheme = this->m_SelectedTheme = ThemeType::Light;

    //! connect signals and slots
    connect(ui->selectFontBtn,SIGNAL(clicked()),this,SLOT(OnSelectFontButtonClicked()));
	connect(ui->themeComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(OnChangeTheme(int)));

	//! initialize with DARK theme
	if (!QVariant(ApplicationPreferences::GetInstance()->GetFileAvailability()).toBool())
	{
		ui->themeComboBox->setCurrentIndex(ThemeType::Dark);
		this->OnOkay();
	}

}

AppearancePage::~AppearancePage()
{
    delete ui;
}

void AppearancePage::OnChangeTheme(int theme)
{
	this->m_SelectedTheme = ThemeType(theme);
	if (this->m_SelectedTheme == ThemeType::Dark)
	{
		cbica::Logging(loggerFile, "applying DARK theme");
		std::string qssPath;
#ifdef QT_DEBUG
		qssPath = getCaPTkDataDir() + /*"../etc/"*/ + "/captk.qss";
#else
		qssPath = getCaPTkDataDir() + "../etc/" + CAPTK_STYLESHEET_FILE;
#endif
		QFile f(qssPath.c_str());
		bool success = f.open(QFile::ReadOnly);
		if (success)
		{
			cbica::Logging(loggerFile, "file reading succeeded, qss file path = " + std::string(qssPath.c_str()));
			this->m_SelectedStyleSheet = f.readAll();
		}
		else
		{
			cbica::Logging(loggerFile, "file reading failed, qss file path = " + std::string(qssPath.c_str()));
		}
		f.close();
	}
	else if (this->m_SelectedTheme == ThemeType::Light)
	{
		cbica::Logging(loggerFile, "applying LIGHT theme");
		this->m_SelectedStyleSheet = "";
	}
}

void AppearancePage::OnSelectFontButtonClicked()
{
    bool ok;
    QFont font = QFontDialog::getFont(
                &ok,
                qApp->font(),
                this,
                tr("Pick a font") );
    if( ok )
    {
 		ui->currentFontLabel->setText(font.family());
		this->m_SelectedFont = font;
    }
	else
	{
		ui->currentFontLabel->setText(this->m_PreviousFont.family());
		this->m_SelectedFont = this->m_PreviousFont;
	}
}

void AppearancePage::OnOkay()
{
	//! update previous
	this->m_PreviousFont = this->m_SelectedFont;
	this->m_PreviousStyleSheet = this->m_SelectedStyleSheet;
	this->m_PreviousTheme = this->m_SelectedTheme;

	//! apply selected font and stylesheet
	qApp->setFont(this->m_SelectedFont);
	qApp->setStyleSheet(this->m_SelectedStyleSheet);

	ApplicationPreferences::GetInstance()->SetFont(this->m_SelectedFont.toString());
	ApplicationPreferences::GetInstance()->SetTheme(QVariant::fromValue(this->m_SelectedTheme).toString());
	ApplicationPreferences::GetInstance()->DisplayPreferences();
	
}

void AppearancePage::OnCancel()
{
	//! revert all changes
	ui->currentFontLabel->setText(this->m_PreviousFont.family());
	this->m_SelectedFont = this->m_PreviousFont;
	ui->themeComboBox->setCurrentIndex(this->m_PreviousTheme);

	//! keep previous font, style
	qApp->setFont(this->m_PreviousFont);
	qApp->setStyleSheet(this->m_PreviousStyleSheet);
}

void AppearancePage::Restore()
{
	this->m_SelectedFont.fromString(ApplicationPreferences::GetInstance()->GetFont());
	ui->currentFontLabel->setText(this->m_SelectedFont.family());

	auto metaEnum = QMetaEnum::fromType<ThemeType>();
	this->m_SelectedTheme = static_cast<ThemeType>(metaEnum.keyToValue(
		ApplicationPreferences::GetInstance()->GetTheme().toStdString().c_str()));

	ui->themeComboBox->setCurrentIndex(this->m_SelectedTheme);
	this->OnOkay();
}