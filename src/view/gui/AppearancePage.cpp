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
#include "CaPTkGUIUtils.h"
#include "cbicaLogging.h"

AppearancePage::AppearancePage(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::AppearancePage)
{
    ui->setupUi(this);
    ui->currentFontLabel->setText(qApp->font().family());

	//! default
	this->m_PreviousFont = qApp->font();
	this->m_PreviousStyleSheet = qApp->styleSheet();
	this->m_PreviousTheme = ThemeType::Light;

    //connect signals and slots
    connect(ui->selectFontBtn,SIGNAL(clicked()),this,SLOT(OnSelectFontButtonClicked()));
	connect(ui->themeComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(OnChangeTheme(int)));
}

AppearancePage::~AppearancePage()
{
    delete ui;
}

void AppearancePage::OnChangeTheme(int theme)
{
	this->m_CurrentTheme = ThemeType(theme);
	if (this->m_CurrentTheme == ThemeType::Dark)
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
			this->m_CurrentStyleSheet = f.readAll();
		}
		else
		{
			cbica::Logging(loggerFile, "file reading failed, qss file path = " + std::string(qssPath.c_str()));
		}
		f.close();
	}
	else if (this->m_CurrentTheme == ThemeType::Light)
	{
		cbica::Logging(loggerFile, "applying LIGHT theme");
		this->m_CurrentStyleSheet = "";
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
		qDebug() << " previous font = " << this->m_PreviousFont << endl;
		ui->currentFontLabel->setText(this->m_PreviousFont.family());
		this->m_SelectedFont = this->m_PreviousFont;
	}
}

void AppearancePage::OnOkay()
{
	//! update previous
	this->m_PreviousFont = this->m_SelectedFont;
	this->m_PreviousStyleSheet = this->m_CurrentStyleSheet;
	this->m_PreviousTheme = this->m_CurrentTheme;

	qDebug() << " previous font = " << this->m_PreviousFont << endl;
	qDebug() << " previous style = " << this->m_PreviousStyleSheet << endl;

	//! apply font and stylesheet
	qApp->setFont(this->m_SelectedFont);
	qApp->setStyleSheet(this->m_CurrentStyleSheet);
}

void AppearancePage::OnCancel()
{
	//! revert all changes
	ui->currentFontLabel->setText(this->m_PreviousFont.family());
	this->m_SelectedFont = this->m_PreviousFont;
	ui->themeComboBox->setCurrentIndex(this->m_PreviousTheme);

	qDebug() << " previous font = " << this->m_PreviousFont << endl;
	qDebug() << " previous style = " << this->m_PreviousStyleSheet << endl;

	qApp->setFont(this->m_PreviousFont);
	qApp->setStyleSheet(this->m_PreviousStyleSheet);
}