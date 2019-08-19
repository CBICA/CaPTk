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
	ThemeType t = ThemeType(theme);
	if (t == ThemeType::Dark)
	{
		qDebug() << t << endl;
		std::string qssPath;
#ifdef QT_DEBUG
		qssPath = getCaPTkDataDir() + /*"../etc/"*/ + "/captk.qss";
#else
		qssPath = getCaPTkDataDir() + "../etc/" + CAPTK_STYLESHEET_FILE;
#endif
		QFile f(qssPath.c_str());
		bool success = f.open(QFile::ReadOnly);
		qDebug() << "success = " << success << endl;
		if (success)
		{
			qDebug() << "reading succeeded, applying DARK theme" << endl;
			cbica::Logging(loggerFile, "qss file path = " + std::string(qssPath.c_str()));
			cbica::Logging(loggerFile, "qss file read succeeded, applying DARK theme");
			this->m_currentStyleSheet = f.readAll();
			//qApp->setStyleSheet(f.readAll());
		}
		else
		{
			qDebug() << "reading failed" << endl;
			cbica::Logging(loggerFile, "qss file path = " + std::string(qssPath.c_str()));
			cbica::Logging(loggerFile, "qss file read failed");
		}
		f.close();
	}
	else if (t == ThemeType::Light)
	{
		qDebug() << t << endl;
		cbica::Logging(loggerFile, "applying LIGHT theme");
		this->m_currentStyleSheet = "";
		//qApp->setStyleSheet("");
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
        qDebug() << "font           : " << font;
        qDebug() << "font weight    : " << font.weight();
        qDebug() << "font family    : " << font.family();
        qDebug() << "font style     : " << font.style();  //  StyleNormal = 0, StyleItalic = 1, StyleOblique = 2
        qDebug() << "font pointSize : " << font.pointSize();
		ui->currentFontLabel->setText(font.family());
		this->m_selectedFont = font;
    }
}


void AppearancePage::OnOkay()
{
	qDebug() << "AppearancePage::OnOkay()" << endl;
	qApp->setFont(this->m_selectedFont);
	qDebug() << "selected font = " << this->m_selectedFont << endl;
	qApp->setStyleSheet(this->m_currentStyleSheet);
}

void AppearancePage::OnCancel()
{
	qDebug() << "AppearancePage::OnCancel()" << endl;
}