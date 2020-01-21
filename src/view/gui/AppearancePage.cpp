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
#include <QTimer>

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
        ApplicationPreferences::GetInstance()->SetFont(this->m_SelectedFont.toString());
        ApplicationPreferences::GetInstance()->SetTheme(QVariant::fromValue(this->m_SelectedTheme).toString());
        ApplicationPreferences::GetInstance()->DisplayPreferences();
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
    //! apply selected font and stylesheet for preview
	qApp->setFont(this->m_SelectedFont);
	qApp->setStyleSheet(this->m_SelectedStyleSheet);

    //! Bring up a confirmation box and query user to keep or discard changes
    QMessageBox* confirmationBox = new QMessageBox(this);
    confirmationBox->setText("You are seeing a preview of your selected changes.");
    confirmationBox->setInformativeText("If you do not confirm within 5 seconds, these changes will automatically revert.");
    confirmationBox->setStandardButtons(QMessageBox::Apply | QMessageBox::Discard);
    confirmationBox->setIcon(QMessageBox::NoIcon);
    int timeout_in_milliseconds = 5000; // set a reasonable timeout before we auto-cancel
    QTimer* confirmationTimer = new QTimer(this);
    confirmationTimer->setSingleShot(true);
    connect(confirmationTimer, SIGNAL(timeout()), confirmationBox, SLOT(reject()));
    confirmationTimer->start(timeout_in_milliseconds);
    confirmationBox->exec();
    confirmationTimer->stop();
    auto response = confirmationBox->result();

    if(response == QMessageBox::Apply) { // If anything else happens, discard changes
        ApplicationPreferences::GetInstance()->SetFont(this->m_SelectedFont.toString());
        ApplicationPreferences::GetInstance()->SetTheme(QVariant::fromValue(this->m_SelectedTheme).toString());
        ApplicationPreferences::GetInstance()->DisplayPreferences();
        //! update previous settings to reflect the change
        this->m_PreviousFont = this->m_SelectedFont;
        this->m_PreviousStyleSheet = this->m_SelectedStyleSheet;
        this->m_PreviousTheme = this->m_SelectedTheme;
    }
    else {
        //! re-apply the previous font and stylesheet to the application
        OnCancel();
    }

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
    ApplicationPreferences::GetInstance()->DisplayPreferences();
    //this->OnOkay();
}
