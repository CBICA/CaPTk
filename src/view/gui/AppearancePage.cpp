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
#include <QMessageBox>
#include <QTimer>
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

    // default
	this->m_PreviousFont = this->m_SelectedFont = qApp->font();
	this->m_PreviousStyleSheet = this->m_SelectedStyleSheet = qApp->styleSheet();
	this->m_PreviousTheme = this->m_SelectedTheme = ThemeType::Light;

    // connect signals and slots
    connect(ui->selectFontBtn,SIGNAL(clicked()),this,SLOT(OnSelectFontButtonClicked()));
	connect(ui->themeComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(OnChangeTheme(int)));

    // initialize with DARK theme
	if (!QVariant(ApplicationPreferences::GetInstance()->GetUserPreferencesAvailability()).toBool())
	{
		ui->themeComboBox->setCurrentIndex(ThemeType::Dark);
        this->ApplySelectedAppearance();

		//update appearance preferences in application preferences
		this->UpdateAppearancePreferences();
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
    // apply selected font and stylesheet
    ApplySelectedAppearance();

    if (GetConfirmationFromUser()) {

        // update appearance preferences to be serialized later
        UpdateAppearancePreferences();

        // update previous
        this->m_PreviousFont = this->m_SelectedFont;
        this->m_PreviousStyleSheet = this->m_SelectedStyleSheet;
        this->m_PreviousTheme = this->m_SelectedTheme;
    } else {
        // if the user doesn't confirm, revert all changes
        OnCancel();
    }

}



void AppearancePage::OnCancel()
{
    // revert all changes
    ui->currentFontLabel->setText(this->m_PreviousFont.family());
    this->m_SelectedFont = this->m_PreviousFont;
    this->m_SelectedTheme = this->m_PreviousTheme;
    ui->themeComboBox->setCurrentIndex(this->m_PreviousTheme);

    // refresh appearance to match
    ApplySelectedAppearance();
}

void AppearancePage::Restore()
{
    // we restore only if user preferences are available
    if (QVariant(ApplicationPreferences::GetInstance()->GetUserPreferencesAvailability()).toBool())
    {
        // get font/theme from preferences, update the appearance-page selections to reflect this
        this->m_SelectedFont.fromString(ApplicationPreferences::GetInstance()->GetFont());
        ui->currentFontLabel->setText(this->m_SelectedFont.family());

        auto metaEnum = QMetaEnum::fromType<ThemeType>();
        this->m_SelectedTheme = static_cast<ThemeType>(metaEnum.keyToValue(
                                                           ApplicationPreferences::GetInstance()->GetTheme().toStdString().c_str()));

        ui->themeComboBox->setCurrentIndex(this->m_SelectedTheme);

        // now apply the font/theme
        this->ApplySelectedAppearance();
    }
}

bool AppearancePage::GetConfirmationFromUser()
{
    // Bring up a confirmation box and query user to keep or discard changes
    QMessageBox* confirmationBox = new QMessageBox(this);
    QFont defaultFont;
    QString systemDefault = defaultFont.defaultFamily(); // try to grab a system default font

    // Hard-coded style sheet for this message box so that it is guaranteed to be readable regardless of the chosen settings
    QString tempStyleSheet = "* { color : black; background-color : white; font-size : 16px;"
                             " font : " + systemDefault + " }";

    confirmationBox->setText("You are seeing a preview of your selected changes. Do you wish to keep them?");
    confirmationBox->setInformativeText("If you do not confirm within 10 seconds, these changes will automatically revert.");
    confirmationBox->setStyleSheet(tempStyleSheet);
    confirmationBox->setStandardButtons(QMessageBox::Apply | QMessageBox::Cancel);
    confirmationBox->setIcon(QMessageBox::NoIcon);

    // wait some time and then return the response value (integer)
    int timeout_in_milliseconds = 10000; // set a reasonable timeout before we auto-cancel
    QTimer* confirmationTimer = new QTimer(this);
    confirmationTimer->setSingleShot(true);
    connect(confirmationTimer, SIGNAL(timeout()), confirmationBox, SLOT(reject()));
    confirmationTimer->start(timeout_in_milliseconds);
    confirmationBox->exec();
    confirmationTimer->stop();
    auto response = confirmationBox->result(); // grab the response value

    delete confirmationTimer; // clean up timer and box
    delete confirmationBox;

    if (response == QMessageBox::Apply) {
        return true;
    }
    else {
        return false;
    }
}

void AppearancePage::ApplySelectedAppearance()
{
    qApp->setFont(m_SelectedFont);
    qApp->setStyleSheet(m_SelectedStyleSheet);
}

void AppearancePage::UpdateAppearancePreferences() {
    ApplicationPreferences::GetInstance()->SetFont(this->m_SelectedFont.toString());
    ApplicationPreferences::GetInstance()->SetTheme(QVariant::fromValue(this->m_SelectedTheme).toString());
}
