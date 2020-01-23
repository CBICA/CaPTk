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
    this->m_PreviousTheme = this->m_SelectedTheme = ThemeType::Dark; // make sure the default matches what we initialize to

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
    //! Apply
    OnChangeTheme(m_SelectedTheme);
    RefreshAppToSelectedAppearance();

}

AppearancePage::~AppearancePage()
{
    delete ui;
}

void AppearancePage::OnChangeTheme(int theme)
{
    // OnChangeTheme updates the instance variables to contain the selected theme and style sheet but does not apply it.
    m_PreviousTheme = m_SelectedTheme;
    m_PreviousStyleSheet = m_SelectedStyleSheet;
	this->m_SelectedTheme = ThemeType(theme);
	if (this->m_SelectedTheme == ThemeType::Dark)
	{
		cbica::Logging(loggerFile, "applying DARK theme");
		std::string qssPath;
#ifdef QT_DEBUG
		qssPath = getCaPTkDataDir() + /*"../etc/"*/ + "/captk.qss";
#else
        qssPath = getCaPTkDataDir() + /*"../etc/"*/ + CAPTK_STYLESHEET_FILE;
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
    qApp->setFont(m_SelectedFont);
    OnChangeTheme(m_SelectedTheme);
    qApp->setStyleSheet(m_SelectedStyleSheet);

    //! Bring up a confirmation box and query user to keep or discard changes

    auto response = RaiseConfirmationDialog();

    if(response == QMessageBox::Apply) { //! Configuration only sets permanently if the user explicitly confirms
        FinalizeSelectedChanges();
    }
    else {
        //! re-apply the previous font and stylesheet to the application
        UndoSelectedChanges();
        RefreshAppToSelectedAppearance();
    }

}

void AppearancePage::OnCancel()
{
	//! revert all changes
    UndoSelectedChanges();
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


void AppearancePage::UndoSelectedChanges() {
    //! reset ui elements and instance variables to reflect previous state
    ui->currentFontLabel->setText(this->m_PreviousFont.family());
    this->m_SelectedFont = this->m_PreviousFont;
    this->m_SelectedTheme = this->m_PreviousTheme;
    this->m_SelectedStyleSheet = this->m_PreviousStyleSheet;
    ui->themeComboBox->setCurrentIndex(this->m_PreviousTheme);
}
void AppearancePage::FinalizeSelectedChanges() {
    //! set ApplicationPreferences to the desired settings
    ApplicationPreferences::GetInstance()->SetFont(this->m_SelectedFont.toString());
    ApplicationPreferences::GetInstance()->SetTheme(QVariant::fromValue(this->m_SelectedTheme).toString());
    //! update previous settings to reflect the change
    this->m_PreviousFont = this->m_SelectedFont;
    this->m_PreviousStyleSheet = this->m_SelectedStyleSheet;
    this->m_PreviousTheme = this->m_SelectedTheme;
}

void AppearancePage::RefreshAppToSelectedAppearance() {
    qApp->setFont(this->m_SelectedFont);
    qApp->setStyleSheet(this->m_SelectedStyleSheet);
}

int AppearancePage::RaiseConfirmationDialog() {

    //! Bring up a confirmation box and query user to keep or discard changes
    QMessageBox* confirmationBox = new QMessageBox(this);
    QFont defaultFont;
    QString systemDefault = defaultFont.defaultFamily();
    // Hard-coded style sheet for this message box so that it is guaranteed to be readable regardless of the chosen settings
    QString tempStyleSheet = "* { color : black; background-color : white; font-size : 16px;"
                             " font : " + systemDefault + " }";
    confirmationBox->setText("You are seeing a preview of your selected changes. Do you wish to keep them?");
    confirmationBox->setInformativeText("If you do not confirm within 10 seconds, these changes will automatically revert.");
    confirmationBox->setStyleSheet(tempStyleSheet);
    confirmationBox->setStandardButtons(QMessageBox::Apply | QMessageBox::Cancel);
    confirmationBox->setIcon(QMessageBox::NoIcon);

    //! wait some time and then return the response value (integer)
    int timeout_in_milliseconds = 10000; // set a reasonable timeout before we auto-cancel
    QTimer* confirmationTimer = new QTimer(this);
    confirmationTimer->setSingleShot(true);
    connect(confirmationTimer, SIGNAL(timeout()), confirmationBox, SLOT(reject()));
    confirmationTimer->start(timeout_in_milliseconds);
    confirmationBox->exec();
    confirmationTimer->stop();
    auto response = confirmationBox->result(); // grab the response value
    return response;

}
