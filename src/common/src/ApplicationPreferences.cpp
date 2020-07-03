#include "ApplicationPreferences.h"
#include <QSettings>
#include <QFile>
#include <QDebug>
//#include "cbicaLogging.h"
//#include <CaPTkDefines.h>

ApplicationPreferences* ApplicationPreferences::m_Instance = nullptr;
QMutex ApplicationPreferences::m_Mutex;

ApplicationPreferences* ApplicationPreferences::GetInstance()
{
	if (m_Instance == nullptr)
	{
		m_Mutex.lock();
		m_Instance = new ApplicationPreferences();
		m_Mutex.unlock();
	}
	return m_Instance;
}

void ApplicationPreferences::SetFont(QString font)
{
    this->m_Font = font;
}

QString ApplicationPreferences::GetFont() const
{
    return m_Font;
}

void ApplicationPreferences::SetTheme(QString theme)
{
    this->m_Theme = theme;
}

QString ApplicationPreferences::GetTheme() const
{
    return m_Theme;
}

void ApplicationPreferences::SerializePreferences()
{
    QSettings appSettings(QSettings::IniFormat,QSettings::UserScope,
        "UPenn", "CaPTk");

	appSettings.beginGroup("Appearance");
	appSettings.setValue("Font", this->m_Font);
	appSettings.setValue("Theme", this->m_Theme);
	appSettings.endGroup();

    appSettings.beginGroup("User-Installed-Applications");
    appSettings.beginGroup("Libra");
    appSettings.setValue("Libra_Download_Started",this->m_LibraDownloadStarted);
    appSettings.setValue("Libra_Download_Finished",this->m_LibraDownloadFinished);
    appSettings.setValue("Libra_Extraction_Started",this->m_LibraExtractionStarted);
    appSettings.setValue("Libra_Extraction_Finished",this->m_LibraExtractionFinished);
    appSettings.endGroup();
    appSettings.endGroup();

    // qDebug() << " status = " << appSettings.status();
    //cbica::Logging(loggerFile, "ApplicationPreferences::SerializePreferences status: " + QVariant::fromValue(appSettings.status()).toString().toStdString() );
}

void ApplicationPreferences::DeSerializePreferences()
{
    QSettings appSettings(QSettings::IniFormat, QSettings::UserScope,
        "UPenn", "UPenn");
	QString filename = appSettings.fileName();
    std::string fname = filename.toStdString();
	if (QFile(filename).exists())
	{
		this->SetUserPreferencesAvailability(QVariant(true).toString());
		appSettings.beginGroup("Appearance");
		this->SetFont(appSettings.value("Font").toString());
		this->SetTheme(appSettings.value("Theme").toString());
		appSettings.endGroup();

		appSettings.beginGroup("User-Installed-Applications");
		appSettings.beginGroup("Libra");
		this->SetLibraDownloadStartedStatus(appSettings.value("Libra_Download_Started").toString());
		this->SetLibraDownloadFinishedStatus(appSettings.value("Libra_Download_Finished").toString());
		this->SetLibraExtractionStartedStatus(appSettings.value("Libra_Extraction_Started").toString());
		this->SetLibraExtractionFinishedStatus(appSettings.value("Libra_Extraction_Finished").toString());
		appSettings.endGroup();
		appSettings.endGroup();
	}
	else
		this->SetUserPreferencesAvailability(QVariant(false).toString());
}

void ApplicationPreferences::DisplayPreferences()
{
	qDebug() << " ApplicationPreferences::DisplayPreferences() " << endl;
	qDebug() << " font = " << this->m_Font << endl;
	qDebug() << " theme = " << this->m_Theme << endl;
	qDebug() << " Libra_Download_Started = " << this->m_LibraDownloadStarted << endl;
	qDebug() << " Libra_Download_Finished = " << this->m_LibraDownloadFinished << endl;
	qDebug() << " Libra_Extraction_Started = " << this->m_LibraExtractionStarted << endl;
	qDebug() << " Libra_Extraction_Finished = " << this->m_LibraExtractionFinished << endl;

}

void ApplicationPreferences::SetUserPreferencesAvailability(QString available)
{
	this->m_UserPreferencesAvailability = available;
}

QString ApplicationPreferences::GetUserPreferencesAvailability() const
{
    return m_UserPreferencesAvailability;
}

void ApplicationPreferences::SetLibraDownloadStartedStatus(QString status)
{
	this->m_LibraDownloadStarted = status;
}

QString ApplicationPreferences::GetLibraDownloadStartedStatus()
{
	return this->m_LibraDownloadStarted;
}

void ApplicationPreferences::SetLibraDownloadFinishedStatus(QString status)
{
	this->m_LibraDownloadFinished = status;
}

QString ApplicationPreferences::GetLibraDownloadFinishedStatus()
{
	return this->m_LibraDownloadFinished;
}

void ApplicationPreferences::SetLibraExtractionStartedStatus(QString status)
{
	this->m_LibraExtractionStarted = status;
}

QString ApplicationPreferences::GetLibraExtractionStartedStatus()
{
	return this->m_LibraExtractionStarted;
}

void ApplicationPreferences::SetLibraExtractionFinishedStatus(QString status)
{
	this->m_LibraExtractionFinished = status;
}

QString ApplicationPreferences::GetLibraExtractionFinishedStatus()
{
	return this->m_LibraExtractionFinished;
}