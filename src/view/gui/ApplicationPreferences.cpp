#include "ApplicationPreferences.h"
#include <QSettings>
#include <QFile>
#include <QDebug>
#include "cbicaLogging.h"
#include <CaPTkDefines.h>

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

    qDebug() << " status = " << appSettings.status();
    cbica::Logging(loggerFile, "ApplicationPreferences::SerializePreferences status: " + QVariant::fromValue(appSettings.status()).toString().toStdString() );
}

void ApplicationPreferences::DeSerializePreferences()
{
    QSettings appSettings(QSettings::IniFormat, QSettings::UserScope,
		"UPenn", "CaPTk");
	QString filename = appSettings.fileName();
    std::string fname = filename.toStdString();
	if (QFile(filename).exists())
	{
		this->SetFileAvailability(QVariant(true).toString());
		appSettings.beginGroup("Appearance");
		this->SetFont(appSettings.value("Font").toString());
		this->SetTheme(appSettings.value("Theme").toString());
		appSettings.endGroup();
	}
	else
		this->SetFileAvailability(QVariant(false).toString());
}

void ApplicationPreferences::DisplayPreferences()
{
	qDebug() << " ApplicationPreferences::DisplayPreferences() " << endl;
	qDebug() << " font = " << this->m_Font << endl;
	qDebug() << " theme = " << this->m_Theme << endl;
}

void ApplicationPreferences::SetFileAvailability(QString available)
{
	this->m_FileAvailability = available;
}

QString ApplicationPreferences::GetFileAvailability() const
{
	return m_FileAvailability;
}
