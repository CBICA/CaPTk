#include "ApplicationPreferences.h"
#include <QSettings>
#include <QDebug>

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
	QSettings appSettings(QSettings::IniFormat,QSettings::SystemScope,
		"UPenn", "CaPTk");
	appSettings.beginGroup("Appearance");
	appSettings.setValue("Font", this->m_Font);
	appSettings.setValue("Theme", this->m_Theme);
	appSettings.endGroup();
	qDebug() << " filename = " << appSettings.fileName() << endl;
}

void ApplicationPreferences::DisplayPreferences()
{
	qDebug() << " font = " << this->m_Font << endl;
	qDebug() << " theme = " << this->m_Theme << endl;
}