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

	//we iterate over the map and serialize settings for each app in map
	if (!this->m_UserInstallationSettings.isEmpty())
	{
		appSettings.beginGroup("User-Installed-Applications");
		QMapIterator<QString, UserInstallationStatus> itr(this->m_UserInstallationSettings);
		while (itr.hasNext())
		{
			itr.next();
			appSettings.beginGroup(itr.key());
			UserInstallationStatus appInstallsettings = itr.value();
			appSettings.setValue("DownloadStarted", appInstallsettings.DownloadStarted);
			appSettings.setValue("DownloadFinished", appInstallsettings.DownloadFinished);
			appSettings.setValue("ExtractionStarted", appInstallsettings.ExtractionStarted);
			appSettings.setValue("ExtractionFinished", appInstallsettings.ExtractionFinished);
			appSettings.endGroup();
		}
		appSettings.endGroup();
	}
}

void ApplicationPreferences::DeSerializePreferences()
{
    QSettings appSettings(QSettings::IniFormat, QSettings::UserScope,
        "UPenn", "CaPTk");
	QString filename = appSettings.fileName();
    std::string fname = filename.toStdString();
	if (QFile(filename).exists())
	{
		this->SetUserPreferencesAvailability(QVariant(true).toString());
		appSettings.beginGroup("Appearance");
		this->SetFont(appSettings.value("Font").toString());
		this->SetTheme(appSettings.value("Theme").toString());
		appSettings.endGroup();

		//we enter the 'user-installed-applications' group
		//then we iterate over all the groups and  de-serialize settings for each app in group
		appSettings.beginGroup("User-Installed-Applications");
		QStringList groups = appSettings.childGroups();
		foreach(QString group, groups)
		{
			qDebug() << " group: " << group << endl;
			appSettings.beginGroup(group);
			UserInstallationStatus appInstallsettings;
			QStringList keys = appSettings.childKeys();
			QStringListIterator itr(keys);
			while (itr.hasNext())
			{
				QString key = itr.next();
				QString value = appSettings.value(key).toString();
				qDebug() << "key: " << key << " value: " << value << endl;

				if (!key.compare("DownloadStarted", Qt::CaseSensitivity::CaseInsensitive))
					appInstallsettings.DownloadStarted = value;
				if (!key.compare("DownloadFinished", Qt::CaseSensitivity::CaseInsensitive))
					appInstallsettings.DownloadFinished = value;
				if (!key.compare("ExtractionStarted", Qt::CaseSensitivity::CaseInsensitive))
					appInstallsettings.ExtractionStarted = value;
				if (!key.compare("ExtractionFinished", Qt::CaseSensitivity::CaseInsensitive))
					appInstallsettings.ExtractionFinished = value;

			}
			this->m_UserInstallationSettings[group] = appInstallsettings;
			appSettings.endGroup();
		}
	}
	else
		this->SetUserPreferencesAvailability(QVariant(false).toString());
}

void ApplicationPreferences::DisplayPreferences()
{
	qDebug() << " ApplicationPreferences::DisplayPreferences() " << endl;
	qDebug() << " font = " << this->m_Font << endl;
	qDebug() << " theme = " << this->m_Theme << endl;

	//we iterate over the map and write out settings for each app in map
	if (!this->m_UserInstallationSettings.isEmpty())
	{
		//appSettings.beginGroup("User-Installed-Applications");
		QMapIterator<QString, UserInstallationStatus> itr(this->m_UserInstallationSettings);
		while (itr.hasNext())
		{
			itr.next();
			qDebug() << " group: " << itr.key();
			UserInstallationStatus appInstallsettings = itr.value();

			qDebug() << "DownloadStarted " << appInstallsettings.DownloadStarted;
			qDebug() << "DownloadFinished " << appInstallsettings.DownloadFinished;
			qDebug() << "ExtractionStarted " << appInstallsettings.ExtractionStarted;
			qDebug() << "ExtractionFinished " << appInstallsettings.ExtractionFinished;

		}
	}

}

void ApplicationPreferences::AddApplication(QString app)
{
	if (!this->m_UserInstallationSettings.contains(app))
	{
		UserInstallationStatus appInstallSettings;
		this->m_UserInstallationSettings.insert(app, appInstallSettings);
	}
}

ApplicationPreferences::UserInstallationStatus ApplicationPreferences::GetUserInstallationSettingsForApp(QString app)
{
	if (this->m_UserInstallationSettings.contains(app))
	{
		return this->m_UserInstallationSettings[app];
	}
}

void ApplicationPreferences::SetUserPreferencesAvailability(QString available)
{
	this->m_UserPreferencesAvailability = available;
}

QString ApplicationPreferences::GetUserPreferencesAvailability() const
{
    return m_UserPreferencesAvailability;
}

void ApplicationPreferences::SetDownloadStartedStatus(QString app, QString status)
{
	//this->m_LibraDownloadStarted = status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		this->m_UserInstallationSettings[app].DownloadStarted = status;
	}
}

QString ApplicationPreferences::GetDownloadStartedStatus(QString app)
{
	//return this->m_LibraDownloadStarted;
	QString status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		status = this->m_UserInstallationSettings[app].DownloadStarted;
	}
	return status;
}

void ApplicationPreferences::SetDownloadFinishedStatus(QString app, QString status)
{
	//this->m_LibraDownloadFinished = status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		this->m_UserInstallationSettings[app].DownloadFinished = status;
	}
}

QString ApplicationPreferences::GetDownloadFinishedStatus(QString app)
{
	//return this->m_LibraDownloadFinished;
	QString status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		status = this->m_UserInstallationSettings[app].DownloadFinished;
	}
	return status;
}

void ApplicationPreferences::SetExtractionStartedStatus(QString app, QString status)
{
	//this->m_LibraExtractionStarted = status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		this->m_UserInstallationSettings[app].ExtractionStarted = status;
	}
}

QString ApplicationPreferences::GetExtractionStartedStatus(QString app)
{
	//return this->m_LibraExtractionStarted;
	QString status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		status = this->m_UserInstallationSettings[app].ExtractionStarted;
	}
	return status;
}

void ApplicationPreferences::SetExtractionFinishedStatus(QString app, QString status)
{
	//this->m_LibraExtractionFinished = status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		this->m_UserInstallationSettings[app].ExtractionFinished = status;
	}
}

QString ApplicationPreferences::GetExtractionFinishedStatus(QString app)
{
	//return this->m_LibraExtractionFinished;
	QString status;
	if (this->m_UserInstallationSettings.contains(app))
	{
		status = this->m_UserInstallationSettings[app].ExtractionFinished;
	}
	return status;
}