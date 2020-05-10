#include "StandaloneApps.h"

#include <QSettings>
#include <QFile>
#include <QDebug>
#include <CaPTkDefines.h>
#include "cbicaLogging.h"
#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "fAppDownloadDialog.h"

StandaloneApps* StandaloneApps::m_Instance = nullptr;
QMutex StandaloneApps::m_Mutex;

StandaloneApps* StandaloneApps::GetInstance()
{
	if (m_Instance == nullptr)
	{
		m_Mutex.lock();
		m_Instance = new StandaloneApps();
      // m_Instance.SetName(appName);
		m_Mutex.unlock();
	}
	return m_Instance;
}

// void StandaloneApps::SetName(QString appName)
// {
//     this->m_AppName = appName;
// }

// QString StandaloneApps::GetName() const
// {
//     return m_AppName;
// }

void StandaloneApps::SetAction(QString action)
{
    this->m_Action = action;
}

QString StandaloneApps::GetAction() const
{
    return m_Action;
}

void StandaloneApps::SetStatus(QString status)
{
    this->m_Status = status;
}

QString StandaloneApps::GetStatus() const
{
    return m_Status;
}

void StandaloneApps::StoreAppSetting(QString action, QString status, QString appName)
{
    QSettings appSettings(QSettings::IniFormat,QSettings::UserScope,
		"UPenn", "CaPTk");

	appSettings.beginGroup(appName);
	appSettings.setValue("Action", action);
	appSettings.setValue("Status", status);
	appSettings.endGroup();

    qDebug() << " status = " << appSettings.status();
    cbica::Logging(loggerFile, "StandaloneApps::SetStatus status: " + QVariant::fromValue(appSettings.status()).toString().toStdString() );
}

void StandaloneApps::RetreiveAppSetting(QString appName)
{
    QSettings appSettings(QSettings::IniFormat, QSettings::UserScope,
		"UPenn", "CaPTk");

	QString filename = appSettings.fileName();
	std::string fname = filename.toStdString();

	if (QFile(filename).exists())
	{
		this->SetFileAvailability(QVariant(true).toString());

		appSettings.beginGroup(appName);
		this->SetAction(appSettings.value("Action").toString());
		this->SetStatus(appSettings.value("Status").toString());
		appSettings.endGroup();
	}
	else
		this->SetFileAvailability(QVariant(false).toString());
}

void StandaloneApps::Debug()
{
	qDebug() << " StandaloneApps::Debug() " << endl;
	qDebug() << " action = " << this->m_Action << endl;
	qDebug() << " status = " << this->m_Status << endl;
}

void StandaloneApps::SetFileAvailability(QString available)
{
	this->m_FileAvailability = available;
}

QString StandaloneApps::GetFileAvailability() const
{
	return m_FileAvailability;
}
