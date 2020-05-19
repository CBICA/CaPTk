#ifndef STANDALONEAPPS_H
#define STANDALONEAPPS_H

#include <QObject>
#include <QMutex>
#include <QVariant>
#include <QThread>

#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "QZipReader.h"

class StandaloneApps : public QObject
{
public:
	static StandaloneApps* GetInstance();

	//! setters/getters
	// void SetName(QString appName);
	// QString GetName() const;

	void SetAction(QString action);
	QString GetAction() const;

	void SetStatus(QString status);
	QString GetStatus() const;

	void SetFileAvailability(QString available);
	QString GetFileAvailability() const;

	void StoreAppSetting(QString action, QString status, QString appName);
	void RetreiveAppSetting(QString appName);

	//! print preferences(for debugging purposes)
	void Debug(QString step);

private:
	//! constructor/desctrucor
	StandaloneApps() = default;
	~StandaloneApps() = default;

	Q_DISABLE_COPY(StandaloneApps)

	//! ivars
	static StandaloneApps* m_Instance;
	static QMutex m_Mutex;

	// QString m_AppName;
	QString m_Action;
	QString m_Status;
	// QString m_DownloadStatus;
	// QString m_ExtractStatus;
	QString m_FileAvailability = QVariant("false").toString();
};

#endif 

class ASyncExtract : public QThread
{
	Q_OBJECT
	void run() override {
	QString result;

	StandaloneApps* stlapps = StandaloneApps::GetInstance();

	QZipReader zr(this->fullPath);
	// stlapps->StoreAppSetting("Extract", "Start", appName);

	// stlapps->RetreiveAppSetting(appName);
	// stlapps->Debug("Extraction start");

	emit resultReady(result);
	}

public:
	static ASyncExtract* GetInstance()
	{
		if (m_Instance == nullptr)
		{
			m_Mutex.lock();
			m_Instance = new ASyncExtract();
			// m_Instance.SetName(appName);
			m_Mutex.unlock();
		}
		return m_Instance;
	}
	void setFullPath(QString fullPath) {
		this->fullPath = fullPath;
	}

private:
	ASyncExtract() = default;
	~ASyncExtract() = default;
	
	static ASyncExtract* m_Instance;
	static QMutex m_Mutex;

	QString fullPath;

signals:
    void resultReady(const QString &s);
};