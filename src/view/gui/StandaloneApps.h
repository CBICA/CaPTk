#ifndef STANDALONEAPPS_H
#define STANDALONEAPPS_H

#include <QObject>
#include <QMutex>
#include <QVariant>

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