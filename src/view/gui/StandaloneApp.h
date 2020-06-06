#ifndef StandaloneApp_H
#define StandaloneApp_H

#include <QObject>
#include <QMutex>
#include <QThread>

#include "fAppDownloadDialog.h"
#include "yaml-cpp/node/node.h"

class StandaloneApp : public QObject
{
	Q_OBJECT

public:
	//! constructor/desctrucor
	StandaloneApp() = default;
	~StandaloneApp() = default;

	//! setters/getters
	void SetName(QString appName);
	QString GetName() const;

	std::string getStandaloneApp(QString appName);

private:

	Q_DISABLE_COPY(StandaloneApp)

	//! ivars
	// static StandaloneApp* m_Instance;
	// static QMutex m_Mutex;

	fAppDownloadDialog appDownloadDialog;
	QString m_AppName;
	
	YAML::Node m_appDownloadConfigs = YAML::LoadFile(getCaPTkDataDir() + "/appsDownloadConfigs.yaml");

	void appDownload();
	void startDownload();
	void cancelDownload();
	void startUnzip(QString fullPath, QString extractPath);
	void doneUnzip();
};

class ASyncExtract : public QThread
{
	Q_OBJECT
	void run() override {
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("true").toString());

		if (QFile::exists(this->fullPath))
		{
			QZipReader zr(this->fullPath);
			bool ret = zr.extractAll(this->extractPath);

			if(ret)
			{
				ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
				//after extraction remove the zip
				bool successfullyremoved = QFile::remove(this->fullPath.toStdString().c_str());
			}
			else
			{
				ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
			}
		}
		else
		{
			ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("false").toString());
		}

		qDebug() << "Extraction done in background" << this->fullPath << endl;

		//serialize only once
		ApplicationPreferences::GetInstance()->SerializePreferences();

		emit resultReady(this->appName);
	}	

public:
	ASyncExtract() = default;
	~ASyncExtract() = default;

	void setFullPath(QString fullPath) {
		this->fullPath = fullPath;
	}

	void setExtractPath(QString extractPath) {
		this->extractPath = extractPath;
	}

	void setAppName(QString appName) {
		this->appName = appName;
	}

private:
	Q_DISABLE_COPY(ASyncExtract)

	QString fullPath;
	QString extractPath;
	QString appName;

signals:
    void resultReady(QString appName);
};

#endif 
