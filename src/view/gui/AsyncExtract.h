#ifndef ASYNCEXTRACT_H
#define ASYNCEXTRACT_H

#include <QObject>
#include <QThread>
#include <QFile>

#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "QZipReader.h"
#include "ApplicationPreferences.h"

class ASyncExtract : public QThread
{
	Q_OBJECT
	void run() override {

		// StandaloneApp* stlapps = StandaloneApp::GetInstance();

		// QZipReader zr(this->fullPath);
		// stlapps->StoreAppSetting("Extract", "Start", this->appName);

		// stlapps->RetreiveAppSetting(this->appName);
		// stlapps->Debug("Extraction start");

		// bool extracted = zr.extractAll(this->extractPath);

		// stlapps->RetreiveAppSetting(appName);
    	// stlapps->Debug("Extraction done in background");
		ApplicationPreferences::GetInstance()->SetLibraExtractionStartedStatus(QVariant("true").toString());

		if (QFile::exists(this->fullPath))
		{
			QZipReader zr(this->fullPath);
			bool ret = zr.extractAll(this->extractPath);

			if(ret)
			{
				ApplicationPreferences::GetInstance()->SetLibraExtractionFinishedStatus(QVariant("true").toString());
				//after extraction remove the zip
				bool successfullyremoved = QFile::remove(fullPath.c_str());
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

	QString fullPath;
	QString extractPath;
	QString appName;

signals:
    void resultReady(QString appName);
};

#endif