#ifndef ASYNCEXTRACT_H
#define ASYNCEXTRACT_H

#include <QObject>
#include <QThread>
#include <QFile>
#include <QDebug>

#include "CaPTkUtils.h"
#include "CaPTkGUIUtils.h"
#include "QZipReader.h"
#include "ApplicationPreferences.h"

class ExtractInBackground : public QThread
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
				//after extraction remove the zip
				bool successfullyremoved = QFile::remove(this->fullPath.toStdString().c_str());
			}
		}

		qDebug() << "Extraction done in background" << this->fullPath << endl;

		//serialize only once
		ApplicationPreferences::GetInstance()->SerializePreferences();
		ApplicationPreferences::GetInstance()->DisplayPreferences();

		emit resultReady(this->appName);
	}	

public:
	ExtractInBackground() = default;
	~ExtractInBackground() = default;

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
	Q_DISABLE_COPY(ExtractInBackground)

	QString fullPath;
	QString extractPath;
	QString appName;

signals:
    void resultReady(QString appName);
};

#endif