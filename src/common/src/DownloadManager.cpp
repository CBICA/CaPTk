#include "DownloadManager.h"
#include <QTextStream>
#include <QDebug>
#include "CaPTkDefines.h"
#include "cbicaLogging.h"

using namespace std;

DownloadManager::DownloadManager(QObject *parent)
	: QObject(parent)
{
	this->m_downloadInProgress = false;
}

void DownloadManager::append(const QStringList &urls)
{
	for (const QString &urlAsString : urls)
		append(QUrl::fromEncoded(urlAsString.toLocal8Bit()));

	if (downloadQueue.isEmpty() && !this->m_downloadInProgress)
		QTimer::singleShot(0, this, &DownloadManager::finished);
}

void DownloadManager::append(const QUrl &url)
{
	if (downloadQueue.isEmpty() && !this->isDownloadInProgress())
		QTimer::singleShot(0, this, &DownloadManager::startNextDownload);

	downloadQueue.enqueue(url);
	++totalCount;
}

QString DownloadManager::saveFileName(const QUrl &url)
{
	QString path = url.path();
	QString basename = QFileInfo(path).fileName();

	return basename;
}

void DownloadManager::setFilename(const QString name)
{
	this->filesavename.enqueue(name);
}

void DownloadManager::startNextDownload()
{
	if (downloadQueue.isEmpty()) 
	{
		emit finished();
		return;
	}

	QUrl url = downloadQueue.dequeue();

	this->fname = this->filesavename.dequeue();
	output.setFileName(this->fname);
	if (!output.open(QIODevice::WriteOnly)) 
	{
		//log
		QString errormsg = QString("Problem opening save file '%1' for download '%2': %3\n").arg(this->fname).arg(url.toEncoded().constData()).arg(output.errorString());
		cbica::Logging(loggerFile, errormsg.toStdString());

		startNextDownload();
		return;                 // skip this download
	}

	this->m_downloadInProgress = true;
	QNetworkRequest request(url);
	currentDownload = manager.get(request);
	connect(currentDownload, &QNetworkReply::downloadProgress,
		this, &DownloadManager::downloadProgress);
	connect(currentDownload, &QNetworkReply::finished,
		this, &DownloadManager::downloadFinished);
	connect(currentDownload, &QNetworkReply::readyRead,
		this, &DownloadManager::downloadReadyRead);

	// log
	cbica::Logging(loggerFile, " Starting Download: " + std::string(url.toEncoded().constData()));
	downloadTimer.start();
}

void DownloadManager::downloadProgress(qint64 bytesReceived, qint64 bytesTotal)
{
	QString basename = QFileInfo(this->fname).fileName();

	// calculate the download speed
	double speed = bytesReceived * 1000.0 / downloadTimer.elapsed();
	QString unit;
	if (speed < 1024) {
		unit = "bytes/sec";
	}
	else if (speed < 1024 * 1024) {
		speed /= 1024;
		unit = "kB/s";
	}
	else {
		speed /= 1024 * 1024;
		unit = "MB/s";
	}

	//download message: Downloading: samplefile.zip at speed 2.1MB/s
	QString downloadspeed(QString::fromLatin1("%1 %2").arg(speed, 3, 'f', 1).arg(unit));
	QString downloadmsg = QString("Downloading: %1 at ").arg(basename) + downloadspeed;

	emit progress(bytesReceived, downloadmsg.toStdString(),bytesTotal);
}

void DownloadManager::downloadFinished()
{
	output.close();
	this->m_downloadInProgress = false;
	if (currentDownload->error()) 
	{
		// download failed - log
		QString errormsg = QString("Download Failed: %1").arg(currentDownload->errorString());
		cbica::Logging(loggerFile, errormsg.toStdString());

		//failed msg on GUI
		emit progress(0, "Download Failed!", 100);

		//remove failed output
		output.remove();
	}
	else {
		// let's check if it was actually a redirect
		if (isHttpRedirect()) 
		{
			reportRedirect();
			output.remove();
		}
		else 
		{
			cbica::Logging(loggerFile, "Download succeeded.");
			++downloadedCount;
		}
	}

	//reset progress to zero after finishing
	emit progress(0,"Download Finished!",100);

	currentDownload->deleteLater();

	//start next
	startNextDownload();
}

void DownloadManager::downloadReadyRead()
{
	output.write(currentDownload->readAll());
}

bool DownloadManager::isHttpRedirect() const
{
	int statusCode = currentDownload->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
	return statusCode == 301 || statusCode == 302 || statusCode == 303
		|| statusCode == 305 || statusCode == 307 || statusCode == 308;
}

void DownloadManager::reportRedirect()
{
	int statusCode = currentDownload->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
	QUrl requestUrl = currentDownload->request().url();

	//log
	QString requestmsg = QString("Request: %1 was redirected with code %2").arg(requestUrl.toDisplayString()).arg(statusCode);
	cbica::Logging(loggerFile, requestmsg.toStdString());

	QVariant target = currentDownload->attribute(QNetworkRequest::RedirectionTargetAttribute);
	if (!target.isValid())
		return;
	QUrl redirectUrl = target.toUrl();
	if (redirectUrl.isRelative())
		redirectUrl = requestUrl.resolved(redirectUrl);

	//log
	QString redirectmsg = QString("Redirected to: %1").arg(redirectUrl.toDisplayString());
	cbica::Logging(loggerFile, redirectmsg.toStdString());

}

bool DownloadManager::isDownloadInProgress() const
{
	return this->m_downloadInProgress;
}
