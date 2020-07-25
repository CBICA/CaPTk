#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>
#include <QMutex>
#include <QVariant>

class ApplicationPreferences : public QObject
{
public:
	struct UserInstallation
	{
		//app download extract settings
		QString DownloadStarted = "false";
		QString DownloadFinished = "false";
		QString ExtractionStarted = "false";
		QString ExtractionFinished = "false";
	};

	static ApplicationPreferences* GetInstance();

	//! setters/getters
	void SetFont(QString font);
	QString GetFont() const;

	void SetTheme(QString theme);
	QString GetTheme() const;

	void SetUserPreferencesAvailability(QString available);
	QString GetUserPreferencesAvailability() const;

	//void SetLibraDownloadStartedStatus(QString status);
	//QString GetLibraDownloadStartedStatus();

	//void SetLibraDownloadFinishedStatus(QString status);
	//QString GetLibraDownloadFinishedStatus();

	//void SetLibraExtractionStartedStatus(QString status);
	//QString GetLibraExtractionStartedStatus();

	//void SetLibraExtractionFinishedStatus(QString status);
	//QString GetLibraExtractionFinishedStatus();

    //! Serialize
	void SerializePreferences();

	//! De-Serliaze
	void DeSerializePreferences();

	//! print preferences(for debugging purposes)
	void DisplayPreferences();

private:
	//! constructor/desctrucor
	ApplicationPreferences() = default;
	~ApplicationPreferences() = default;

	Q_DISABLE_COPY(ApplicationPreferences)

	//! ivars
	static ApplicationPreferences* m_Instance;
	static QMutex m_Mutex;

	QString m_Font;
	QString m_Theme;
	QString m_UserPreferencesAvailability = QVariant("false").toString();

	//app download extract settings
	QString m_LibraDownloadStarted;
	QString m_LibraDownloadFinished;
	QString m_LibraExtractionStarted;
	QString m_LibraExtractionFinished;
};

#endif // APPLICATIONPREFERENCES_H