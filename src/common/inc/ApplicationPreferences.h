#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>
#include <QMutex>
#include <QVariant>

class ApplicationPreferences : public QObject
{
public:
	struct UserInstallationStatus
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

	void SetDownloadStartedStatus(QString app, QString status);
	QString GetDownloadStartedStatus(QString app);

	void SetDownloadFinishedStatus(QString app, QString status);
	QString GetDownloadFinishedStatus(QString app);

	void SetExtractionStartedStatus(QString app, QString status);
	QString GetExtractionStartedStatus(QString app);

	void SetExtractionFinishedStatus(QString app, QString status);
	QString GetExtractionFinishedStatus(QString app);

    //! Serialize
	void SerializePreferences();

	//! De-Serliaze
	void DeSerializePreferences();

	//! print preferences(for debugging purposes)
	void DisplayPreferences();

	//! Add a new user installed application
	void AddApplication(QString app);

	//! Get user installation settings
	UserInstallationStatus GetUserInstallationSettingsForApp(QString app);

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
	QMap<QString, UserInstallationStatus> m_UserInstallationSettings;

	//app download extract settings
	//QString m_LibraDownloadStarted;
	//QString m_LibraDownloadFinished;
	//QString m_LibraExtractionStarted;
	//QString m_LibraExtractionFinished;
};

#endif // APPLICATIONPREFERENCES_H