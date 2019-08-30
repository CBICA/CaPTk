#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>
#include <QMutex>
#include <QVariant>

class ApplicationPreferences : public QObject
{
public:
	static ApplicationPreferences* GetInstance();

	//! setters/getters
    void SetFont(QString font);
    QString GetFont() const;

    void SetTheme(QString theme);
    QString GetTheme() const;

	void SetFileAvailability(QString available);
	QString GetFileAvailability() const;

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
	QString m_FileAvailability = QVariant("false").toString();
};

#endif // APPLICATIONPREFERENCES_H
