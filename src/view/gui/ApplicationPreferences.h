#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>
#include <QMutex>

class ApplicationPreferences : public QObject
{
public:
	static ApplicationPreferences* GetInstance();

	//! setters/getters
    void SetFont(QString font);
    QString GetFont() const;

    void SetTheme(QString theme);
    QString GetTheme() const;

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
};

#endif // APPLICATIONPREFERENCES_H
