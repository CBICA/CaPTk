#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>

class ApplicationPreferences : public QObject
{
public:

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
    QString m_Font;
    QString m_Theme;
};

#endif // APPLICATIONPREFERENCES_H
