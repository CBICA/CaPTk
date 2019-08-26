#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>

class ApplicationPreferences
{
public:

    void SetFont(QString font);
    QString GetFont();

    void SetTheme(QString theme);
    QString GetTheme();

private:

    QString m_font;
    QString m_Theme;
};

#endif // APPLICATIONPREFERENCES_H
