#ifndef APPLICATIONPREFERENCES_H
#define APPLICATIONPREFERENCES_H

#include <QObject>

class ApplicationPreferences : public QObject
{
public:

    void SetFont(QString font);
    QString GetFont() const;

    void SetTheme(QString theme);
    QString GetTheme() const;

private:

    QString m_Font;
    QString m_Theme;
};

#endif // APPLICATIONPREFERENCES_H
