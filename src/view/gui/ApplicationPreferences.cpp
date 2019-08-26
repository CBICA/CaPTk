#include "ApplicationPreferences.h"

void AppSettings::SetFont(QString font)
{
    this->m_font = font;
}

QString AppSettings::GetFont()
{
    return m_font;
}

void AppSettings::SetTheme(QString theme)
{
    this->m_Theme = theme;
}

QString AppSettings::GetTheme()
{
    return m_Theme;
}

