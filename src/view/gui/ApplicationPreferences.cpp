#include "ApplicationPreferences.h"

void ApplicationPreferences::SetFont(QString font)
{
    this->m_font = font;
}

QString ApplicationPreferences::GetFont()
{
    return m_font;
}

void ApplicationPreferences::SetTheme(QString theme)
{
    this->m_Theme = theme;
}

QString ApplicationPreferences::GetTheme()
{
    return m_Theme;
}

