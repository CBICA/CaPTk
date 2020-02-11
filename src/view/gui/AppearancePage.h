/**
\file  AppearancePage.h

\brief Declaration of AppearancePage class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/
#ifndef APPEARANCEPAGE_H
#define APPEARANCEPAGE_H

#include <QWidget>
#include "IPreferencePage.h"

namespace Ui {
class AppearancePage;
}

class AppearancePage : public QWidget, public IPreferencePage
{
    Q_OBJECT

public:
	enum ThemeType
	{
		Light = 0,
		Dark
	};

    explicit AppearancePage(QWidget *parent = nullptr);
    ~AppearancePage();

	Q_ENUM(ThemeType)

public slots:
	void OnSelectFontButtonClicked();
	void OnChangeTheme(int);

    /**
     * @brief OnOkay is called when the Appearance Page dialog is accepted. It prompts for user confirmation.
     */
	void OnOkay() override;

    /**
     * @brief OnCancel is called when the Appearance Page dialog is rejected. It resets the appearance to its previous configuration.
     */
	void OnCancel() override;
	void Restore() override;

private:
    Ui::AppearancePage *ui;

	//! ivars for handling cancel cases
    QFont m_SelectedFont, m_PreviousFont;
	QString m_SelectedStyleSheet, m_PreviousStyleSheet;
	ThemeType m_SelectedTheme, m_PreviousTheme;

    /**
     * @brief ApplySelectedAppearance applies whatever the currently selected font and stylesheet is to the global qApp.
     */
    void ApplySelectedAppearance();

    /**
     * @brief GetConfirmationFromUser raises a style-independent message box asking the user if they wish to keep their changes.
     * @return  true if the user presses "Apply", false if "Cancel" or if 10 seconds pass
     */
    bool GetConfirmationFromUser();

    /**
     * @brief SetApplicationPreferences saves the current selections to the global ApplicationPreferences (allowing them to be saved to disk later).
     */
    void UpdateAppearancePreferences();
};

#endif // APPEARANCEPAGE_H
