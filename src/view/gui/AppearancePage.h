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
	void OnOkay() override;
	void OnCancel() override;

private:
    Ui::AppearancePage *ui;
    QFont m_selectedFont;
	QString m_currentStyleSheet;
};

#endif // APPEARANCEPAGE_H
