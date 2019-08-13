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

namespace Ui {
class AppearancePage;
}

class AppearancePage : public QWidget
{
    Q_OBJECT

public:
    explicit AppearancePage(QWidget *parent = nullptr);
    ~AppearancePage();

public slots:
	void OnSelectFontButtonClicked();

private:
    Ui::AppearancePage *ui;
    QFont m_selectFont;
};

#endif // APPEARANCEPAGE_H
