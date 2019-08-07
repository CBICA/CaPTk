/**
\file  ui_AppearancePage.h

\brief Declaration of AppearancePage class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef UI_APPEARANCEPAGE_H
#define UI_APPEARANCEPAGE_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_AppearancePage
{
public:
    QGridLayout *gridLayout;
    QLabel *label;
    QLabel *label_2;
    QComboBox *comboBox;
    QPushButton *pushButton;
    QLabel *label_3;

    void setupUi(QWidget *AppearancePage)
    {
        if (AppearancePage->objectName().isEmpty())
            AppearancePage->setObjectName(QStringLiteral("AppearancePage"));
        AppearancePage->resize(400, 300);
        gridLayout = new QGridLayout(AppearancePage);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        label = new QLabel(AppearancePage);
        label->setObjectName(QStringLiteral("label"));
        label->setMaximumSize(QSize(100, 16777215));

        gridLayout->addWidget(label, 0, 1, 1, 1);

        label_2 = new QLabel(AppearancePage);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        comboBox = new QComboBox(AppearancePage);
        comboBox->setObjectName(QStringLiteral("comboBox"));

        gridLayout->addWidget(comboBox, 1, 1, 1, 1);

        pushButton = new QPushButton(AppearancePage);
        pushButton->setObjectName(QStringLiteral("pushButton"));
        pushButton->setMaximumSize(QSize(50, 16777215));

        gridLayout->addWidget(pushButton, 0, 0, 1, 1);

        label_3 = new QLabel(AppearancePage);
        label_3->setObjectName(QStringLiteral("label_3"));

        gridLayout->addWidget(label_3, 0, 2, 1, 1);


        retranslateUi(AppearancePage);

        QMetaObject::connectSlotsByName(AppearancePage);
    } // setupUi

    void retranslateUi(QWidget *AppearancePage)
    {
        AppearancePage->setWindowTitle(QApplication::translate("AppearancePage", "Form", nullptr));
        label->setText(QApplication::translate("AppearancePage", "Current Font:", nullptr));
        label_2->setText(QApplication::translate("AppearancePage", "Theme", nullptr));
        pushButton->setText(QApplication::translate("AppearancePage", "Font", nullptr));
        label_3->setText(QApplication::translate("AppearancePage", "text", nullptr));
    } // retranslateUi

};

namespace Ui {
    class AppearancePage: public Ui_AppearancePage {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_APPEARANCEPAGE_H
