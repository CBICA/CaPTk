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
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_AppearancePage
{
public:
    QGridLayout *gridLayout;
    QLabel *currentFontLabel;
    QLabel *label;
    QLabel *label_2;
    QComboBox *themeComboBox;
    QPushButton *selectFontBtn;
    QSpacerItem *verticalSpacer;

    void setupUi(QWidget *AppearancePage)
    {
        if (AppearancePage->objectName().isEmpty())
            AppearancePage->setObjectName(QStringLiteral("AppearancePage"));
        AppearancePage->resize(400, 300);
        gridLayout = new QGridLayout(AppearancePage);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        currentFontLabel = new QLabel(AppearancePage);
        currentFontLabel->setObjectName(QStringLiteral("currentFontLabel"));

        gridLayout->addWidget(currentFontLabel, 0, 2, 1, 1);

        label = new QLabel(AppearancePage);
        label->setObjectName(QStringLiteral("label"));
        label->setMaximumSize(QSize(100, 16777215));

        gridLayout->addWidget(label, 0, 1, 1, 1);

        label_2 = new QLabel(AppearancePage);
        label_2->setObjectName(QStringLiteral("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        themeComboBox = new QComboBox(AppearancePage);
        themeComboBox->addItem(QString());
        themeComboBox->addItem(QString());
        themeComboBox->setObjectName(QStringLiteral("themeComboBox"));

        gridLayout->addWidget(themeComboBox, 1, 1, 1, 1);

        selectFontBtn = new QPushButton(AppearancePage);
        selectFontBtn->setObjectName(QStringLiteral("selectFontBtn"));
        selectFontBtn->setMaximumSize(QSize(75, 16777215));

        gridLayout->addWidget(selectFontBtn, 0, 0, 1, 1);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);
        gridLayout->addItem(verticalSpacer, 2, 1, 1, 1);

        retranslateUi(AppearancePage);

        QMetaObject::connectSlotsByName(AppearancePage);
    } // setupUi

    void retranslateUi(QWidget *AppearancePage)
    {
        AppearancePage->setWindowTitle(QApplication::translate("AppearancePage", "Form", nullptr));
        currentFontLabel->setText(QApplication::translate("AppearancePage", "text", nullptr));
        label->setText(QApplication::translate("AppearancePage", "Current Font:", nullptr));
        label_2->setText(QApplication::translate("AppearancePage", "Theme", nullptr));
        themeComboBox->setItemText(0, QApplication::translate("AppearancePage", "Light", nullptr));
        themeComboBox->setItemText(1, QApplication::translate("AppearancePage", "Dark", nullptr));

        selectFontBtn->setText(QApplication::translate("AppearancePage", "Select Font", nullptr));
    } // retranslateUi

};

namespace Ui {
    class AppearancePage: public Ui_AppearancePage {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_APPEARANCEPAGE_H
