///////////////////////////////////////////////////////////////////////////////////////
// fPopulationAtlasDialog.h
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/sbia/software/license.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef ui_fLungFieldDialog_H
#define ui_fLungFieldDialog_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_fLungFieldDialog
{
public:
	QVBoxLayout *verticalLayout;
	QHBoxLayout *horizontalLayout_2;
	QLabel *maskFileLabel;
	QLineEdit *maskFileName;
	QPushButton *maskFileBrowseButton;
	QHBoxLayout *horizontalLayout_4;
	QLabel *outputDirLabel;
	QLineEdit *outputDirName;
	QPushButton *outputDirBrowseButton;
	QHBoxLayout *horizontalLayout_3;
	QSpacerItem *horizontalSpacer;
	QPushButton *okButton;
	QPushButton *cancelButton;
	QSpacerItem *horizontalSpacer_2;

	void setupUi(QDialog *fLungFieldDialog)
	{
		if (fLungFieldDialog->objectName().isEmpty())
			fLungFieldDialog->setObjectName(QStringLiteral("fLungFieldDialog"));
		fLungFieldDialog->resize(397, 151);
		verticalLayout = new QVBoxLayout(fLungFieldDialog);
		verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
		horizontalLayout_2 = new QHBoxLayout();
		horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
		maskFileLabel = new QLabel(fLungFieldDialog);
		maskFileLabel->setObjectName(QStringLiteral("maskFileLabel"));
		maskFileLabel->setMinimumSize(QSize(37, 0));

		horizontalLayout_2->addWidget(maskFileLabel);

		maskFileName = new QLineEdit(fLungFieldDialog);
		maskFileName->setObjectName(QStringLiteral("maskFileName"));
		maskFileName->setMaximumSize(QSize(173, 16777215));

		horizontalLayout_2->addWidget(maskFileName);

		maskFileBrowseButton = new QPushButton(fLungFieldDialog);
		maskFileBrowseButton->setObjectName(QStringLiteral("maskFileBrowseButton"));

		horizontalLayout_2->addWidget(maskFileBrowseButton);


		verticalLayout->addLayout(horizontalLayout_2);

		horizontalLayout_4 = new QHBoxLayout();
		horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
		outputDirLabel = new QLabel(fLungFieldDialog);
		outputDirLabel->setObjectName(QStringLiteral("outputDirLabel"));

		horizontalLayout_4->addWidget(outputDirLabel);

		outputDirName = new QLineEdit(fLungFieldDialog);
		outputDirName->setObjectName(QStringLiteral("outputDirName"));

		horizontalLayout_4->addWidget(outputDirName);

		outputDirBrowseButton = new QPushButton(fLungFieldDialog);
		outputDirBrowseButton->setObjectName(QStringLiteral("outputDirBrowseButton"));

		horizontalLayout_4->addWidget(outputDirBrowseButton);


		verticalLayout->addLayout(horizontalLayout_4);

		horizontalLayout_3 = new QHBoxLayout();
		horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
		horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

		horizontalLayout_3->addItem(horizontalSpacer);

		okButton = new QPushButton(fLungFieldDialog);
		okButton->setObjectName(QStringLiteral("okButton"));

		horizontalLayout_3->addWidget(okButton);

		cancelButton = new QPushButton(fLungFieldDialog);
		cancelButton->setObjectName(QStringLiteral("cancelButton"));

		horizontalLayout_3->addWidget(cancelButton);

		horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

		horizontalLayout_3->addItem(horizontalSpacer_2);


		verticalLayout->addLayout(horizontalLayout_3);


		retranslateUi(fLungFieldDialog);

		QMetaObject::connectSlotsByName(fLungFieldDialog);
	} // setupUi

	void retranslateUi(QDialog *fLungFieldDialog)
	{
		fLungFieldDialog->setWindowTitle(QApplication::translate("fLungFieldDialog", "Lung Field", nullptr));
		maskFileLabel->setText(QApplication::translate("fLungFieldDialog", "Mask File", nullptr));
		maskFileBrowseButton->setText(QApplication::translate("fLungFieldDialog", "Browse", nullptr));
		outputDirLabel->setText(QApplication::translate("fLungFieldDialog", "Output Directory", nullptr));
		outputDirBrowseButton->setText(QApplication::translate("fLungFieldDialog", "Browse", nullptr));
		okButton->setText(QApplication::translate("fLungFieldDialog", "Ok", nullptr));
		cancelButton->setText(QApplication::translate("fLungFieldDialog", "Cancel", nullptr));
	} // retranslateUi

};

namespace Ui {
	class fLungFieldDialog : public Ui_fLungFieldDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif 





