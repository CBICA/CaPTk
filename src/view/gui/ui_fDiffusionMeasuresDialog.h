#ifndef ui_fDiffusionEstimator_H
#define ui_fDiffusionEstimator_H

#include <QtCore/QVariant>
// #include <QtGui/QAction>
// #include <QtGui/QApplication>
// #include <QtGui/QButtonGroup>
// #include <QtGui/QCheckBox>
// #include <QtGui/QDialog>
// #include <QtGui/QFrame>
// #include <QtGui/QGridLayout>
// #include <QtGui/QHBoxLayout>
// #include <QtGui/QHeaderView>
// #include <QtGui/QLabel>
// #include <QtGui/QLineEdit>
// #include <QtGui/QPushButton>
// #include <QtGui/QSpacerItem>
// #include <QtGui/QSpinBox>
// NEW CHANGES
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QSpinBox>
#include <QGroupBox>
#include <QtWidgets/QRadioButton>

QT_BEGIN_NAMESPACE

class ui_fDiffusionEstimator
{
public:
	QGridLayout *gridLayout;
	QVBoxLayout * verticalLayout;


	QGroupBox *inputGroupBox;
	QGridLayout *inputGridLayout;

	QLabel		*inputImageLabel;
	QLineEdit	*inputImageName;
	QPushButton *inputImageButton;

	QLabel		*inputMaskLabel;
	QLineEdit	*inputMaskName;
	QPushButton *inputMaskButton;

	QLabel		*inputBvalLabel;
	QLineEdit	*inputBvalName;
	QPushButton *inputBvalButton;

	QLabel		*inputBvecLabel;
	QLineEdit	*inputBvecName;
	QPushButton *inputBvecButton;

	QGroupBox	*outputGroupBox;
	QGridLayout *outputGridLayout;
	QLineEdit	*outputImageName;
	QLabel		*outputImageLabel;
	QPushButton *outputImageButton;

	QPushButton *confirmButton;
	QPushButton *cancelButton;

	QCheckBox* m_ax;
	QCheckBox* m_fa;
	QCheckBox* m_rad;
	QCheckBox* m_tr;
	QCheckBox* m_bzero;

	QHBoxLayout * horizontalLayout;

  QLabel *longRunningWarning;

	void setupUi(QDialog *fDiffusionEstimator)
	{
		if (fDiffusionEstimator->objectName().isEmpty())
			fDiffusionEstimator->setObjectName(QString::fromUtf8("fDiffusionEstimator"));
    //fDiffusionEstimator->setWindowModality(Qt::NonModal);
    //fDiffusionEstimator->setModal(true);
		fDiffusionEstimator->resize(400, 200); // needs to be screenSize dependent
		QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
		sizePolicy.setHorizontalStretch(0);
		sizePolicy.setVerticalStretch(0);
		sizePolicy.setHeightForWidth(fDiffusionEstimator->sizePolicy().hasHeightForWidth());
		fDiffusionEstimator->setSizePolicy(sizePolicy);
		fDiffusionEstimator->setMinimumSize(QSize(0, 0));

		fDiffusionEstimator->setModal(true);
		gridLayout = new QGridLayout(fDiffusionEstimator);
		gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

		//input
		inputGroupBox = new QGroupBox(fDiffusionEstimator);
		inputGroupBox->setTitle(QString::fromStdString("Input Data"));

		inputGridLayout = new QGridLayout(inputGroupBox);
		inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

		inputImageLabel = new QLabel(inputGroupBox);
		sizePolicy.setHeightForWidth(inputImageLabel->sizePolicy().hasHeightForWidth());
		inputImageLabel->setSizePolicy(sizePolicy);

		inputImageName = new QLineEdit("");
		inputImageName->setObjectName(QString::fromUtf8("inputImageName"));
		sizePolicy.setHeightForWidth(inputImageName->sizePolicy().hasHeightForWidth());
		inputImageName->setSizePolicy(sizePolicy);
		inputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

		inputImageButton = new QPushButton(inputGroupBox);
		inputImageButton->setObjectName(QString::fromUtf8("inputImageButton"));
		inputImageButton->setText(QString("Browse"));


		inputMaskLabel = new QLabel(inputGroupBox);
		sizePolicy.setHeightForWidth(inputMaskLabel->sizePolicy().hasHeightForWidth());
		inputMaskLabel->setSizePolicy(sizePolicy);

		inputMaskName = new QLineEdit("");
		inputMaskName->setObjectName(QString::fromUtf8("inputMaskName"));
		sizePolicy.setHeightForWidth(inputMaskName->sizePolicy().hasHeightForWidth());
		inputMaskName->setSizePolicy(sizePolicy);
		inputMaskName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

		inputMaskButton = new QPushButton(inputGroupBox);
		inputMaskButton->setObjectName(QString::fromUtf8("inputMaskButton"));
		inputMaskButton->setText(QString("Browse"));




		inputBvalLabel = new QLabel(inputGroupBox);
		sizePolicy.setHeightForWidth(inputBvalLabel->sizePolicy().hasHeightForWidth());
		inputBvalLabel->setSizePolicy(sizePolicy);

		inputBvalName = new QLineEdit("");
		inputBvalName->setObjectName(QString::fromUtf8("inputBvalName"));
		sizePolicy.setHeightForWidth(inputBvalName->sizePolicy().hasHeightForWidth());
		inputBvalName->setSizePolicy(sizePolicy);
		inputBvalName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

		inputBvalButton = new QPushButton(inputGroupBox);
		inputBvalButton->setObjectName(QString::fromUtf8("inputBvalButton"));
		inputBvalButton->setText(QString("Browse"));




		inputBvecLabel = new QLabel(inputGroupBox);
		sizePolicy.setHeightForWidth(inputBvecLabel->sizePolicy().hasHeightForWidth());
		inputBvecLabel->setSizePolicy(sizePolicy);

		inputBvecName = new QLineEdit("");
		inputBvecName->setObjectName(QString::fromUtf8("inputBvecName"));
		sizePolicy.setHeightForWidth(inputBvecName->sizePolicy().hasHeightForWidth());
		inputBvecName->setSizePolicy(sizePolicy);
		inputBvecName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

		inputBvecButton = new QPushButton(inputGroupBox);
		inputBvecButton->setObjectName(QString::fromUtf8("inputBvecButton"));
		inputBvecButton->setText(QString("Browse"));





		m_ax = new QCheckBox("Axial Diffusivity");
		m_ax->setEnabled(true);
		m_fa = new QCheckBox("Fractional Anisotropy");
		m_fa->setEnabled(true);
		m_rad = new QCheckBox("Radial Diffusivity");
		m_rad->setEnabled(true);
		m_tr  = new QCheckBox("Apparent Diffusion Coefficient");
		m_tr->setEnabled(true);
		m_bzero = new QCheckBox("Extract averaged B0 image");
		m_bzero->setEnabled(true);


		inputGridLayout->addWidget(inputImageLabel, 0, 0, 1, 1);
		inputGridLayout->addWidget(inputImageName, 0, 1, 1, 1);
		inputGridLayout->addWidget(inputImageButton, 0, 2, 1, 1);

		inputGridLayout->addWidget(inputMaskLabel, 1, 0, 1, 1);
		inputGridLayout->addWidget(inputMaskName, 1, 1, 1, 1);
		inputGridLayout->addWidget(inputMaskButton, 1, 2, 1, 1);

		inputGridLayout->addWidget(inputBvalLabel, 2, 0, 1, 1);
		inputGridLayout->addWidget(inputBvalName, 2, 1, 1, 1);
		inputGridLayout->addWidget(inputBvalButton, 2, 2, 1, 1);

		inputGridLayout->addWidget(inputBvecLabel, 3, 0, 1, 1);
		inputGridLayout->addWidget(inputBvecName, 3, 1, 1, 1);
		inputGridLayout->addWidget(inputBvecButton, 3, 2, 1, 1);


		inputGridLayout->addWidget(m_fa, 4, 0, 1, 1);
		inputGridLayout->addWidget(m_rad, 5, 0, 1, 1);
		inputGridLayout->addWidget(m_ax, 6, 0, 1, 1);
		inputGridLayout->addWidget(m_tr, 7, 0, 1, 1);
		inputGridLayout->addWidget(m_bzero, 8, 0, 1, 1);

		// output
		outputGroupBox = new QGroupBox(fDiffusionEstimator);
		outputGroupBox->setTitle(QString::fromStdString("Output"));

		outputGridLayout = new QGridLayout(outputGroupBox);
		outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

		//outputImageLabel = new QLabel(outputGroupBox);
		//sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
		//outputImageLabel->setSizePolicy(sizePolicy);

		outputImageName = new QLineEdit("");
		outputImageName->setObjectName(QString::fromUtf8("outputImageName"));
		sizePolicy.setHeightForWidth(outputImageName->sizePolicy().hasHeightForWidth());
		outputImageName->setSizePolicy(sizePolicy);
		outputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

		outputImageButton = new QPushButton(outputGroupBox);
		outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
		outputImageButton->setText(QString("Browse"));

    longRunningWarning = new QLabel(outputGroupBox);
    sizePolicy.setHeightForWidth(longRunningWarning->sizePolicy().hasHeightForWidth());
    longRunningWarning->setSizePolicy(sizePolicy);
    longRunningWarning->setAlignment(Qt::AlignRight);
    longRunningWarning->setText("NOTE: CaPTk will not let you interact with the UI while this application runs.");

		//outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
		outputGridLayout->addWidget(outputImageName, 0, 0, 1, 1);
		outputGridLayout->addWidget(outputImageButton, 0, 1, 1, 1);
    outputGridLayout->addWidget(longRunningWarning, 1, 0, 1, 2);

		// put the layout in perspective
		gridLayout->addWidget(inputGroupBox, 1, 0, 1, 2);
		gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


		confirmButton = new QPushButton(fDiffusionEstimator);
		confirmButton->setObjectName(QString::fromUtf8("confirm"));
		//confirmButton->setIcon(ButtonIcon);
		//confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

		cancelButton = new QPushButton(fDiffusionEstimator);
    cancelButton->setObjectName(QString::fromUtf8("Cancel"));
		//cancelButton->setIcon(ButtonIcon);
		//cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

		gridLayout->addWidget(confirmButton, 3, 0, 1, 1);
		gridLayout->addWidget(cancelButton, 3, 1, 1, 1);

		retranslateUi(fDiffusionEstimator);

		QMetaObject::connectSlotsByName(fDiffusionEstimator);
	} // setupUi

	void retranslateUi(QDialog *fDiffusionEstimator)
	{
		// fDiffusionEstimator->setWindowTitle(QApplication::translate("fDiffusionEstimator", "Diffusion Derivatives", 0, QApplication::UnicodeUTF8));
		// confirmButton->setText(QApplication::translate("fDiffusionEstimator", "Confirm", 0, QApplication::UnicodeUTF8));
		// cancelButton->setText(QApplication::translate("fDiffusionEstimator", "Cancel", 0, QApplication::UnicodeUTF8));
		// m_fa->setText(QApplication::translate("fDiffusionEstimator", "Fractional Anisotropy", 0, QApplication::UnicodeUTF8));
		// m_ax->setText(QApplication::translate("fDiffusionEstimator", "Axial Diffusivity", 0, QApplication::UnicodeUTF8));
		// m_rad->setText(QApplication::translate("fDiffusionEstimator", "Radial Diffusivity", 0, QApplication::UnicodeUTF8));
		// m_tr->setText(QApplication::translate("fDiffusionEstimator", "Apparent Diffusion Coefficient", 0, QApplication::UnicodeUTF8));
      //
		// inputImageLabel->setText(QApplication::translate("fDiffusionEstimator", "DWI Image:", 0, QApplication::UnicodeUTF8));
		// inputMaskLabel->setText(QApplication::translate("fDiffusionEstimator", "Mask:", 0, QApplication::UnicodeUTF8));
		// inputBvalLabel->setText(QApplication::translate("fDiffusionEstimator", "Bval File:", 0, QApplication::UnicodeUTF8));
		// inputBvecLabel->setText(QApplication::translate("fDiffusionEstimator", "Bvec File:", 0, QApplication::UnicodeUTF8));
		// NEW CHANGES
		fDiffusionEstimator->setWindowTitle(QApplication::translate("fDiffusionEstimator", "Diffusion Derivatives", 0));
		confirmButton->setText(QApplication::translate("fDiffusionEstimator", "Confirm", 0));
		cancelButton->setText(QApplication::translate("fDiffusionEstimator", "Cancel", 0));
		m_fa->setText(QApplication::translate("fDiffusionEstimator", "Fractional Anisotropy", 0));
		m_ax->setText(QApplication::translate("fDiffusionEstimator", "Axial Diffusivity", 0));
		m_rad->setText(QApplication::translate("fDiffusionEstimator", "Radial Diffusivity", 0));
		m_tr->setText(QApplication::translate("fDiffusionEstimator", "Apparent Diffusion Coefficient", 0));
		m_bzero->setText(QApplication::translate("fDiffusionEstimator", "Extract averaged b0 image", 0));

		inputImageLabel->setText(QApplication::translate("fDiffusionEstimator", "DWI Image:", 0));
		inputMaskLabel->setText(QApplication::translate("fDiffusionEstimator", "Mask:", 0));
		inputBvalLabel->setText(QApplication::translate("fDiffusionEstimator", "Bval File:", 0));
		inputBvecLabel->setText(QApplication::translate("fDiffusionEstimator", "Bvec File:", 0));
	} // retranslateUi
};

namespace Ui {
	class fDiffusionEstimator : public ui_fDiffusionEstimator {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fDiffusionEstimator_H
