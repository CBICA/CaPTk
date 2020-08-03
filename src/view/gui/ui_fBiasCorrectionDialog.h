#ifndef ui_fBiasCorrectionDialog_H
#define ui_fBiasCorrectionDialog_H

#include <QtCore/QVariant>

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

class ui_fBiasCorrectionDialog
{
public:
    QGridLayout* gridLayout;
    QVBoxLayout* verticalLayout;

    QGroupBox* optionsGroupBox;
    QGridLayout* optionsGridLayout;
    QLabel* options_splineOrderLabel;
    QLineEdit* options_splineOrderName;
    QLabel* options_otsuBinsLabel;
    QLineEdit* options_otsuBinsName;
    QLabel* options_maxIterationsLabel;
    QLineEdit* options_maxIterationsName;
    QLabel* options_fittingLevelsLabel;
    QLineEdit* options_fittingLevelsName;
    QLabel* options_filterNoiseLabel;
    QLineEdit* options_filterNoiseName;
    QLabel* options_fwhmLabel;
    QLineEdit* options_fwhmName;

    // radio buttons
    QGroupBox* modeGroupBox;
    QGridLayout* modeGridLayout;
    QButtonGroup* mode_radioButtons;
    QRadioButton* mode_N3Button;
    QLabel* mode_N3Label;
    QRadioButton* mode_N4Button;
    QLabel* mode_N4Label;


    QGroupBox* outputGroupBox;
    QGridLayout* outputGridLayout;
    QLineEdit* outputImageName;
    QLabel* outputImageLabel;
    QPushButton* outputImageButton;

    QPushButton* confirmButton;
    QPushButton* cancelButton;

    QFrame* line_3;

    QHBoxLayout* horizontalLayout;

    QCheckBox* wholeImageThresholdCheckBoxBox;

    void setupUi(QDialog* fBiasCorrectorDialog)
    {

        if (fBiasCorrectorDialog->objectName().isEmpty())
            fBiasCorrectorDialog->setObjectName(QString::fromUtf8("fBiasCorrectorDialog"));
        //fBiasCorrectorDialog->setWindowModality(Qt::NonModal);
        fBiasCorrectorDialog->resize(256, 256); // needs to be screenSize dependent
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(fBiasCorrectorDialog->sizePolicy().hasHeightForWidth());
        fBiasCorrectorDialog->setSizePolicy(sizePolicy);
        fBiasCorrectorDialog->setMinimumSize(QSize(0, 0));

        //fBiasCorrectorDialog->setModal(true);
        gridLayout = new QGridLayout(fBiasCorrectorDialog);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));

        modeGroupBox = new QGroupBox(fBiasCorrectorDialog);
        mode_radioButtons = new QButtonGroup(modeGroupBox);
        modeGroupBox->setTitle(QString::fromStdString("modeGridLayout"));
        modeGridLayout = new QGridLayout(modeGroupBox);
        mode_N3Label = new QLabel();
        mode_N3Button = new QRadioButton();
        mode_N4Label = new QLabel();
        mode_N4Button = new QRadioButton();
        int gridRowCounter = 0;

        modeGridLayout->addWidget(mode_N3Label, gridRowCounter, 0);
        modeGridLayout->addWidget(mode_N3Button, gridRowCounter, 1);

        mode_N3Label->setText("N3");
        sizePolicy.setHeightForWidth(mode_N3Label->sizePolicy().hasHeightForWidth());
        mode_N3Label->setSizePolicy(sizePolicy);
        mode_N3Label->setObjectName(QString::fromUtf8("mode_N3Label"));
        mode_N3Label->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


        gridRowCounter++;
        modeGridLayout->addWidget(mode_N4Label, gridRowCounter, 0);
        modeGridLayout->addWidget(mode_N4Button, gridRowCounter, 1);
        
        mode_N4Label->setText("N4");
        sizePolicy.setHeightForWidth(mode_N4Label->sizePolicy().hasHeightForWidth());
        mode_N4Label->setSizePolicy(sizePolicy);
        mode_N4Label->setObjectName(QString::fromUtf8("mode_N4Label"));
        mode_N4Label->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);


        mode_radioButtons->addButton(mode_N3Button, 0);
        mode_radioButtons->addButton(mode_N4Button, 1);
        mode_radioButtons->setExclusive(true);


        // options
        optionsGroupBox = new QGroupBox(fBiasCorrectorDialog);
        optionsGroupBox->setTitle(QString::fromStdString("Options"));

        optionsGridLayout = new QGridLayout(optionsGroupBox);
        optionsGridLayout->setObjectName(QString::fromUtf8("optionsGridLayout"));
        gridRowCounter = 0;

        options_splineOrderLabel = new QLabel(optionsGroupBox);
        sizePolicy.setHeightForWidth(options_splineOrderLabel->sizePolicy().hasHeightForWidth());
        options_splineOrderLabel->setSizePolicy(sizePolicy);
        options_splineOrderLabel->setText("Spline order");
        options_splineOrderName = new QLineEdit("");
        options_splineOrderName->setObjectName(QString::fromUtf8("options_splineOrderName"));
        sizePolicy.setHeightForWidth(options_splineOrderName->sizePolicy().hasHeightForWidth());
        options_splineOrderName->setSizePolicy(sizePolicy);
        options_splineOrderName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);
        options_splineOrderName->setText("");

        optionsGridLayout->addWidget(options_splineOrderLabel, gridRowCounter, 0);
        optionsGridLayout->addWidget(options_splineOrderName, gridRowCounter, 1);

        options_otsuBinsLabel = new QLabel(optionsGroupBox);
        sizePolicy.setHeightForWidth(options_otsuBinsLabel->sizePolicy().hasHeightForWidth());
        options_otsuBinsLabel->setSizePolicy(sizePolicy);
        options_otsuBinsLabel->setText("Otsu bins");
        options_otsuBinsName = new QLineEdit("");
        options_otsuBinsName->setObjectName(QString::fromUtf8("options_otsuBinsName"));
        sizePolicy.setHeightForWidth(options_otsuBinsName->sizePolicy().hasHeightForWidth());
        options_otsuBinsName->setSizePolicy(sizePolicy);
        options_otsuBinsName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);
        options_otsuBinsName->setText("");

        gridRowCounter++;
        optionsGridLayout->addWidget(options_otsuBinsLabel, gridRowCounter, 0);
        optionsGridLayout->addWidget(options_otsuBinsName, gridRowCounter, 1);

        options_maxIterationsLabel = new QLabel(optionsGroupBox);
        sizePolicy.setHeightForWidth(options_maxIterationsLabel->sizePolicy().hasHeightForWidth());
        options_maxIterationsLabel->setSizePolicy(sizePolicy);
        options_maxIterationsLabel->setText("Max Iterations");
        options_maxIterationsName = new QLineEdit("");
        options_maxIterationsName->setObjectName(QString::fromUtf8("options_maxIterationsName"));
        sizePolicy.setHeightForWidth(options_maxIterationsName->sizePolicy().hasHeightForWidth());
        options_maxIterationsName->setSizePolicy(sizePolicy);
        options_maxIterationsName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);
        options_maxIterationsName->setText("");

        gridRowCounter++;
        optionsGridLayout->addWidget(options_maxIterationsLabel, gridRowCounter, 0);
        optionsGridLayout->addWidget(options_maxIterationsName, gridRowCounter, 1);

        options_fittingLevelsLabel = new QLabel(optionsGroupBox);
        sizePolicy.setHeightForWidth(options_fittingLevelsLabel->sizePolicy().hasHeightForWidth());
        options_fittingLevelsLabel->setSizePolicy(sizePolicy);
        options_fittingLevelsLabel->setText("Fitting levels");
        options_fittingLevelsName = new QLineEdit("");
        options_fittingLevelsName->setObjectName(QString::fromUtf8("options_fittingLevelsName"));
        sizePolicy.setHeightForWidth(options_fittingLevelsName->sizePolicy().hasHeightForWidth());
        options_fittingLevelsName->setSizePolicy(sizePolicy);
        options_fittingLevelsName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);
        options_fittingLevelsName->setText("");

        gridRowCounter++;
        optionsGridLayout->addWidget(options_fittingLevelsLabel, gridRowCounter, 0);
        optionsGridLayout->addWidget(options_fittingLevelsName, gridRowCounter, 1);

        options_filterNoiseLabel = new QLabel(optionsGroupBox);
        sizePolicy.setHeightForWidth(options_filterNoiseLabel->sizePolicy().hasHeightForWidth());
        options_filterNoiseLabel->setSizePolicy(sizePolicy);
        options_filterNoiseLabel->setText("Filter noise");
        options_filterNoiseName = new QLineEdit("");
        options_filterNoiseName->setObjectName(QString::fromUtf8("options_filterNoiseName"));
        sizePolicy.setHeightForWidth(options_filterNoiseName->sizePolicy().hasHeightForWidth());
        options_filterNoiseName->setSizePolicy(sizePolicy);
        options_filterNoiseName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);
        options_filterNoiseName->setText("");
        options_filterNoiseName->setValidator(new QIntValidator(0, 100, optionsGroupBox));

        gridRowCounter++;
        optionsGridLayout->addWidget(options_filterNoiseLabel, gridRowCounter, 0);
        optionsGridLayout->addWidget(options_filterNoiseName, gridRowCounter, 1);

        options_fwhmLabel = new QLabel(optionsGroupBox);
        sizePolicy.setHeightForWidth(options_fwhmLabel->sizePolicy().hasHeightForWidth());
        options_fwhmLabel->setSizePolicy(sizePolicy);
        options_fwhmLabel->setText("FWHM");
        options_fwhmName = new QLineEdit("");
        options_fwhmName->setObjectName(QString::fromUtf8("options_fwhmName"));
        sizePolicy.setHeightForWidth(options_fwhmName->sizePolicy().hasHeightForWidth());
        options_fwhmName->setSizePolicy(sizePolicy);
        options_fwhmName->setAlignment(Qt::AlignRight | Qt::AlignTrailing | Qt::AlignVCenter);
        options_fwhmName->setText("");

        gridRowCounter++;
        optionsGridLayout->addWidget(options_fwhmLabel, gridRowCounter, 0);
        optionsGridLayout->addWidget(options_fwhmName, gridRowCounter, 1);


        // output
        outputGroupBox = new QGroupBox(fBiasCorrectorDialog);
        outputGroupBox->setTitle(QString::fromStdString("Output"));

        outputGridLayout = new QGridLayout(outputGroupBox);
        outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

        outputImageLabel = new QLabel(outputGroupBox);
        sizePolicy.setHeightForWidth(outputImageLabel->sizePolicy().hasHeightForWidth());
        outputImageLabel->setSizePolicy(sizePolicy);

        outputImageName = new QLineEdit("");
        outputImageName->setObjectName(QString::fromUtf8("outputImageName"));
        sizePolicy.setHeightForWidth(outputImageName->sizePolicy().hasHeightForWidth());
        outputImageName->setSizePolicy(sizePolicy);
        outputImageName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        outputImageButton = new QPushButton(outputGroupBox);
        outputImageButton->setObjectName(QString::fromUtf8("outputImageButton"));
        //outputImageButton->setIcon(ButtonIcon);
        //outputImageButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent
        outputImageButton->setText(QString("Browse"));

        outputGridLayout->addWidget(outputImageLabel, 0, 0, 1, 1);
        outputGridLayout->addWidget(outputImageName, 1, 0, 1, 1);
        outputGridLayout->addWidget(outputImageButton, 1, 1, 1, 1);

        // put the layout in perspective
        gridLayout->addWidget(modeGroupBox, 0, 0, 1, 2);
        gridLayout->addWidget(optionsGroupBox, 1, 0, 1, 2);
        gridLayout->addWidget(outputGroupBox, 2, 0, 1, 2);


        confirmButton = new QPushButton(fBiasCorrectorDialog);
        confirmButton->setObjectName(QString::fromUtf8("confirm"));
        //confirmButton->setIcon(ButtonIcon);
        //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

        cancelButton = new QPushButton(fBiasCorrectorDialog);
        cancelButton->setObjectName(QString::fromUtf8("Cancel"));
        //cancelButton->setIcon(ButtonIcon);
        //cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent

        gridLayout->addWidget(confirmButton, 4, 0, 1, 1);
        gridLayout->addWidget(cancelButton, 4, 1, 1, 1);

        retranslateUi(fBiasCorrectorDialog);

        QMetaObject::connectSlotsByName(fBiasCorrectorDialog);
    } // setupUi

    void retranslateUi(QDialog* fBiasCorrectorDialog)
    {
        fBiasCorrectorDialog->setWindowTitle(QApplication::translate("fBiasCorrectorDialog", "Bias Correction", 0));
        confirmButton->setText(QApplication::translate("fBiasCorrectorDialog", "Confirm", 0));
        cancelButton->setText(QApplication::translate("fBiasCorrectorDialog", "Cancel", 0));

    } // retranslateUi

};

namespace Ui {
    class fBiasCorrectionDialog : public ui_fBiasCorrectionDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fBiasCorrectionDialog_H
