#ifndef ui_fRegistrationDialog_H
#define ui_fRegistrationDialog_H

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
#include <QComboBox>

QT_BEGIN_NAMESPACE

class ui_fRegistrationDialog
{
public:
    QGridLayout * gridLayout;
    QVBoxLayout * verticalLayout;
    QHBoxLayout *buttonGridLayout;

    QGroupBox *modeGroupBox;
    QGroupBox *registrationGroupBox;
    QGroupBox *inputGroupBox;
    QGroupBox *matrixGroupBox;
    QGroupBox *outputGroupBox;
    QGroupBox *buttonGroupBox;

    QGridLayout *outputGridLayout;
    QGridLayout *inputGridLayout;
    QGridLayout *matrixGridLayout;
    QGridLayout *modeGridLayout;
    QGridLayout *registrationGridLayout;

    QLineEdit *fixedFileName;
    QLineEdit *iterations;
    QLineEdit *nccRadii;
    QLineEdit *movingFileName1;
    QLineEdit *movingFileName2;
    QLineEdit *movingFileName3;
    QLineEdit *movingFileName4;
    QLineEdit *movingFileName5;

    QLabel *fixedFileLabel;
    QLabel *movingFileLabel1;
    QLabel *movingFileLabel2;
    QLabel *movingFileLabel3;
    QLabel *movingFileLabel4;
    QLabel *movingFileLabel5;
    QLabel *options_NCC_radii;
    QLabel *options_iterations;

    QPushButton *fixedFileButton;
    QPushButton *movingFileButton1;
    QPushButton *movingFileButton2;
    QPushButton *movingFileButton3;
    QPushButton *movingFileButton4;
    QPushButton *movingFileButton5;
    QPushButton *addMoreButton;

    QLabel *matrixRadioLabel;
    QCheckBox *matrixRadioButton;

    QLineEdit *matrixFileName;
    QLabel *matrixFileLabel;
    QPushButton *matrixFileButton;

    QLineEdit *matrixName1;
    QLabel *matrixLabel1;
    QPushButton *matrixButton1;

    QLineEdit *matrixName2;
    QLabel *matrixLabel2;
    QPushButton *matrixButton2;

    QLineEdit *matrixName3;
    QLabel *matrixLabel3;
    QPushButton *matrixButton3;

    QLineEdit *matrixName4;
    QLabel *matrixLabel4;
    QPushButton *matrixButton4;

    QLineEdit *matrixName5;
    QLabel *matrixLabel5;
    QPushButton *matrixButton5;
    
    QLineEdit *movingFileOutputName1;
    QLabel *movingFileOutputLabel1;
    QPushButton *movingFileOutputButton1;

    QLineEdit *movingFileOutputName2;
    QLabel *movingFileOutputLabel2;
    QPushButton *movingFileOutputButton2;

    QLineEdit *movingFileOutputName3;
    QLabel *movingFileOutputLabel3;
    QPushButton *movingFileOutputButton3;

    QLineEdit *movingFileOutputName4;
    QLabel *movingFileOutputLabel4;
    QPushButton *movingFileOutputButton4;

    QLineEdit *movingFileOutputName5;
    QLabel *movingFileOutputLabel5;
    QPushButton *movingFileOutputButton5;

    QLabel *options_mode_label;

    QComboBox *options_MetricSelector;
    QLabel *options_MetricSelector_label;

    QLabel *options_registration_label;
    QComboBox *options_registrationSelector;

    QRadioButton *option_registrationMode;
    QRadioButton *option_transformationMode;

    QRadioButton *options_AFFINE_selected;
    QRadioButton *options_RIGID_selected;
    QRadioButton *options_DEFORMABLE_selected;

    QPushButton *confirmButton;
    QPushButton *cancelButton;
    QPushButton *resetButton;


    void setupUi(QDialog *fRegistrationDialog)
    {

        if (fRegistrationDialog->objectName().isEmpty())
            fRegistrationDialog->setObjectName(QString::fromUtf8("fRegistrationDialog"));
        //fRegistrationDialog->setWindowModality(Qt::NonModal);
        fRegistrationDialog->resize(300, 300); // needs to be screenSize dependent 
        QSizePolicy sizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(fRegistrationDialog->sizePolicy().hasHeightForWidth());
        fRegistrationDialog->setSizePolicy(sizePolicy);
        fRegistrationDialog->setMinimumSize(QSize(0, 0));

        /*------------------Layout---------------------------------------*/

        //fRegistrationDialog->setModal(true);
        gridLayout = new QGridLayout(fRegistrationDialog);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));        

        modeGroupBox = new QGroupBox(fRegistrationDialog);
        modeGroupBox->setTitle(QString::fromStdString("Registration Mode:"));

        matrixGroupBox = new QGroupBox(fRegistrationDialog);
        matrixGroupBox->setTitle(QString::fromStdString("Transformation Matrix:"));

        registrationGroupBox = new QGroupBox(fRegistrationDialog);
        registrationGroupBox->setTitle(QString::fromStdString("Options:"));

        inputGroupBox = new QGroupBox(fRegistrationDialog);
        inputGroupBox->setTitle(QString::fromStdString("Input"));

        outputGroupBox = new QGroupBox(fRegistrationDialog);
        outputGroupBox->setTitle(QString::fromStdString("Output"));

        buttonGroupBox = new QGroupBox(fRegistrationDialog);
        buttonGroupBox->setTitle(QString::fromStdString("Actions"));      

        buttonGridLayout = new QHBoxLayout(buttonGroupBox);
        buttonGridLayout->setObjectName(QString::fromUtf8("buttonGridLayout"));  

        modeGridLayout = new QGridLayout(modeGroupBox);
        modeGridLayout->setObjectName(QString::fromUtf8("modeGridLayout"));

        matrixGridLayout = new QGridLayout(matrixGroupBox);
        matrixGridLayout->setObjectName(QString::fromUtf8("matrixGridLayout"));

        registrationGridLayout = new QGridLayout(registrationGroupBox);
        registrationGridLayout->setObjectName(QString::fromUtf8("registrationGridLayout"));

        inputGridLayout = new QGridLayout(inputGroupBox);
        inputGridLayout->setObjectName(QString::fromUtf8("inputGridLayout"));

        outputGridLayout = new QGridLayout(outputGroupBox);
        outputGridLayout->setObjectName(QString::fromUtf8("outputGridLayout"));

        int gridRowCounter = 0;        

        /*-----------------------Registration Mode----------------------*/
        options_registration_label = new QLabel(modeGroupBox);
        options_registration_label->setObjectName("registrationMode");
        options_registration_label->setSizePolicy(sizePolicy);

        options_AFFINE_selected = new QRadioButton(modeGroupBox);
        options_AFFINE_selected->setText("Affine");

        options_RIGID_selected = new QRadioButton(modeGroupBox);
        options_RIGID_selected->setText("Rigid");
        
        options_DEFORMABLE_selected = new QRadioButton(modeGroupBox);
        options_DEFORMABLE_selected->setText("Deformable");

        gridRowCounter++;
        modeGridLayout->addWidget(options_registration_label, gridRowCounter, 0, 1, 1);
        modeGridLayout->addWidget(options_AFFINE_selected, gridRowCounter, 1, 1, 1);
        modeGridLayout->addWidget(options_RIGID_selected, gridRowCounter, 2, 1, 1);
        modeGridLayout->addWidget(options_DEFORMABLE_selected, gridRowCounter, 3, 1, 1);

        /*---------------------Metrics---------------------*/
        options_MetricSelector_label = new QLabel(registrationGroupBox);
        options_MetricSelector_label->setObjectName("metric");
        options_MetricSelector_label->setSizePolicy(sizePolicy);
        options_MetricSelector = new QComboBox(registrationGroupBox);
        QSizePolicy sizePolicy8(QSizePolicy::Minimum, QSizePolicy::Fixed);
        options_MetricSelector->setSizePolicy(sizePolicy);
        //options_T1Selector->setMaximumSize(QSize(100, 1000));
        options_MetricSelector->insertItem(0, "NMI(Default)");
        options_MetricSelector->insertItem(1, "MI");
        options_MetricSelector->insertItem(2, "NCC");
        options_MetricSelector->insertItem(3, "SSD");
        options_MetricSelector->setToolTip("Set the type Metric for image registration");

        gridRowCounter++;
        registrationGridLayout->addWidget(options_MetricSelector_label, gridRowCounter, 0, 1, 1);
        registrationGridLayout->addWidget(options_MetricSelector, gridRowCounter, 1, 1, 1);

        /*------------------NCC Radius ------------------*/
        options_NCC_radii = new QLabel(registrationGroupBox);
        options_NCC_radii->setObjectName("nccRadiusLabel");
        options_NCC_radii->setSizePolicy(sizePolicy);
        QRegExp rx("^([0-9]|[1-9][0-9]|100)x([0-9]|[1-9][0-9]|100)x([0-9]|[1-9][0-9]|100)");
        QValidator *validator = new QRegExpValidator(rx);
        nccRadii = new QLineEdit("5x5x5");
        nccRadii->setObjectName("nccRadius");
        nccRadii->setValidator(validator);
        nccRadii->setToolTip("Patch Radius (Eg: 2x2x2)");
        nccRadii->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        gridRowCounter++;
        registrationGridLayout->addWidget(options_NCC_radii, gridRowCounter, 0, 1, 1);
        registrationGridLayout->addWidget(nccRadii, gridRowCounter, 1, 1, 1);

        /*------------------Iterations--------------------------------*/
        options_iterations = new QLabel(registrationGroupBox);
        options_iterations->setObjectName("iterationLabel");
        options_iterations->setSizePolicy(sizePolicy);
        QRegExp ra("^([0-9]|[1-9][0-9]|100)x([0-9]|[1-9][0-9]|100)x([0-9]|[1-9][0-9]|100)");
        QValidator *val = new QRegExpValidator(ra);
        iterations = new QLineEdit("100x50x5");        
        iterations->setObjectName("iterations");
        iterations->setValidator(val); 
        iterations->setToolTip("Iterations per resolution: Low x Medium x High");
        iterations->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        gridRowCounter++;
        registrationGridLayout->addWidget(options_iterations, gridRowCounter, 0, 1, 1);
        registrationGridLayout->addWidget(iterations, gridRowCounter, 1, 1, 1);

        /*--------------------Fixed Image--------------------------------*/
        fixedFileLabel = new QLabel(registrationGroupBox);
        sizePolicy.setHeightForWidth(fixedFileLabel->sizePolicy().hasHeightForWidth());
        fixedFileLabel->setSizePolicy(sizePolicy);

        fixedFileName = new QLineEdit(registrationGroupBox);
        fixedFileName->setObjectName(QString::fromUtf8("Fixed File"));
        fixedFileName->setToolTip("Fixed Image file. Eg: fixed.nii.gz");
        sizePolicy.setHeightForWidth(fixedFileName->sizePolicy().hasHeightForWidth());
        fixedFileName->setSizePolicy(sizePolicy);
        fixedFileName->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        fixedFileButton = new QPushButton(registrationGroupBox);
        fixedFileButton->setObjectName(QString::fromUtf8("fixedFileButton"));
        fixedFileButton->setText("Browse");

        gridRowCounter++;
        registrationGridLayout->addWidget(fixedFileLabel, gridRowCounter, 0, 1, 1);
        registrationGridLayout->addWidget(fixedFileName, gridRowCounter, 1, 1, 1);
        registrationGridLayout->addWidget(fixedFileButton, gridRowCounter, 2, 1, 1);

        /*----------------------Matrix checkbox--------------------*/
        matrixRadioLabel = new QLabel(registrationGroupBox);
        sizePolicy.setHeightForWidth(matrixRadioLabel->sizePolicy().hasHeightForWidth());
        matrixRadioLabel->setSizePolicy(sizePolicy);

        matrixRadioButton = new QCheckBox(registrationGroupBox);
        matrixRadioButton->setObjectName(QString::fromUtf8("matrixRadioButton"));
        matrixRadioButton->setToolTip("Generate transformation matrix");
        sizePolicy.setHeightForWidth(matrixRadioButton->sizePolicy().hasHeightForWidth());
        matrixRadioButton->setSizePolicy(sizePolicy);

        gridRowCounter++;

        registrationGridLayout->addWidget(matrixRadioLabel, gridRowCounter, 0, 1, 1);
        registrationGridLayout->addWidget(matrixRadioButton, gridRowCounter, 1, 1, 1);


        /*---------------------Add more button-------------------*/
        addMoreButton = new QPushButton(inputGroupBox);
        addMoreButton->setObjectName(QString::fromUtf8("addMoreButton"));
        addMoreButton->setText("+");
        addMoreButton->setMaximumWidth(30);

        gridRowCounter++;

        inputGridLayout->addWidget(addMoreButton, gridRowCounter, 0, 1, 1);


        /*---------------------Moving Image 1-------------------*/
        movingFileLabel1 = new QLabel(inputGroupBox);
        sizePolicy.setHeightForWidth(movingFileLabel1->sizePolicy().hasHeightForWidth());
        movingFileLabel1->setSizePolicy(sizePolicy);

        movingFileName1 = new QLineEdit(inputGroupBox);
        movingFileName1->setObjectName(QString::fromUtf8("movingInputImage"));
        sizePolicy.setHeightForWidth(movingFileName1->sizePolicy().hasHeightForWidth());
        movingFileName1->setSizePolicy(sizePolicy);
        movingFileName1->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileButton1 = new QPushButton(inputGroupBox);
        movingFileButton1->setObjectName(QString::fromUtf8("movingInputFileButton"));
        movingFileButton1->setText("Browse");        

        gridRowCounter++;

        inputGridLayout->addWidget(movingFileLabel1, gridRowCounter, 0, 1, 1);
        inputGridLayout->addWidget(movingFileName1, gridRowCounter, 1, 1, 1);
        inputGridLayout->addWidget(movingFileButton1, gridRowCounter, 2, 1, 1);        


        /*---------------------Moving Image 2-------------------*/
        movingFileLabel2 = new QLabel(inputGroupBox);
        sizePolicy.setHeightForWidth(movingFileLabel2->sizePolicy().hasHeightForWidth());
        movingFileLabel2->setSizePolicy(sizePolicy);

        movingFileName2 = new QLineEdit(inputGroupBox);
        movingFileName2->setObjectName(QString::fromUtf8("movingInputImage"));
        sizePolicy.setHeightForWidth(movingFileName2->sizePolicy().hasHeightForWidth());
        movingFileName2->setSizePolicy(sizePolicy);
        movingFileName2->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileButton2 = new QPushButton(inputGroupBox);
        movingFileButton2->setObjectName(QString::fromUtf8("movingInputFileButton"));
        movingFileButton2->setText("Browse");

        gridRowCounter++;

        inputGridLayout->addWidget(movingFileLabel2, gridRowCounter, 0, 1, 1);
        inputGridLayout->addWidget(movingFileName2, gridRowCounter, 1, 1, 1);
        inputGridLayout->addWidget(movingFileButton2, gridRowCounter, 2, 1, 1);


        /*---------------------Moving Image 3-------------------*/
        movingFileLabel3 = new QLabel(inputGroupBox);
        sizePolicy.setHeightForWidth(movingFileLabel3->sizePolicy().hasHeightForWidth());
        movingFileLabel3->setSizePolicy(sizePolicy);

        movingFileName3 = new QLineEdit(inputGroupBox);
        movingFileName3->setObjectName(QString::fromUtf8("movingInputImage"));
        sizePolicy.setHeightForWidth(movingFileName3->sizePolicy().hasHeightForWidth());
        movingFileName3->setSizePolicy(sizePolicy);
        movingFileName3->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileButton3 = new QPushButton(inputGroupBox);
        movingFileButton3->setObjectName(QString::fromUtf8("movingInputFileButton"));
        movingFileButton3->setText("Browse");

        gridRowCounter++;

        inputGridLayout->addWidget(movingFileLabel3, gridRowCounter, 0, 1, 1);
        inputGridLayout->addWidget(movingFileName3, gridRowCounter, 1, 1, 1);
        inputGridLayout->addWidget(movingFileButton3, gridRowCounter, 2, 1, 1);

        /*---------------------Moving Image 4-------------------*/
        movingFileLabel4 = new QLabel(inputGroupBox);
        sizePolicy.setHeightForWidth(movingFileLabel4->sizePolicy().hasHeightForWidth());
        movingFileLabel4->setSizePolicy(sizePolicy);

        movingFileName4 = new QLineEdit(inputGroupBox);
        movingFileName4->setObjectName(QString::fromUtf8("movingInputImage"));
        sizePolicy.setHeightForWidth(movingFileName4->sizePolicy().hasHeightForWidth());
        movingFileName4->setSizePolicy(sizePolicy);
        movingFileName4->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileButton4 = new QPushButton(inputGroupBox);
        movingFileButton4->setObjectName(QString::fromUtf8("movingInputFileButton"));
        movingFileButton4->setText("Browse");

        gridRowCounter++;

        inputGridLayout->addWidget(movingFileLabel4, gridRowCounter, 0, 1, 1);
        inputGridLayout->addWidget(movingFileName4, gridRowCounter, 1, 1, 1);
        inputGridLayout->addWidget(movingFileButton4, gridRowCounter, 2, 1, 1);

        /*---------------------Moving Image 5-------------------*/
        movingFileLabel5 = new QLabel(inputGroupBox);
        sizePolicy.setHeightForWidth(movingFileLabel5->sizePolicy().hasHeightForWidth());
        movingFileLabel5->setSizePolicy(sizePolicy);

        movingFileName5 = new QLineEdit(inputGroupBox);
        movingFileName5->setObjectName(QString::fromUtf8("movingInputImage"));
        sizePolicy.setHeightForWidth(movingFileName5->sizePolicy().hasHeightForWidth());
        movingFileName5->setSizePolicy(sizePolicy);
        movingFileName5->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileButton5 = new QPushButton(inputGroupBox);
        movingFileButton5->setObjectName(QString::fromUtf8("movingInputFileButton"));
        movingFileButton5->setText("Browse");

        gridRowCounter++;

        inputGridLayout->addWidget(movingFileLabel5, gridRowCounter, 0, 1, 1);
        inputGridLayout->addWidget(movingFileName5, gridRowCounter, 1, 1, 1);
        inputGridLayout->addWidget(movingFileButton5, gridRowCounter, 2, 1, 1);
       
        /*----------------------Matrix 1----------------------*/
        matrixLabel1 = new QLabel(matrixGroupBox);
        sizePolicy.setHeightForWidth(matrixLabel1->sizePolicy().hasHeightForWidth());
        matrixLabel1->setSizePolicy(sizePolicy);

        matrixName1 = new QLineEdit(matrixGroupBox);
        matrixName1->setObjectName(QString::fromUtf8("matrixFileName1"));
        sizePolicy.setHeightForWidth(matrixName1->sizePolicy().hasHeightForWidth());
        matrixName1->setSizePolicy(sizePolicy);
        matrixName1->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        matrixButton1 = new QPushButton(matrixGroupBox);
        matrixButton1->setObjectName(QString::fromUtf8("matrixButtonName1"));
        matrixButton1->setText(QString("Browse"));

        gridRowCounter++;

        matrixGridLayout->addWidget(matrixLabel1, gridRowCounter, 0, 1, 1);
        matrixGridLayout->addWidget(matrixName1, gridRowCounter, 1, 1, 1);
        matrixGridLayout->addWidget(matrixButton1, gridRowCounter, 2, 1, 1);

        /*----------------------Matrix 2----------------------*/
        matrixLabel2 = new QLabel(matrixGroupBox);
        sizePolicy.setHeightForWidth(matrixLabel2->sizePolicy().hasHeightForWidth());
        matrixLabel2->setSizePolicy(sizePolicy);

        matrixName2 = new QLineEdit(matrixGroupBox);
        matrixName2->setObjectName(QString::fromUtf8("matrixFileName2"));
        sizePolicy.setHeightForWidth(matrixName2->sizePolicy().hasHeightForWidth());
        matrixName2->setSizePolicy(sizePolicy);
        matrixName2->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        matrixButton2 = new QPushButton(matrixGroupBox);
        matrixButton2->setObjectName(QString::fromUtf8("matrixButtonName2"));
        matrixButton2->setText(QString("Browse"));
        gridRowCounter++;

        matrixGridLayout->addWidget(matrixLabel2, gridRowCounter, 0, 1, 1);
        matrixGridLayout->addWidget(matrixName2, gridRowCounter, 1, 1, 1);
        matrixGridLayout->addWidget(matrixButton2, gridRowCounter, 2, 1, 1);

        /*----------------------Matrix 3----------------------*/
        matrixLabel3 = new QLabel(matrixGroupBox);
        sizePolicy.setHeightForWidth(matrixLabel3->sizePolicy().hasHeightForWidth());
        matrixLabel3->setSizePolicy(sizePolicy);

        matrixName3 = new QLineEdit(matrixGroupBox);
        matrixName3->setObjectName(QString::fromUtf8("matrixFileName3"));
        sizePolicy.setHeightForWidth(matrixName3->sizePolicy().hasHeightForWidth());
        matrixName3->setSizePolicy(sizePolicy);
        matrixName3->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        matrixButton3 = new QPushButton(matrixGroupBox);
        matrixButton3->setObjectName(QString::fromUtf8("matrixButtonName3"));
        matrixButton3->setText(QString("Browse"));
        gridRowCounter++;

        matrixGridLayout->addWidget(matrixLabel3, gridRowCounter, 0, 1, 1);
        matrixGridLayout->addWidget(matrixName3, gridRowCounter, 1, 1, 1);
        matrixGridLayout->addWidget(matrixButton3, gridRowCounter, 2, 1, 1);

        /*----------------------Matrix 4----------------------*/
        matrixLabel4 = new QLabel(matrixGroupBox);
        sizePolicy.setHeightForWidth(matrixLabel4->sizePolicy().hasHeightForWidth());
        matrixLabel4->setSizePolicy(sizePolicy);

        matrixName4 = new QLineEdit(matrixGroupBox);
        matrixName4->setObjectName(QString::fromUtf8("matrixFileName4"));
        sizePolicy.setHeightForWidth(matrixName4->sizePolicy().hasHeightForWidth());
        matrixName4->setSizePolicy(sizePolicy);
        matrixName4->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        matrixButton4 = new QPushButton(matrixGroupBox);
        matrixButton4->setObjectName(QString::fromUtf8("matrixButtonName3"));
        matrixButton4->setText(QString("Browse"));
        gridRowCounter++;

        matrixGridLayout->addWidget(matrixLabel4, gridRowCounter, 0, 1, 1);
        matrixGridLayout->addWidget(matrixName4, gridRowCounter, 1, 1, 1);
        matrixGridLayout->addWidget(matrixButton4, gridRowCounter, 2, 1, 1);

        /*----------------------Matrix 5----------------------*/
        matrixLabel5 = new QLabel(matrixGroupBox);
        sizePolicy.setHeightForWidth(matrixLabel5->sizePolicy().hasHeightForWidth());
        matrixLabel5->setSizePolicy(sizePolicy);

        matrixName5 = new QLineEdit(matrixGroupBox);
        matrixName5->setObjectName(QString::fromUtf8("matrixFileName5"));
        sizePolicy.setHeightForWidth(matrixName5->sizePolicy().hasHeightForWidth());
        matrixName5->setSizePolicy(sizePolicy);
        matrixName5->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        matrixButton5 = new QPushButton(matrixGroupBox);
        matrixButton5->setObjectName(QString::fromUtf8("matrixButtonName5"));
        matrixButton5->setText(QString("Browse"));
        gridRowCounter++;

        matrixGridLayout->addWidget(matrixLabel5, gridRowCounter, 0, 1, 1);
        matrixGridLayout->addWidget(matrixName5, gridRowCounter, 1, 1, 1);
        matrixGridLayout->addWidget(matrixButton5, gridRowCounter, 2, 1, 1);        

        /*----------------------Output Image 1----------------------*/
        movingFileOutputLabel1 = new QLabel(outputGroupBox);
        sizePolicy.setHeightForWidth(movingFileOutputLabel1->sizePolicy().hasHeightForWidth());
        movingFileOutputLabel1->setSizePolicy(sizePolicy);

        movingFileOutputName1 = new QLineEdit(outputGroupBox);
        movingFileOutputName1->setObjectName(QString::fromUtf8("outputImageName1"));
        sizePolicy.setHeightForWidth(movingFileOutputName1->sizePolicy().hasHeightForWidth());
        movingFileOutputName1->setSizePolicy(sizePolicy);
        movingFileOutputName1->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileOutputButton1 = new QPushButton(outputGroupBox);
        movingFileOutputButton1->setObjectName(QString::fromUtf8("outputImageButton1"));
        movingFileOutputButton1->setText(QString("Browse"));

        gridRowCounter++;

        outputGridLayout->addWidget(movingFileOutputLabel1, gridRowCounter, 0, 1, 1);
        outputGridLayout->addWidget(movingFileOutputName1, gridRowCounter, 1, 1, 1);
        outputGridLayout->addWidget(movingFileOutputButton1, gridRowCounter, 2, 1, 1);

        /*----------------------Output Image 2----------------------*/
        movingFileOutputLabel2 = new QLabel(outputGroupBox);
        sizePolicy.setHeightForWidth(movingFileOutputLabel2->sizePolicy().hasHeightForWidth());
        movingFileOutputLabel2->setSizePolicy(sizePolicy);

        movingFileOutputName2 = new QLineEdit(outputGroupBox);
        movingFileOutputName2->setObjectName(QString::fromUtf8("outputImageName2"));
        sizePolicy.setHeightForWidth(movingFileOutputName2->sizePolicy().hasHeightForWidth());
        movingFileOutputName2->setSizePolicy(sizePolicy);
        movingFileOutputName2->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileOutputButton2 = new QPushButton(outputGroupBox);
        movingFileOutputButton2->setObjectName(QString::fromUtf8("outputImageButton2"));
        movingFileOutputButton2->setText(QString("Browse"));
        gridRowCounter++;

        outputGridLayout->addWidget(movingFileOutputLabel2, gridRowCounter, 0, 1, 1);
        outputGridLayout->addWidget(movingFileOutputName2, gridRowCounter, 1, 1, 1);
        outputGridLayout->addWidget(movingFileOutputButton2, gridRowCounter, 2, 1, 1);

        /*----------------------Output Image 3----------------------*/
        movingFileOutputLabel3 = new QLabel(outputGroupBox);
        sizePolicy.setHeightForWidth(movingFileOutputLabel3->sizePolicy().hasHeightForWidth());
        movingFileOutputLabel3->setSizePolicy(sizePolicy);

        movingFileOutputName3 = new QLineEdit(outputGroupBox);
        movingFileOutputName3->setObjectName(QString::fromUtf8("outputImageName3"));
        sizePolicy.setHeightForWidth(movingFileOutputName3->sizePolicy().hasHeightForWidth());
        movingFileOutputName3->setSizePolicy(sizePolicy);
        movingFileOutputName3->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileOutputButton3 = new QPushButton(outputGroupBox);
        movingFileOutputButton3->setObjectName(QString::fromUtf8("outputImageButton3"));
        movingFileOutputButton3->setText(QString("Browse"));
        gridRowCounter++;

        outputGridLayout->addWidget(movingFileOutputLabel3, gridRowCounter, 0, 1, 1);
        outputGridLayout->addWidget(movingFileOutputName3, gridRowCounter, 1, 1, 1);
        outputGridLayout->addWidget(movingFileOutputButton3, gridRowCounter, 2, 1, 1);

        /*----------------------Output Image 4----------------------*/
        movingFileOutputLabel4 = new QLabel(outputGroupBox);
        sizePolicy.setHeightForWidth(movingFileOutputLabel4->sizePolicy().hasHeightForWidth());
        movingFileOutputLabel4->setSizePolicy(sizePolicy);

        movingFileOutputName4 = new QLineEdit(outputGroupBox);
        movingFileOutputName4->setObjectName(QString::fromUtf8("outputImageName4"));
        sizePolicy.setHeightForWidth(movingFileOutputName4->sizePolicy().hasHeightForWidth());
        movingFileOutputName4->setSizePolicy(sizePolicy);
        movingFileOutputName4->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileOutputButton4 = new QPushButton(outputGroupBox);
        movingFileOutputButton4->setObjectName(QString::fromUtf8("outputImageButton3"));
        movingFileOutputButton4->setText(QString("Browse"));
        gridRowCounter++;

        outputGridLayout->addWidget(movingFileOutputLabel4, gridRowCounter, 0, 1, 1);
        outputGridLayout->addWidget(movingFileOutputName4, gridRowCounter, 1, 1, 1);
        outputGridLayout->addWidget(movingFileOutputButton4, gridRowCounter, 2, 1, 1);

        /*----------------------Output Image 5----------------------*/
        movingFileOutputLabel5 = new QLabel(outputGroupBox);
        sizePolicy.setHeightForWidth(movingFileOutputLabel5->sizePolicy().hasHeightForWidth());
        movingFileOutputLabel5->setSizePolicy(sizePolicy);

        movingFileOutputName5 = new QLineEdit(outputGroupBox);
        movingFileOutputName5->setObjectName(QString::fromUtf8("outputImageName5"));
        sizePolicy.setHeightForWidth(movingFileOutputName5->sizePolicy().hasHeightForWidth());
        movingFileOutputName5->setSizePolicy(sizePolicy);
        movingFileOutputName5->setAlignment(Qt::AlignCenter | Qt::AlignTrailing | Qt::AlignVCenter);

        movingFileOutputButton5 = new QPushButton(outputGroupBox);
        movingFileOutputButton5->setObjectName(QString::fromUtf8("outputImageButton5"));
        movingFileOutputButton5->setText(QString("Browse"));
        gridRowCounter++;

        outputGridLayout->addWidget(movingFileOutputLabel5, gridRowCounter, 0, 1, 1);
        outputGridLayout->addWidget(movingFileOutputName5, gridRowCounter, 1, 1, 1);
        outputGridLayout->addWidget(movingFileOutputButton5, gridRowCounter, 2, 1, 1);

        // put the layout in perspective        
        gridLayout->addWidget(modeGroupBox, 0, 0, 1, 2);
        gridLayout->addWidget(registrationGroupBox, 1, 0, 1, 2);        
        gridLayout->addWidget(inputGroupBox, 2, 0, 1, 2);
        gridLayout->addWidget(matrixGroupBox, 3, 0, 1, 2);
        gridLayout->addWidget(outputGroupBox, 4, 0, 1, 2);
        gridLayout->addWidget(buttonGroupBox, 5, 0, 1, 2);

        /*------------------------Buttons-------------------------------*/
        confirmButton = new QPushButton(buttonGroupBox);
        confirmButton->setObjectName(QString::fromUtf8("confirm"));
        //confirmButton->setIcon(ButtonIcon);
        //confirmButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
        confirmButton->setText(QString("Register"));

        cancelButton = new QPushButton(buttonGroupBox);
        cancelButton->setObjectName(QString::fromUtf8("Cancel"));
        ////cancelButton->setIcon(ButtonIcon);
        ////cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
        cancelButton->setText(QString("Cancel"));

        resetButton = new QPushButton(buttonGroupBox);
        resetButton->setObjectName(QString::fromUtf8("Reset"));
        ////cancelButton->setIcon(ButtonIcon);
        ////cancelButton->setIconSize(QSize(20, 20)); // needs to be screenSize dependent 
        resetButton->setText(QString("Reset"));

        buttonGridLayout->addWidget(confirmButton);
        buttonGridLayout->addWidget(resetButton);
        buttonGridLayout->addWidget(cancelButton);

        retranslateUi(fRegistrationDialog);

        QMetaObject::connectSlotsByName(fRegistrationDialog);
    } // setupUi

    void retranslateUi(QDialog *fRegistrationDialog) {

        fRegistrationDialog->setWindowTitle(QApplication::translate("fRegistrationDialog", "Registration", 0));
        matrixRadioButton->setToolTip("Generate transformation matix");
        movingFileButton1->setToolTip("First moving image.");
        movingFileButton2->setToolTip("Second moving image.");
        movingFileButton3->setToolTip("Third moving image.");
        movingFileButton4->setToolTip("Fourth moving image.");
        movingFileButton5->setToolTip("Fifth moving image.");
        movingFileOutputButton1->setToolTip("First output image.");
        movingFileOutputButton2->setToolTip("Second output image.");
        movingFileOutputButton3->setToolTip("Third output image.");
        movingFileOutputButton4->setToolTip("Fourth output image.");
        movingFileOutputButton5->setToolTip("Fifth output image.");        
        matrixButton1->setToolTip("First transformation matrix.");
        matrixButton2->setToolTip("Second transformation matrix.");
        matrixButton3->setToolTip("Third transformation matrix.");
        matrixButton4->setToolTip("Fourth transformation matrix.");
        matrixButton5->setToolTip("Fifth transformation matrix.");
        options_MetricSelector_label->setText("Metrics:");        
        fixedFileLabel->setText("Fixed Image:");
        options_registration_label->setText("Registration:");
        options_NCC_radii->setText("Enter radii:");
        options_iterations->setText("Iterations:");
        movingFileLabel1->setText("Moving Image 1: ");
        movingFileLabel2->setText("Moving Image 2:");
        movingFileLabel3->setText("Moving Image 3:");
        movingFileLabel4->setText("Moving Image 4:");
        movingFileLabel5->setText("Moving Image 5:");
        matrixRadioLabel->setText("Generate Matrix");
        matrixLabel1->setText("Matrix 1:");
        matrixLabel2->setText("Matrix 2:");
        matrixLabel3->setText("Matrix 3:");
        matrixLabel4->setText("Matrix 4:");
        matrixLabel5->setText("Matrix 5:");
        movingFileOutputLabel1->setText("Output Image 1:");
        movingFileOutputLabel2->setText("Output Image 2:");
        movingFileOutputLabel3->setText("Output Image 3:");
        movingFileOutputLabel4->setText("Output Image 4:");
        movingFileOutputLabel5->setText("Output Image 5:");

    }
    // retranslateUi

};

namespace Ui {
    class fRegistrationDialog : public ui_fRegistrationDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // ui_fRegistrationDialog_H




