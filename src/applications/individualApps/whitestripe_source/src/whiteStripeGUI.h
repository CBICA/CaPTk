/**
\file  whiteStripeGUI.h

\brief Declaration of the WhiteStripe GUI class

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2017 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html

*/

#ifndef WhiteStripeGUI_H
#define WhiteStripeGUI_H
#include <QMainWindow>
#include <QtGui>
#include <QWidget>
#include <QVTKWidget.h>
#include <itkImageToVTKImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkN4BiasFieldCorrectionImageFilter.h>
#include "itkOtsuThresholdImageFilter.h"

#include <vtkImageViewer2.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>

#include <vtkLookupTable.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkThresholdPoints.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h> 
#include <vtkImageMapToColors.h>
#include <vtkProperty.h> 

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkImageRegistrationMethod.h>

#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkPen.h>
#include <vtkLookupTable.h>
#include <vtkInteractorStyleImage.h>
#include <vtkInteractorStyleRubberBandZoom.h>
#include <vtkImageShiftScale.h>
#include <vtkImageActor.h>

#include "utils.h"
#include "whiteStripe.h"
/**
\class HistWidget
\brief a simple QWidget based widget for displaying  histogram
*/
class HistWidget : public QWidget
{
	Q_OBJECT
private: 
	QVTKWidget* widget;
	vtkSmartPointer<vtkContextView> m_view;
	vtkSmartPointer<vtkTable> m_table;
	vtkSmartPointer<vtkFloatArray> m_arrX;
	vtkSmartPointer<vtkChartXY> m_chart;

public:
	explicit HistWidget(QWidget *parent): QWidget(parent, Qt::Window)
	{

		widget = new QVTKWidget(this);
		widget->setMinimumSize(QSize(256, 256));
		m_view = vtkSmartPointer<vtkContextView>::New();
		m_view->SetInteractor(widget->GetInteractor());
		widget->SetRenderWindow(m_view->GetRenderWindow());
		m_view->GetRenderer()->SetBackground(0.3, 0.3, 0.3);
		m_view->GetInteractor()->Initialize();
		QHBoxLayout *mainLayout = new QHBoxLayout;
		mainLayout->addWidget(widget);
		setLayout(mainLayout);
		setWindowTitle("Histogram visualization");
		m_table = vtkSmartPointer<vtkTable>::New();
		m_chart = vtkSmartPointer<vtkChartXY>::New();
		m_view->GetScene()->AddItem(m_chart);

	}
	void plotVirticalLine(float x=0.0f, float height=10.0f, string name= "Line")
	{
		// Create a table with some points in it
		vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
		vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
		arrX->SetName("X Axis");
		table->AddColumn(arrX);
		vtkSmartPointer<vtkFloatArray> arrC = vtkSmartPointer<vtkFloatArray>::New();
		arrC->SetName(name.c_str());
		table->AddColumn(arrC);


		table->SetNumberOfRows(2);
		table->SetValue(0, 0, x);
		table->SetValue(0, 1, 0);
		table->SetValue(1, 0, x);
		table->SetValue(1, 1, height);


		vtkPlot *line = m_chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION >= 6
		line->SetInputData(table, 0, 1);
#else
		line->SetInput(table, 0, 1);
#endif
		line->SetColor(255, 0, 0, 255);
		line->SetWidth(2.0);
	
	}
	void setAxis(const vector<float>& axis, const int colCount, const string& axisName = "X Axis")
	{
		vtkSmartPointer<vtkFloatArray> arrayAxis = vtkSmartPointer<vtkFloatArray>::New();
		arrayAxis->SetName(axisName.c_str());
		m_table->AddColumn(arrayAxis);
    for (size_t i = 0; i < (uint)colCount; i++)
		{
			m_table->AddColumn(vtkSmartPointer<vtkFloatArray>::New());
		}
		m_table->SetNumberOfRows(axis.size());
    for (size_t i = 0; i < axis.size(); ++i)
		{
			m_table->SetValue(i, 0, axis[i]);
		}
		return;
	}
	void addColumn(const vector<float>& curve, const string& name, int colId, Scalar color = Scalar(0,255,0,255))
	{
		if (colId >= m_table->GetNumberOfColumns())
		{
			return;//Axis is not set properly
		}
		m_table->GetColumn(colId)->SetName(name.c_str());
		for (size_t i = 0; i < curve.size(); ++i)
		{
			m_table->SetValue(i, colId, curve[i]);
		}
		vtkPlot *linePlot = m_chart->AddPlot(vtkChart::LINE);
#if VTK_MAJOR_VERSION >= 6
		linePlot->SetInputData(m_table, 0, colId);
#else
		linePlot->SetInput(m_table, 0, colId);
#endif
		linePlot->SetColor(color.val[0], color.val[1], color.val[2], color.val[3]);
		linePlot->SetWidth(1.0);
	}
};

/**
\class WhiteStripeGUI
\brief a quick UI prototype for using WhiteStripe algorithm
Reference:
@article{Shinohara20149,
title = "Statistical normalization techniques for magnetic resonance imaging ",
journal = "NeuroImage: Clinical ",
volume = "6",
pages = "9 - 19",
year = "2014",
issn = "2213-1582",
}
*/
class WhiteStripeGUI : public QMainWindow {
    Q_OBJECT
public:
	WhiteStripeGUI(QWidget *parent) : QMainWindow(parent)
	{
		createMenus();
		QHBoxLayout* mainLayout = new QHBoxLayout();
		mainLayout->addLayout(this->createViews());
		mainLayout->addLayout(this->creatSideBar());

		QWidget *centralWidget = new QWidget();
		centralWidget->setLayout(mainLayout);
		setCentralWidget(centralWidget);
		this->setWindowTitle(tr("WhiteStripe"));
		this->log("WhiteStripe� launch Success!");

		m_hWdg = NULL;
		m_vtkImg = NULL;
		m_vtkImgMask = NULL;
		m_itkImg = NULL;
    m_itkImg_mask = NULL;

	}


private slots:
  void batchItemDoubleClicked(QListWidgetItem* item);
  void updateDisplay();
  void open(QString fileName = QString());
  void addBatchFiles();
  void save();
  void imgSliderValuechanged(int val);
  void runAlg();
  void runCore();
  void runBatch();
  void displayHistChart();
  void displayAboutAct();
  void unloadImage();
  void updateWindowLevel();
  void toggleMaskClicked();
  void setWindowLevel();
private:
  vtkSmartPointer<vtkImageData> toVtkImg(ImageType::Pointer image);
  vtkImageData* getDisplayableVtkImg(vtkImageData* vtkImg, bool mask = true);
  void displayMask(vtkImageData* image);
  void updateDisplay3D();
  void openMri();
  void createMenus();
  void CreateDummyImage(vtkSmartPointer<vtkImageData>& image);
  QGridLayout* createViews();
  QVBoxLayout* creatSideBar();
  QScrollArea* createBatchTab(const int lblWidth);
  QScrollArea* createAlgTab(const int lblWidth);
  float  parseProgress(string original);
  void log(string msg, bool  bold = false, string colour = "");
  void maxDblClick(size_t index);
  bool eventFilter(QObject *obj, QEvent *event);

	//Settings fields 
	QLineEdit* m_txt_inFile;
	QLineEdit* m_txt_outFile;
	QLineEdit* m_txt_twsWidth;
	QLineEdit* m_txt_sliceStartZ;
	QLineEdit* m_txt_sliceStopZ;
	QLineEdit* m_txt_tissuesMax;
	QLineEdit* m_txt_smoothMax;
	QLineEdit* m_txt_smoothDelta;
	QLineEdit* m_txt_histSize;
	QComboBox* m_cmb_bT1;

	double m_window = 255;
	double m_level = 127.5;

	//Member GUI elements
	QListWidget* m_batchFileList;
	QLineEdit* m_BatchAppendName;
	QGridLayout* m_gridLayout;
	QLabel* m_progressLable;
	QTextEdit* m_logOutput;
	QProgressBar* m_progressBar;
	QCheckBox* m_toggleMask;

	HistWidget* m_hWdg;
	WhiteStripe m_alg;
	vtkSmartPointer<vtkImageData> m_vtkImg;
	vtkSmartPointer<vtkImageData> m_vtkImgMask;
  ImageType::Pointer m_itkImg, m_itkImg_mask;

	vector<QWidget*> m_views;
	vector<vtkImageViewer2*> m_sliceViewers;
	vector<vtkRenderer*> m_renderers;
	vector<QSlider*> m_sliders;
};

#endif
