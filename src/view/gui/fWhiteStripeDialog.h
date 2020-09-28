///////////////////////////////////////////////////////////////////////////////////////
// fWhiteStripeObj.h
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
// License Agreement: https://www.med.upenn.edu/cbica/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#ifndef _fWhiteStripeDialog_h_
#define _fWhiteStripeDialog_h_


//#include "CAPTk.h"
#include "ui_fWhiteStripeDialog.h"

#include "QWidget"
#include "QVTKOpenGLWidget.h"
#include <QScopedPointer>
#include "vtkSmartPointer.h"
#include "vtkChartXY.h"
#include "vtkTable.h"
#include "vtkPlot.h"
#include "vtkFloatArray.h"
#include "vtkContextView.h"
#include "vtkContextScene.h"
#include "vtkRenderer.h"
#include "vtkGenericOpenGLRenderWindow.h"

#include "opencv2/core/core.hpp"

#include "CaPTkEnums.h"

/**
\class HistWidget
\brief a simple QWidget based widget for displaying  histogram
*/
class HistWidget : public QWidget
{
  Q_OBJECT
private:
  QScopedPointer<QVTKOpenGLWidget> widget;
  vtkSmartPointer<vtkContextView> m_view;
  vtkSmartPointer<vtkTable> m_table;
  vtkSmartPointer<vtkFloatArray> m_arrX;
  vtkSmartPointer<vtkChartXY> m_chart;

public:
  explicit HistWidget(QWidget *parent) : QWidget(parent, Qt::Window)
  {
    widget.reset(new QVTKOpenGLWidget(this));
    widget->setMinimumSize(QSize(256, 256));
    m_view = vtkSmartPointer<vtkContextView>::New();
    auto m_renWin = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    auto m_iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    m_renWin->SetInteractor(m_iren);
    widget->SetRenderWindow(m_renWin);
    m_view->SetInteractor(widget->GetInteractor()); // crash is happening here with reference to QVTKOpenGLWidget
    widget->SetRenderWindow(m_view->GetRenderWindow());
    m_view->GetRenderer()->SetBackground(0.3, 0.3, 0.3);
    m_view->GetInteractor()->Initialize();
    QHBoxLayout *mainLayout = new QHBoxLayout;
    mainLayout->addWidget(widget.data());
    setLayout(mainLayout);
    setWindowTitle("Histogram visualization");
    m_table = vtkSmartPointer<vtkTable>::New();
    m_chart = vtkSmartPointer<vtkChartXY>::New();
    m_view->GetScene()->AddItem(m_chart);

  }
  void plotVerticalLine(float x = 0.0f, float height = 10.0f, std::string name = "Line")
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
  void setAxis(const std::vector<float>& axis, const int colCount, const std::string& axisName = "X Axis")
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
  void addColumn(const std::vector<float>& curve, const std::string& name, int colId, cv::Scalar color = cv::Scalar(0, 255, 0, 255))
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
\class fWhiteStripeObj

\brief This class controls the elements in the WhiteStripe dialog
*/
class fWhiteStripeObj : public QDialog, private Ui::fWhiteStripeObj
{
  Q_OBJECT

public:
  fWhiteStripeObj();
  ~fWhiteStripeObj();
  int mode;
  bool skullStrippedImage = true, t1Image = true;

  void SetCurrentImagePath(const QString &inputPath)
  {
    mInputPathName = inputPath;
    outputImageName->setText(mInputPathName + "whiteStripe_output.nii.gz");
  }

  void SetImageModality(int modality)
  {
    switch (modality)
    {
    case CAPTK::ImageModalityType::IMAGE_TYPE_T1:
      //options_T1selected->setChecked(true);
      //options_T2selected->setChecked(false);
      options_T1Selector->setCurrentIndex(0);
      break;
    case CAPTK::ImageModalityType::IMAGE_TYPE_T2:
      //options_T2selected->setChecked(true);
      //options_T1selected->setChecked(false);
      options_T1Selector->setCurrentIndex(1);
      break;
    default:
      // no longer needed since this is handled from the app layer itself
      //ShowErrorMessage("Modalities other than T1 or T2 are not supported in WhiteStripe."); 
      break;
    }
    return;
  }

  //void checkSkullStripOption()
  //{
  //  if (skullStrippedImage)
  //  {
  //    options_skullStrippedImage->setChecked(true);
  //    options_axialSlicing->setChecked(false);
  //  }
  //  else
  //  {
  //    options_skullStrippedImage->setChecked(false);
  //    options_axialSlicing->setChecked(true);
  //  }
  //}

  //void checkT1SelectedOption()
  //{
  //  if (t1Image)
  //  {
  //    options_T1selected->setChecked(true);
  //    options_T2selected->setChecked(false);
  //  }
  //  else
  //  {
  //    options_T2selected->setChecked(true);
  //    options_T1selected->setChecked(false);
  //  }
  //}

  QString mInputPathName;


  public slots:
  void CancelButtonPressed();
  void ConfirmButtonPressed();

  void SelectOutputImage();

  //void T1ImageSelected();
  //void T2ImageSelected();
  void T1ComboBoxInvoked(int index);

  void SetSkullStrippedImage();
  void SetAxialSlicingNeeded();

signals:
  
  void RunWhiteStripe(double twsWidth, int sliceStartZ, int sliceStopZ, int tissuesMax, double smoothMax, double smoothDelta, int histSize,
    bool T1Image, const std::string outputFileName);
};


#endif
