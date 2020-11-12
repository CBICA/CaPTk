/**
\file  PerfusionAlignment.cpp

\brief The source file containing the PerfusionAlignment class, used to calculate PSR, RCBV, and PH
Author: Saima Rathore
Library Dependecies: ITK 4.7+ <br>

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html

*/


#include "PerfusionAlignment.h"

#include "vtkTable.h"
#include "vtkFloatArray.h"
#include "vtkPlot.h"
#include "vtkSmartPointer.h"
#include "vtkChartXY.h"
#include "vtkJPEGWriter.h"
#include "vtkRenderWindow.h"
#include "vtkContextView.h"
#include "vtkWindowToImageFilter.h"
#include "vtkContextScene.h"
#include "vtkRenderer.h"
#include "vtkLookupTable.h"
#include "vtkPlotPoints.h"
#include "vtkAxis.h"

//#include <vtkChartXY.h>
//#include <vtkTable.h>
//#include <vtkFloatArray.h>
//#include <vtkPlot.h>

//void PerfusionAlignment::Createchart(vtkChartXY* chart, std::map<std::string, std::vector<float>*>* map)
//{
//	for (std::map<std::string, std::vector<float>*>::iterator itr = map->begin(); itr != map->end(); ++itr)
//	{
//		vtkNew< vtkTable> table;
//		vtkNew<vtkFloatArray> a1, a2;
//		a1->SetName("X");
//		a2->SetName(itr->first.c_str());
//		table->AddColumn(a1);
//		table->AddColumn(a2);
//		table->SetNumberOfRows(itr->second->size());
//
//		for (vtkIdType i = 0; i < itr->second->size(); i++)
//		{
//			//fill table
//			table->SetValue(i, 0, i);
//			table->SetValue(i, 1, itr->second->at(i));
//		}
//
//		//add to chart
//		vtkPlot *line = chart->AddPlot(vtkChart::POINTS);
//		line->SetInputData(table, 0, 1);
//		line->SetColor(0, 255.0, 0, 255);
//		line->SetWidth(1.0);
//	}
//}
void PerfusionAlignment::SaveChart(const std::string &fileToSave)
{
	vtkNew< vtkChartXY > chart;
	size_t c = 0;
	for (auto itr = m_curves.begin(); itr != m_curves.end(); ++itr, ++c)
	{
		vtkNew<vtkTable> table;
		vtkNew<vtkFloatArray> a1, a2;
		a1->SetName("X");
		a2->SetName(itr->first.c_str());
		table->AddColumn(a1);
		table->AddColumn(a2);
		table->SetNumberOfRows(itr->second.size());

		for (vtkIdType i = 0; i < itr->second.size(); i++)
		{
			//fill table
			table->SetValue(i, 0, i);
			table->SetValue(i, 1, itr->second[i]);
		}

		chart->SetShowLegend(true);
		chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Timepoints");
		chart->GetAxis(vtkAxis::LEFT)->SetTitle("Mean Intensity Inside ROI");

		//add to chart
		auto line = chart->AddPlot(vtkChart::LINE);
		line->SetInputData(table, 0, 1);
		line->SetWidth(1.0);
		switch (c)
		{
		case 0:
		{
			line->SetColor(255, 0, 0, 255);
			vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CROSS);
			break;
		}
		case 1:
		{
			line->SetColor(0, 255, 0, 255);
			vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
			break;
		}
		case 2:
		{
			line->SetColor(0, 0, 255, 255);
			vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::SQUARE);
			break;
		}
		case 3:
		{
			line->SetColor(0, 255, 255, 255);
			vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::DIAMOND);
			break;
		}
		default:
			break;
		}
	}

	vtkNew<vtkContextView> pView;
	auto temp = pView->GetScene();
	pView->GetScene()->AddItem(chart.GetPointer());
	pView->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
	vtkNew<vtkRenderWindow> renderWindow;
	renderWindow->AddRenderer(pView->GetRenderer());
	renderWindow->SetSize(600, 600);
	//renderWindow->SetDPI(600);
	renderWindow->OffScreenRenderingOn();
	renderWindow->Render();
	vtkNew<vtkWindowToImageFilter> windowToImageFilter;
	windowToImageFilter->SetInput(renderWindow.Get());

	vtkNew<vtkJPEGWriter> writer;
	writer->SetFileName(fileToSave.c_str());
	writer->SetInputConnection(windowToImageFilter->GetOutputPort());
	writer->Write();
}


std::string PerfusionAlignment::ReadMetaDataTags(std::string filepaths, std::string tag)
{
	DicomMetadataReader dcmMetaReader;
	dcmMetaReader.SetFilePath(filepaths);
	bool readstatus = dcmMetaReader.ReadMetaData();
	std::string label;
	std::string value;

	if (readstatus)
	{
		bool tagFound = dcmMetaReader.GetTagValue(tag, label, value);
		if (tagFound)
			return value;
		else
			return "";
	}
	else
	{
		return "";
	}
}

std::vector<double> PerfusionAlignment::GetInterpolatedCurve(std::vector<double> averagecurve, double timeinseconds, double totaltimeduration)
{
	return averagecurve;
}
void PerfusionAlignment::GetParametersFromTheCurve(std::vector<double> curve, double &base, double&drop, double &maxval, double &minval)
{
	std::vector<double> CER;
	//curve = vq1;
	//CER = curve * 100 / (mean(curve(4:8)) - min(curve));
	//MAX(i) = mean(curve(4:8));
	//MIN(i) = min(curve);
	//Base1(i) = mean(curve(4:8));

	//Base2(i) = mean(CER(4:8));
	//CER = CER - (mean(CER(4:8))) + 300;

	//Drop = find(CER == min(CER));
	//SP_Drop(i) = Drop;
	//CER = CER(Drop - 17:Drop + 36);
	//plot(CER);      hold on;
	//CERData(i, :) = CER;

	base = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;
	maxval = *max_element(curve.begin(), curve.end());
	minval = *min_element(curve.begin(), curve.end());
	drop = std::min_element(curve.begin(), curve.end()) - curve.begin();

	//for (int index = 0; index < curve.size(); index++)
	//  CER.push_back((curve[index] * 100) / (((curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5) - min));

	//maxval = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;
	//minval = min;
	//base = (curve[3] + curve[4] + curve[5] + curve[6] + curve[7]) / 5;

	//double average_new_cer = (CER[3] + CER[4] + CER[5] + CER[6] + CER[7]) / 5;
	//for (int index = 0; index < CER.size(); index++)
	//  CER[index] = CER[index] - average_new_cer + 300;

	//drop = std::min_element(CER.begin(), CER.end()) - CER.begin();
}