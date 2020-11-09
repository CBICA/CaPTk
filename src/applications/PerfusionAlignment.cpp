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

#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkPlot.h>

void PerfusionAlignment::Createchart(vtkChartXY* chart, std::map<std::string, std::vector<float>*>* map)
{
	for (std::map<std::string, std::vector<float>*>::iterator itr = map->begin(); itr != map->end(); ++itr)
	{
		vtkNew< vtkTable> table;
		vtkNew<vtkFloatArray> a1, a2;
		a1->SetName("X");
		a2->SetName(itr->first.c_str());
		table->AddColumn(a1);
		table->AddColumn(a2);
		table->SetNumberOfRows(itr->second->size());

		for (vtkIdType i = 0; i < itr->second->size(); i++)
		{
			//fill table
			table->SetValue(i, 0, i);
			table->SetValue(i, 1, itr->second->at(i));
		}

		//add to chart
		vtkPlot *line = chart->AddPlot(vtkChart::POINTS);
		line->SetInputData(table, 0, 1);
		line->SetColor(0, 255.0, 0, 255);
		line->SetWidth(1.0);
	}
}