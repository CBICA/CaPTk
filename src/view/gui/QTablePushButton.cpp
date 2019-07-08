/**
\file  QTablePushButton.cpp

\brief Implementation of QTablePushButton class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "QTablePushButton.h"

QTablePushButton::QTablePushButton() : QPushButton()
{
  m_item = NULL;
  m_index = 0;
  connect(this, SIGNAL(clicked()), this, SLOT(clicked()));
}

void QTablePushButton::clicked()
{
  emit clickedInto(m_item);
}

void QTablePushButton::setIndex(int index)
{
  m_index = index;
}

void QTablePushButton::setItem(QTableWidgetItem* item)
{
  m_item = item;
}
