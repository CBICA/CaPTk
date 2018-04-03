/**
\file  QTablePushButton.cpp

\brief Implementation of QTablePushButton class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html

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
  emit clickedInto(m_index);
}
