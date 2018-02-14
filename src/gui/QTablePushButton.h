/**
\file  QTablePushButton.h

\brief Declaration of QTablePushButton class

http://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or http://www.med.upenn.edu/sbia/software/license.html

*/

#ifndef _QTablePushButton_h_
#define _QTablePushButton_h_


#include "CAPTk.h"


class QTablePushButton : public QPushButton
{
  Q_OBJECT

public:
  QTablePushButton();

  void setIndex(int index) {
    m_index = index;
  }
  void setItem(QTableWidgetItem* item) {
    m_item = item;
  }

  public slots:
  void clicked();

signals:
  void clickedInto(QTableWidgetItem* item);
  void clickedInto(int index);

private:
  QTableWidgetItem* m_item;
  int m_index;
};


#endif
