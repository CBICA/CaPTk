/**
\file  QTablePushButton.h

\brief Declaration of QTablePushButton class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _QTablePushButton_h_
#define _QTablePushButton_h_


#include <QPushButton>

class QTableWidgetItem;


class QTablePushButton : public QPushButton
{
  Q_OBJECT

public:
  QTablePushButton();

  void setIndex(int index);
  void setItem(QTableWidgetItem* item);

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
