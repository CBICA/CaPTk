/**
\file ApplicationBase.h

This file holds the declaration of the class ApplicationBase.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/

#pragma once

#include <stdio.h>

//#include "fProgressDialog.h" // used to initialize progress dialog
//#include "CAPTk.h"

/**
\class ApplicationBase 

\brief This class is meant to be used after being inherited by an algorithm that runs in CAPTk.

Once an application class inherits this, they will need to implement the following methods:

1. SetLongRunning - This denotes whether the application inheriting is long runnig or not (approximate run time is >5 minutes)
2. Run - This basically tells the application to start the correlated computation(s).
*/

#ifndef APP_BASE_CAPTK_H
#define APP_BASE_CAPTK_H
#endif

#include <QObject>

class ApplicationBase : public QObject
{
	Q_OBJECT
private:
    //Q_DISABLE_COPY(ApplicationBase)
    
public:
  //! Default Constructor
  ApplicationBase() {};

  //! Default Destructor
  virtual ~ApplicationBase() {};

  void progressUpdate(const int val)
  { 
	  emit signalProgress(val);
  }
  void messageUpdate(const QString val)
  {
	  emit signalMessage(val);
  }
  virtual void Run()
  {

  }

  std::string m_LastError;
signals:
  void signalProgress(int);
  void signalMessage(QString);

public slots:

};
//#endif
