/**
\file ASyncAppManager.h

This file holds the declaration of the class ASyncAppManager.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html
*/

#pragma once

#include <stdio.h>
#include "qwidget.h"
#include "qmessagebox.h"
#include "qthread.h"
#include "qfuture.h"

//#include "CAPTk.h"
#include "ApplicationBase.h" // used for conherence between different applications and the sync manager

/**
\class ASyncAppManager

\brief Class is templated over the Application class

This class is used for running application in an asynchronous manner. The application should have inherited properties from
the class 'BaseClass', which gives ASyncAppManager access to a few pre-defined functions for easy synchronization. This class 
uses QThread for syncing between the GUI and processing class.
*/
template< class TAppType >
class ASyncAppManager : public ApplicationBase
{
public:
  //! The input can only be defined with respect to the Application Type (TAppType)
  explicit ASyncAppManager(TAppType &inputApplication, QWidget * parent = 0) : 
    m_application(inputApplication)
  {
    m_message = "This application should finish quickly. Sit tight!";
  }

  //! Default destructor
  virtual ~ASyncAppManager()
  {

  }

  inline void SetLongRunning(bool m_application.m_longRunning) :
    m_longRunning(m_application.m_longRunning)
  {
    // nothing to do here because the single parameter of interest has already been set
  }

  //! This method calls the TAppType::Run() method, where the main processing happens
  void Run()
  {
    if (m_application.m_longRunning)
    {
      m_message = "This application will take a while. Might we suggest you lock your computer and get yourself some nice tea/coffee?";
    }

    QMessageBox *msgBox = new QMessageBox(this);
    msgBox->setText("Waiting for Application");
    msgBox->setInformativeText(QString(m_message));
    msgBox->setWindowModality(Qt::WindowModal);
    //msgBox->setAttribute(Qt::WA_DeleteOnClose); //makes sure the msgbox is deleted automatically when closed
    msgBox->show();

    // application runs in a different thread while future makes sure that the messageBox stays up while it does
    QFuture<void> future = QtConcurrent::run(m_application.Run()); // most likely will go before msgBox->show()
    future.waitForFinished();
    msgBox->close();
  }

private:
  //! The application which will be called
  TAppType m_application;

  //! Message to display in the message box, defaults to short-running application
  std::string m_message;
};
