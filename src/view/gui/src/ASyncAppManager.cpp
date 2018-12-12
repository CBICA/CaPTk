/**
\file ASyncAppManager.cpp

This file holds the definition of the class ASyncAppManager.

https://www.cbica.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

\author Sarthak Pati

Copyright (c) 2015 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.cbica.upenn.edu/sbia/software/license.html
*/
#include "ASyncAppManager.h"

#include "qmessagebox.h"
#include "qthread.h"
#include "qfuture.h"
//
////template< class TAppType>
//ASyncAppManager::~ASyncAppManager()
//{
//
//}
//
////template< class TAppType>
//void ASyncAppManager::Run()
//{
//  if (m_application.m_longRunning)
//  {
//    m_message = "This application will take a while. Might we suggest you lock your computer and get yourself some nice tea/coffee?";
//  }
//
//  QMessageBox *msgBox = new QMessageBox(this);
//  msgBox->setText("Waiting for application");
//  msgBox->setInformativeText(QString(m_message));
//  msgBox->setWindowModality(Qt::WindowModal);
//  msgBox->show();
//
//  QFuture<void> future = QtConcurrent::run(m_application.Run());
//  future.waitForFinished();
//  msgBox->close();
//}