/**
\file  fProgressDialog.cpp

\brief Implementation of fProgressDialog class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "fProgressDialog.h"

fProgressDialog::fProgressDialog() 
{
  // don't do anything in this case since it is a default contructor
}

fProgressDialog::fProgressDialog(const fProgressDialog& origin)
{
  // don't do anything in this case since it is a default contructor
}

void fProgressDialog::/*fProgressDialog*/Initialize(const std::string &message, const bool show_progress = false)
{
  setupUi(this);
  textLabel->setText(message.c_str());
  if (show_progress) 
  {
    progressBar->show();
  }
  else 
  {
    progressBar->hide();
  }
  this->show();
}

void fProgressDialog::SetText(const std::string &message)
{
  textLabel->setText(message.c_str());
}

void fProgressDialog::AddToText(const std::string &message)
{
  textLabel->setText(QString("%1\n%2").arg(textLabel->text()).arg(message.c_str()));
}

void fProgressDialog::SetProgress(const unsigned int &current, const unsigned int &max)
{
  progressBar->setMaximum(max);
  progressBar->setValue(current);
}
