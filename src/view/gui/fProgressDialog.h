/**
\file  fProgressDialog.h

\brief Declaration of fProgressDialog class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _fProgressDialog_h_
#define _fProgressDialog_h_


//#include "CAPTk.h"
#include "ui_fProgressDialog.h"


/**
\class fProgressDialog

\brief This class controls the elements in the progress dialog
*/
class fProgressDialog : public QDialog, private Ui::fProgressDialog
{
  Q_OBJECT

public:
  fProgressDialog();
  fProgressDialog(const fProgressDialog& origin);
  //fProgressDialog(std::string message, bool show_progress = false);
  ~fProgressDialog() {}
  void Initialize(const std::string &message, const bool show_progress /*= false*/);
  void SetText(const std::string &message);
  void AddToText(const std::string &message);
  void SetProgress(const unsigned int &current, const unsigned int &max);
};


#endif
