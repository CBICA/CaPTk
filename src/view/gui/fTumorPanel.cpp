///////////////////////////////////////////////////////////////////////////////////////
// fTumorPanel.cxx
//
// Copyright (c) 2018. All rights reserved.
// Section of Biomedical Image Analysis
// Center for Biomedical Image Computing and Analytics
// Department of Radiology
// Perelman School of Medicine
// University of Pennsylvania
//
// Contact details: software@cbica.upenn.edu
//
// License Agreement: https://www.med.upenn.edu/cbica/software-agreement.html
///////////////////////////////////////////////////////////////////////////////////////

#include "fTumorPanel.h"
//#include "CAPTk.h"
#include "CaPTkGUIUtils.h"
#include "fMainWindow.h"
#include "SlicerManager.h"
#include "Landmarks.h"

#include "qsignalmapper.h"

class tEventFilter : public QObject
{
public:
  fTumorPanel * mTObj;
  tEventFilter(fTumorPanel * obj) :QObject()
  {
    mTObj = obj;
  };
  ~tEventFilter(){};

  bool eventFilter(QObject* object, QEvent* event)
  {
    if (event->type() == QEvent::KeyPress)
    {
      QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);

      switch (keyEvent->key())
      {
      case Qt::Key_Return:
      {
        mTObj->HandleKeyPressingEventTTable();
        return true;
      }
      case Qt::Key_Down:
      {
        mTObj->HandleDownKeyEventTTable();
        return true;
      }
      case Qt::Key_Up:
      {
        mTObj->HandleUpKeyEventTTable();
        return true;
      }
      case Qt::Key_Delete:
      {
        mTObj->HandleDeleteKeyEventTTable();
        return true;
      }
      default:
      {
        cbica::Logging(loggerFile, "Undefined Key press.");
        //exit(EXIT_FAILURE); // probably not required
        return false;
      }
      }
    }
    else
    {
      return QObject::eventFilter(object, event);
    }
  }
};
class sEventFilter : public QObject
{
public:
  fTumorPanel * mTObj;
  sEventFilter(fTumorPanel * obj) :QObject()
  {
    mTObj = obj;
  };
  ~sEventFilter(){};

  bool eventFilter(QObject* object, QEvent* event)
  {
    if (event->type() == QEvent::KeyPress)
    {
      QKeyEvent *keyEvent = static_cast<QKeyEvent *>(event);

      switch (keyEvent->key())
      {
      case Qt::Key_Return:
      {
        mTObj->HandleKeyPressingEventSTable();
        return true;
      }
      case Qt::Key_Down:
      {
        mTObj->HandleDownKeyEventSTable();
        return true;
      }
      case Qt::Key_Up:
      {
        mTObj->HandleUpKeyEventTTable();
        return true;
      }
      case Qt::Key_Delete:
      {
        mTObj->HandleDeleteKeyEventTTable();
        return true;
      }
      default:
      {
        cbica::Logging(loggerFile, "Undefined Key press.");
        //exit(EXIT_FAILURE); // probably not required
        return false;
      }
      }
    }
    else
    {
      return QObject::eventFilter(object, event);
    }
  }
};
//static fMainWindow* getMainWindow()
//{
//	QWidgetList widgets = qApp->topLevelWidgets();
//	for (QWidgetList::iterator i = widgets.begin(); i != widgets.end(); ++i)
//	{
//		if ((*i)->objectName() == "MainWindow")
//		{
//			return qobject_cast<fMainWindow*>(*i);
//		}
//	}
//	return qobject_cast<fMainWindow*>(qApp->activeWindow());
//}
fTumorPanel::fTumorPanel(QWidget * parent) : QWidget(parent)
{
  initializationFileName = "initializedPoints_all.txt";
  setupUi(this);
  mTissueType = 0;
  mCurrentSPoints = NULL;
  mCurrentTPoints = NULL;
  sTableWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
  sTableWidget->setSelectionMode(QAbstractItemView::SingleSelection);

  sLoadButton->setEnabled(false);
  sSaveButton->setEnabled(false);
  sRemoveButton->setEnabled(false);
  sRemoveAllButton->setEnabled(false);
  tLoadButton->setEnabled(false);
  tSaveButton->setEnabled(false);
  tRemoveButton->setEnabled(false);
  tRemoveAllButton->setEnabled(false);

  connect(sLoadButton, SIGNAL(clicked()), this, SLOT(sLoad()));
  connect(sSaveButton, SIGNAL(clicked()), this, SLOT(sSave()));
  connect(sRemoveButton, SIGNAL(clicked()), this, SLOT(sRemoveSelectedPoints()));
  connect(sRemoveAllButton, SIGNAL(clicked()), this, SLOT(sRemoveAllPoints()));
  connect(tLoadButton, SIGNAL(clicked()), this, SLOT(tLoad()));
  connect(tSaveButton, SIGNAL(clicked()), this, SLOT(tSave()));
  connect(tRemoveButton, SIGNAL(clicked()), this, SLOT(tRemoveSelectedPoints()));
  connect(tRemoveAllButton, SIGNAL(clicked()), this, SLOT(tRemoveAllPoints()));
  connect(sTableWidget, SIGNAL(cellDoubleClicked(int, int)), this, SLOT(sTableDoubleClicked(int, int)));
  connect(tTableWidget, SIGNAL(cellDoubleClicked(int, int)), this, SLOT(tTableDoubleClicked(int, int)));
  connect(HelpButton, SIGNAL(clicked()), this, SLOT(helpClicked()));

  /// genericLabels
  QSignalMapper* signalMapper = new QSignalMapper(this); // needed to enable a signal to pass a value to a slot
  for (size_t i = 0; i < NumberOfTissueTypes; i++)
  {
    connect(m_vectorOfSeedPointLabels[i].radioButton, SIGNAL(toggled(bool)), signalMapper, SLOT(map()));
    signalMapper->setMapping(m_vectorOfSeedPointLabels[i].radioButton, i);
  }
  connect(signalMapper, SIGNAL(mapped(int)), this, SLOT(SetTissueType(int)));
  /// genericLabels

  connect(m_typeRadBtnTumor, SIGNAL(toggled(bool)), this, SLOT(SetSeedType()));
  connect(m_typeRadBtnAllTissues, SIGNAL(toggled(bool)), this, SLOT(SetTissueType()));
  connect(m_typeRadBtnGlister, SIGNAL(toggled(bool)), this, SLOT(SetGLISTRTissueType()));
  connect(m_typeRadBtnPorterPre, SIGNAL(toggled(bool)), this, SLOT(SetPORTRPRETissueType()));
  connect(m_typeRadBtnPorterPost, SIGNAL(toggled(bool)), this, SLOT(SetPORTRPOSTTissueType()));

  tTableWidget->installEventFilter(new tEventFilter(this));
  sTableWidget->installEventFilter(new sEventFilter(this));

  sTableWidget->resizeColumnsToContents();
  sTableWidget->resizeRowsToContents();
  tTableWidget->resizeColumnsToContents();
  tTableWidget->resizeRowsToContents();

  emit SetTissueCounter(mType);
  m_vectorOfSeedPointLabels[CSF].radioButton->setChecked(true);
  SetTissueType(CSF);
  m_typeRadBtnTumor->setChecked(false);
  m_typeRadBtnAllTissues->setChecked(true);
  UpdateCurrentPoints();

}

void fTumorPanel::Clear()
{
  sTableWidget->clearContents();
  tTableWidget->clearContents();
  sLoadButton->setEnabled(false);
  sSaveButton->setEnabled(false);
  sRemoveButton->setEnabled(false);
  sRemoveAllButton->setEnabled(false);
  tLoadButton->setEnabled(false);
  tSaveButton->setEnabled(false);
  tRemoveButton->setEnabled(false);
  tRemoveAllButton->setEnabled(false);
}

void fTumorPanel::sLoad()
{
  m_typeRadBtnTumor->setChecked(true);
  QString file;
  QFileDialog dialog(this, tr("Load Seed Points"), mCurrentPath.c_str(), tr("Seed Points (*.txt)"));
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptOpen);
  dialog.selectFile("init_seed.txt");

  if (dialog.exec() == QDialog::Accepted)
  {
    file = dialog.selectedFiles()[0];

    sLoad(file);
  }
}

void fTumorPanel::sLoad(QString file)
{
  QFile inputFile(file);
  double scan_seeds_info[LANDMARK_MAX_TUMOR_POINTS][4];
  int scan_seeds_num;

  if (inputFile.open(QIODevice::ReadOnly))
  {
    QTextStream in(&inputFile);
    scan_seeds_num = (int)atoi(in.readLine().toStdString().c_str());
    if (scan_seeds_num > LANDMARK_MAX_TUMOR_POINTS)
    {
      ShowErrorMessage("Number of tumor points in the file exceed the max number that can be processed; returning");
      return;
    }
    mCurrentSPoints->Clear();
    for (int i = 0; i < scan_seeds_num; i++)
    {
      std::vector<std::string > values = cbica::stringSplit(in.readLine().toStdString(), " ");
      if (values.size() != 4)
      {
        ShowErrorMessage("The information for tumor point # " + std::to_string(i) + " is incorrect, please check");
        return;
      }
      scan_seeds_info[i][0] = static_cast<double>(std::atof(values[0].c_str()));
      scan_seeds_info[i][1] = static_cast<double>(std::atof(values[1].c_str()));
      scan_seeds_info[i][2] = static_cast<double>(std::atof(values[2].c_str()));
      scan_seeds_info[i][3] = static_cast<double>(std::atof(values[3].c_str()));

#ifdef USE_LANDMARKS_LPS_COORD
      mCurrentSPoints->AddLandmark(scan_seeds_info[i][0], scan_seeds_info[i][1], scan_seeds_info[i][2], scan_seeds_info[i][3] / 2.0, 0, i);
#endif

#ifdef USE_LANDMARKS_RAS_COORD
      mCurrentSPoints->AddLandmark(-scan_seeds_info[i][0], -scan_seeds_info[i][1], scan_seeds_info[i][2], scan_seeds_info[i][3] / 2.0, 0, i);
#endif
    }
  }
  SetCurrentSPoints(mCurrentSPoints);
  emit UpdateRenderWindows();
}

void fTumorPanel::sSave()
{
  m_typeRadBtnTumor->setChecked(true);
  QString file;
  QFileDialog dialog(this, tr("Save Seed Points"), mCurrentPath.c_str(), tr("Seed Points (*.txt)"));
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.selectFile("init_seed.txt");

  int ret = dialog.exec();
  if (ret == QDialog::Accepted)
  {
    file = dialog.selectedFiles()[0];
    //*/
    FILE* fp;

    // check write access
    //if (((_access(file.toStdString().c_str(), 2)) == -1) || ((_access(file.toStdString().c_str(), 6)) == -1))
    //{
    //  ShowErrorMessage("You don't have write access in selected location. Please check.");
    //}

    //QFile fileToWrite(file);
    //if (!fileToWrite.open(QIODevice::WriteOnly | QIODevice::Text))
    //{
    //  printf("Cannot open %s\n", file.toStdString().c_str());
    //  return;
    //}
    //QTextStream out(&fileToWrite);
    //out << mCurrentSPoints->GetNumberOfPoints() << "\n";

#ifdef _WIN32
    fopen_s(&fp, file.toStdString().c_str(), "w");
#else
    fp = fopen(file.toStdString().c_str(), "w");
#endif

    if (fp == NULL)
    {
      printf("Cannot open %s\n", file.toStdString().c_str());
      return;
    }
    fprintf(fp, "%d\n", mCurrentSPoints->GetNumberOfPoints());
    for (unsigned int i = 0; i < mCurrentSPoints->GetNumberOfPoints(); i++)
    {
#ifdef USE_LANDMARKS_LPS_COORD
      fprintf(fp, "%f %f %f %f\n", mCurrentSPoints->mLandmarks[i].coordinates[0], mCurrentSPoints->mLandmarks[i].coordinates[1], mCurrentSPoints->mLandmarks[i].coordinates[2], mCurrentSPoints->mLandmarks[i].radius * 2.0);
      //out << mCurrentSPoints->mLandmarks[i].coordinates[0] << " " << mCurrentSPoints->mLandmarks[i].coordinates[1] << " " << mCurrentSPoints->mLandmarks[i].coordinates[2] << " " << mCurrentSPoints->mLandmarks[i].radius * 2.0 << "\n";
#endif
#ifdef USE_LANDMARKS_RAS_COORD
      fprintf(fp, "%f %f %f %f\n", -mCurrentSPoints->mLandmarks[i].coordinates[0], -mCurrentSPoints->mLandmarks[i].coordinates[1], mCurrentSPoints->mLandmarks[i].coordinates[2], mCurrentSPoints->mLandmarks[i].radius * 2.0);
      //out << -mCurrentSPoints->mLandmarks[i].coordinates[0] << " " << -mCurrentSPoints->mLandmarks[i].coordinates[1] << " " << mCurrentSPoints->mLandmarks[i].coordinates[2] << " " << mCurrentSPoints->mLandmarks[i].radius * 2.0 << "\n";
#endif
    }
    //fileToWrite.close();
    fclose(fp);
  }
}

void fTumorPanel::tLoad()
{
  m_typeRadBtnAllTissues->setChecked(true);
  QString file;
  QFileDialog dialog(this, tr("Load Tissue Points"), mCurrentPath.c_str(), tr("Tissue Points (*.txt)"));
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptOpen);
  dialog.selectFile("init_point.txt");

  if (dialog.exec() == QDialog::Accepted)
  {
    file = dialog.selectedFiles()[0];

    tLoad(file); // call actual reader here
  }
}
//
void fTumorPanel::tLoad(QString file)
{
  QFile inputFile(file);

  unsigned int i;
  int j;

  // check read access
  /*if (((_access(file.toStdString().c_str(), 4)) == -1) || ((_access(file.toStdString().c_str(), 6)) == -1))
  {
  ShowErrorMessage("You don't have read access in selected location. Please check.");
  return;
  }*/

  mCurrentTPoints->Clear();


  if (inputFile.open(QIODevice::ReadOnly))
  {
    QTextStream in(&inputFile);
    while (!in.atEnd())
    {
      std::string currentLine = in.readLine().toStdString();
      if (!currentLine.empty())
      {
        for (size_t k = 0; k < NumberOfTissueTypes; k++)
        {
          if (currentLine == labels_laconic[k])
          {
            // tissue point matched!
            k = NumberOfTissueTypes;
          }
          if ((currentLine != labels_laconic[k]) && (k == NumberOfTissueTypes - 1))
          {
            ShowErrorMessage("Input tissue point file has an incorrect tissue type, cannot load.");
            return;
          }
        }
        std::string id = currentLine;

        std::vector< std::string > values = cbica::stringSplit(in.readLine().toStdString(), " ");

        if (values.size() != 3)
        {
          ShowErrorMessage("Input tissue point file has an incorrect number of coordinates, cannot load.");
          return;
        }

        double fx, fy, fz;
        fx = static_cast<double>(std::atof(values[0].c_str()));
        fy = static_cast<double>(std::atof(values[1].c_str()));
        fz = static_cast<double>(std::atof(values[2].c_str()));

        for (j = 0; j < NumberOfTissueTypes; j++)
        {
          if (id == labels_laconic[j])
          {
            if (j == 0)
            {
              // skip BG
            }
            else
            {
#ifdef USE_LANDMARKS_LPS_COORD
              mCurrentTPoints->AddLandmark(fx, fy, fz, 0, 0, j);
#endif
#ifdef USE_LANDMARKS_RAS_COORD
              mCurrentTPoints->AddLandmark(-fx, -fy, fz, 0, 0, j);
#endif
            }
            break;
          }
        }
        if (j == NumberOfTissueTypes)
        {
          printf("class id is wrong\n");
          mCurrentTPoints->Clear();
          return;
        }
      }
      inputFile.close();
    }
  }
  for (i = 0; i < mCurrentTPoints->GetNumberOfPoints(); i++)
  {
    if (mCurrentTPoints->mLandmarks[i].bValid)
    {
      printf("point %d: %s %f %f %f\n", i, labels_laconic[mCurrentTPoints->mLandmarks[i].id], mCurrentTPoints->mLandmarks[i].coordinates[0], mCurrentTPoints->mLandmarks[i].coordinates[1], mCurrentTPoints->mLandmarks[i].coordinates[2]);
    }
  }
  SetCurrentTPoints(mCurrentTPoints);
  emit UpdateRenderWindows();

}
void fTumorPanel::tSave(QString file)
{
  FILE* fp;
#ifdef _WIN32
  fopen_s(&fp, file.toStdString().c_str(), "w");
#else
  fp = fopen(file.toStdString().c_str(), "w");
#endif

  if (fp == NULL)
  {
    printf("Cannot open %s\n", file.toStdString().c_str());
    return;
  }
  for (int i = 0; i < NumberOfTissueTypes; i++)
  {
    for (unsigned int j = 0; j < mCurrentTPoints->GetNumberOfPoints(); j++)
    {
      if (mCurrentTPoints->mLandmarks[j].id == i && mCurrentTPoints->mLandmarks[j].bValid)
      {
#ifdef USE_LANDMARKS_LPS_COORD
        fprintf(fp, "%s\n%f %f %f\n", labels_laconic[i], mCurrentTPoints->mLandmarks[j].coordinates[0], mCurrentTPoints->mLandmarks[j].coordinates[1], mCurrentTPoints->mLandmarks[j].coordinates[2]);
#endif
#ifdef USE_LANDMARKS_RAS_COORD
        fprintf(fp, "%s\n%f %f %f\n", labels_laconic[i], -mCurrentTPoints->mLandmarks[j].coordinates[0], -mCurrentTPoints->mLandmarks[j].coordinates[1], mCurrentTPoints->mLandmarks[j].coordinates[2]);
#endif
      }
    }
  }
  fclose(fp);
}
void fTumorPanel::tSave()
{
  //m_typeRadBtnAllTissues->setChecked(true);
  // check write access
  //if (((_access(file.toStdString().c_str(), 2)) == -1) || ((_access(file.toStdString().c_str(), 6)) == -1))
  //{
  //  ShowErrorMessage("You don't have write access in selected location. Please check.");
  //}

  // check for saveType_* for three and compare different tissue types
  std::vector< unsigned int > tissuesToCheck;
  /// these "required" tissue types have been decided after discussion with Spyros
  switch (m_saveType)
  {
  case SaveType::GLISTR:
  {
    tissuesToCheck.push_back(CSF);
    tissuesToCheck.push_back(GM);
    tissuesToCheck.push_back(WM);
    tissuesToCheck.push_back(VS);
    //tissuesToCheck.push_back(ED);
    //tissuesToCheck.push_back(NCR);
    tissuesToCheck.push_back(TU);
    tissuesToCheck.push_back(NE); // one of either this or above should be present
    tissuesToCheck.push_back(CB);
    break;
  }
  case SaveType::PORTRPost:
  {
    tissuesToCheck.push_back(CSF);
    tissuesToCheck.push_back(VT); // one of either this or above should be present
    tissuesToCheck.push_back(GM);
    tissuesToCheck.push_back(WM);
    tissuesToCheck.push_back(VS);
    //tissuesToCheck.push_back(ED);
    tissuesToCheck.push_back(CAN);
    tissuesToCheck.push_back(CAE); // one of either this or above should be present
    //tissuesToCheck.push_back(RTN);
    //tissuesToCheck.push_back(RTE);
    break;
  }
  case SaveType::PORTRPre:
  {
    tissuesToCheck.push_back(CSF);
    tissuesToCheck.push_back(VT); // one of either this or above should be present
    tissuesToCheck.push_back(GM);
    tissuesToCheck.push_back(WM);
    tissuesToCheck.push_back(VS);
    //tissuesToCheck.push_back(ED);
    tissuesToCheck.push_back(NCR);
    tissuesToCheck.push_back(TU); // one of either this or above should be present
    break;
  }

  default:
    break;
  }

  if (!tissuesToCheck.empty())
  {
    for (size_t i = 0; i < mCurrentTPoints->mLandmarks.size(); i++)
    {
      if (!tissuesToCheck.empty() && mCurrentTPoints->mLandmarks[i].bValid)
      {
        std::vector< unsigned int >::iterator iterator = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), mCurrentTPoints->mLandmarks[i].id);
        if (iterator != tissuesToCheck.end())
        {
          tissuesToCheck.erase(iterator);
        }
      }
      else
      {
        break;
      }
    }
  }

  bool glistrOptionalsFound = false, portrPreOptionalsFound = false, portrPostOptionalsFound = false;

  if (tissuesToCheck.size() > 0)
  {
    // check for optional tissue types for the specific applications
    std::string tissues_later = "";
    std::vector< unsigned int >::iterator it_1, it_2;
    if (m_saveType == SaveType::GLISTR)
    {
      it_1 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), TU);
      it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), NE);
      if ((it_1 != tissuesToCheck.end()) && (it_2 != tissuesToCheck.end()))
      {
        tissues_later += "  Either TU (Enhancing Tumor) or NE (Non-enhancing Tumor) is required.";
        // ensure the checked tissue is removed
        if (!mCurrentTPoints->mLandmarks.empty())
        {
          if (it_1 != tissuesToCheck.end())
          {
            tissuesToCheck.erase(it_1);
            it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), NE);
          }

          if (it_2 != tissuesToCheck.end())
            tissuesToCheck.erase(it_2);
        }
      }
      else
      {
        if (!mCurrentTPoints->mLandmarks.empty())
          glistrOptionalsFound = true;
      }
    }
    else if (m_saveType > 0)
    {
      it_1 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), CSF);
      it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), VT);
      if ((it_1 != tissuesToCheck.end()) && (it_2 != tissuesToCheck.end()))
      {
        tissues_later += "  Either CSF (Cerebrospinal Fluid) or VT (Ventricular CSF) is required.\n";
        // ensure the checked tissue is removed
        if (!mCurrentTPoints->mLandmarks.empty())
        {
          if (it_1 != tissuesToCheck.end())
          {
            tissuesToCheck.erase(it_1);
            it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), VT);
          }

          if (it_2 != tissuesToCheck.end())
            tissuesToCheck.erase(it_2);
        }
      }
      else
      {
        if (!mCurrentTPoints->mLandmarks.empty())
          portrPostOptionalsFound = true;
      }

      it_1 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), CAN);
      it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), CAE);
      if ((it_1 != tissuesToCheck.end()) && (it_2 != tissuesToCheck.end()))
      {
        tissues_later += "  Either CAN (Cavity Non Enhancing) or CAE (Cavity Enhancing) is required.\n";
        // ensure the checked tissue is removed
        if (!mCurrentTPoints->mLandmarks.empty())
        {
          if (it_1 != tissuesToCheck.end())
          {
            tissuesToCheck.erase(it_1);
            it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), CAE);
          }

          if (it_2 != tissuesToCheck.end())
            tissuesToCheck.erase(it_2);
        }
      }
      else
      {
        if (!mCurrentTPoints->mLandmarks.empty())
          portrPostOptionalsFound = true;
      }

      if (m_saveType == SaveType::PORTRPre)
      {
        it_1 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), NCR);
        it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), TU);
        if ((it_1 != tissuesToCheck.end()) && (it_2 != tissuesToCheck.end()))
        {
          tissues_later += "  Either NCR (Necrosis) or TU (Enhancing Tumor) is required.";
          // ensure the checked tissue is removed
          if (!mCurrentTPoints->mLandmarks.empty())
          {
            if (it_1 != tissuesToCheck.end())
            {
              tissuesToCheck.erase(it_1);
              it_2 = std::find(tissuesToCheck.begin(), tissuesToCheck.end(), TU);
            }

            if (it_2 != tissuesToCheck.end())
              tissuesToCheck.erase(it_2);
          }
        }
        else
        {
          if (!mCurrentTPoints->mLandmarks.empty())
            portrPreOptionalsFound = true;
        }
      }

      if (m_saveType == SaveType::PORTRPost)
      {
        portrPreOptionalsFound = true;
      }
    }

    if (!tissuesToCheck.empty() || !tissues_later.empty())
    {
      std::string tissues = "";
      for (size_t i = 0; i < tissuesToCheck.size(); i++)
      {
        switch (tissuesToCheck[i])
        {
        case CSF:
          if (!portrPostOptionalsFound || !portrPreOptionalsFound)
            tissues += "  CSF (Cerebrospinal Fluid)\n";
          break;
        case GM:
          tissues += "  GM  (Gray Matter)\n";
          break;
        case WM:
          tissues += "  WM  (White Matter)\n";
          break;
        case VS:
          tissues += "  VS  (Vessels)\n";
          break;
        case ED:
          tissues += "  ED  (Edema)\n";
          break;
        case NCR:
          if (!portrPreOptionalsFound)
            tissues += "  NCR (Necrosis)\n";
          break;
        case TU:
          if (!glistrOptionalsFound && !portrPreOptionalsFound)
            tissues += "  TU  (Enhancing Tumor)\n";
          break;
        case NE:
          if (!glistrOptionalsFound)
            tissues += "  NE  (Non-enhancing Tumor)\n";
          break;
        case CB:
          tissues += "  CB  (Cerebellum)\n";
          break;
        case VT:
          if (!portrPostOptionalsFound && !portrPreOptionalsFound)
            tissues += "  VT  (Ventricular CSF)\n";
          break;
        case CAN:
          if (!portrPostOptionalsFound)
            tissues += "  CAN (Non-enhancing Cavity)\n";
          break;
        case CAE:
          if (!portrPostOptionalsFound)
            tissues += "  CAE (Enhancing Cavity)\n";
          break;
        case RTN:
          tissues += "  RTN (Tumor Recurrence)\n";
          break;
        case RTE:
          tissues += "  RTE (Enhanced Tumor Recurrence)\n";
          break;
        default:
          break;
        }
      }

      if (!(tissues.empty() && tissues_later.empty()))
      {
        // this needs to show the text in a better aligned way
        QMessageBox msgBox;
        msgBox.setText("The following required tissue type(s) are either not initialized or aren't valid:");
        if (!tissues.empty())
        {
          tissues = tissues.substr(0, tissues.size() - 1); // remove trailing '\n'
          msgBox.setInformativeText((tissues + "\n\n" + tissues_later).c_str());
        }
        else
        {
          msgBox.setInformativeText(tissues_later.c_str());
        }
        msgBox.setStandardButtons(QMessageBox::Ok);
        msgBox.setDefaultButton(QMessageBox::Ok);
        msgBox.exec();

        return;
      }
    }
  }

  QFileDialog dialog(this, tr("Save Tissue Points"), mCurrentPath.c_str(), tr("Tissue Points (*.txt)"));
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.selectFile(initializationFileName.c_str());
  int ret = dialog.exec();
  QString file;

  if (ret == QDialog::Accepted)
  {
    file = dialog.selectedFiles()[0];
    //*/
    FILE* fp;
#ifdef _WIN32
    fopen_s(&fp, file.toStdString().c_str(), "w");
#else
    fp = fopen(file.toStdString().c_str(), "w");
#endif

    if (fp == NULL)
    {
      printf("Cannot open %s\n", file.toStdString().c_str());
      return;
    }
    for (int i = 0; i < NumberOfTissueTypes; i++)
    {
      for (unsigned int j = 0; j < mCurrentTPoints->GetNumberOfPoints(); j++)
      {
        if (mCurrentTPoints->mLandmarks[j].id == i && mCurrentTPoints->mLandmarks[j].bValid)
        {
#ifdef USE_LANDMARKS_LPS_COORD
          fprintf(fp, "%s\n%f %f %f\n", labels_laconic[i], mCurrentTPoints->mLandmarks[j].coordinates[0], mCurrentTPoints->mLandmarks[j].coordinates[1], mCurrentTPoints->mLandmarks[j].coordinates[2]);
#endif
#ifdef USE_LANDMARKS_RAS_COORD
          fprintf(fp, "%s\n%f %f %f\n", labels_laconic[i], -mCurrentTPoints->mLandmarks[j].coordinates[0], -mCurrentTPoints->mLandmarks[j].coordinates[1], mCurrentTPoints->mLandmarks[j].coordinates[2]);
#endif
        }
      }
    }
    fclose(fp);
  }
}

void fTumorPanel::sRemoveSelectedPoints()
{
  QList<QTableWidgetItem*> items = sTableWidget->selectedItems();
  if (items.empty())
    return;

  int rowIndex = items[0]->row();
  std::vector<int> rowIndices;
  rowIndices.push_back(rowIndex);

  for (int i = 0; i < items.size(); i++)
  {
    if ((std::find(rowIndices.begin(), rowIndices.end(), items[i]->row()) != rowIndices.end()) == false)
      rowIndices.push_back(items[i]->row());
  }
  for (unsigned int j = 0; j < rowIndices.size(); j++)
    sRemoveCurrentIndexPoint(rowIndex);
}

void fTumorPanel::tRemoveCurrentIndexPoint(unsigned int rowIndex)
{
  mCurrentTPoints->RemoveLandmark(rowIndex);
  tTableWidget->removeRow(rowIndex);
  emit UpdateRenderWindows();

  if (rowIndex > mCurrentTPoints->GetNumberOfPoints())
    return;

  tTableWidget->selectRow(rowIndex);
  QTableWidgetItem * item2 = tTableWidget->item(rowIndex, 0);
  tTableWidget->setItemSelected(item2, true);
  UpdateCurrentPoints();
}

void fTumorPanel::sRemoveCurrentIndexPoint(unsigned int rowIndex)
{
  mCurrentSPoints->RemoveLandmark(rowIndex);
  sTableWidget->removeRow(rowIndex);

  emit UpdateRenderWindows();

  if (rowIndex>mCurrentSPoints->GetNumberOfPoints())
    return;

  sTableWidget->selectRow(rowIndex);
  QTableWidgetItem * item2 = sTableWidget->item(rowIndex, 0);
  sTableWidget->setItemSelected(item2, true);
}

void fTumorPanel::tRemoveSelectedPoints()
{
  QList<QTableWidgetItem*> items = tTableWidget->selectedItems();
  if (items.empty())
    return;

  int rowIndex = items[0]->row();
  std::vector<int> rowIndices;
  rowIndices.push_back(rowIndex);

  for (int i = 0; i < items.size(); i++)
  {
    if ((std::find(rowIndices.begin(), rowIndices.end(), items[i]->row()) != rowIndices.end()) == false)
      rowIndices.push_back(items[i]->row());
  }
  for (unsigned int j = 0; j < rowIndices.size(); j++)
    tRemoveCurrentIndexPoint(rowIndex);

  UpdateCurrentPoints();
}

void fTumorPanel::sAddPoint()
{
  sAddPoint(mCurrentSPoints->GetNumberOfPoints() - 1, false);
}

void fTumorPanel::sAddPoint(int landmarksIndex, bool update)
{
  int rowIndex = landmarksIndex;// mCurrentSPoints->GetNumberOfPoints() - 1; //sTableWidget->rowCount() //landmarksIndex;
  float x, y, z;
  // double val;
  bool bValid;
  float r;

  //sTableWidget->setRowCount(rowIndex+1);
  //sTableWidget->setRowHeight(rowIndex, 20);

  if (update == false)
    sTableWidget->setRowCount(rowIndex + 1);


#ifdef USE_LANDMARKS_LPS_COORD
  x = mCurrentSPoints->mLandmarks[landmarksIndex].coordinates[0];
  y = mCurrentSPoints->mLandmarks[landmarksIndex].coordinates[1];
  z = mCurrentSPoints->mLandmarks[landmarksIndex].coordinates[2];
#endif
#ifdef USE_LANDMARKS_RAS_COORD
  x = -mCurrentSPoints->mLandmarks[landmarksIndex].coordinates[0];
  y = -mCurrentSPoints->mLandmarks[landmarksIndex].coordinates[1];
  z = mCurrentSPoints->mLandmarks[landmarksIndex].coordinates[2];
#endif
  // val = mCurrentSPoints->mLandmarks[landmarksIndex].pixel_value;
  bValid = mCurrentSPoints->mLandmarks[landmarksIndex].bValid;
  r = mCurrentSPoints->mLandmarks[landmarksIndex].radius;

  if (bValid)
  {
    QTableWidgetItem* rowCounterItem = new QTableWidgetItem(QString::number(rowIndex));
    rowCounterItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    sTableWidget->setVerticalHeaderItem(rowIndex, rowCounterItem);

    QTableWidgetItem* xItem = new QTableWidgetItem(QString::number(x, 'f', 3));
    //xItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    xItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    sTableWidget->setItem(rowIndex, 0, xItem);

    QTableWidgetItem* yItem = new QTableWidgetItem(QString::number(y, 'f', 3));
    //yItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    yItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    sTableWidget->setItem(rowIndex, 1, yItem);

    QTableWidgetItem* zItem = new QTableWidgetItem(QString::number(z, 'f', 3));
    //zItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    zItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    sTableWidget->setItem(rowIndex, 2, zItem);

    QTableWidgetItem* rItem = new QTableWidgetItem(QString::number(r, 'f', 3));
    //rItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    rItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    sTableWidget->setItem(rowIndex, 3, rItem);
  }
}

void fTumorPanel::tAddPoint(int rowIndex)
{
  float x, y, z;
  // double val;
  bool bValid;

  int idx = mCurrentTPoints->GetNumberOfPoints() - 1;
  tTableWidget->setRowCount(mCurrentTPoints->GetNumberOfPoints());

#ifdef USE_LANDMARKS_LPS_COORD
  x = mCurrentTPoints->mLandmarks[idx].coordinates[0];
  y = mCurrentTPoints->mLandmarks[idx].coordinates[1];
  z = mCurrentTPoints->mLandmarks[idx].coordinates[2];
#endif
#ifdef USE_LANDMARKS_RAS_COORD
  x = -mCurrentTPoints->mLandmarks[idx].coordinates[0];
  y = -mCurrentTPoints->mLandmarks[idx].coordinates[1];
  z = mCurrentTPoints->mLandmarks[idx].coordinates[2];
#endif
  // val = mCurrentTPoints->mLandmarks[idx].pixel_value;
  bValid = mCurrentTPoints->mLandmarks[idx].bValid;

  std::string tissue_label = labels_laconic[mCurrentTPoints->mLandmarks[idx].id];

  if (bValid)
  {

    QTableWidgetItem* labelItem = new QTableWidgetItem(QString::fromStdString(tissue_label));
    //xItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    labelItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    tTableWidget->setItem(idx, 0, labelItem);

    QTableWidgetItem* xItem = new QTableWidgetItem(QString::number(x, 'f', 3));
    //xItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    xItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    tTableWidget->setItem(idx, 1, xItem);

    QTableWidgetItem* yItem = new QTableWidgetItem(QString::number(y, 'f', 3));
    //yItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    yItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    tTableWidget->setItem(idx, 2, yItem);

    QTableWidgetItem* zItem = new QTableWidgetItem(QString::number(z, 'f', 3));
    //zItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    zItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
    tTableWidget->setItem(idx, 3, zItem);

    emit UpdateRenderWindows();
    tTableWidget->scrollToItem(labelItem, QAbstractItemView::PositionAtTop);
  }
  UpdateCurrentPoints();
}

void fTumorPanel::SetCurrentSPoints(Landmarks* lm)
{
  sLoadButton->setEnabled(true);
  sSaveButton->setEnabled(true);
  sRemoveButton->setEnabled(true);
  sRemoveAllButton->setEnabled(true);

  mCurrentSPoints = lm;
  sTableWidget->clearContents();
  for (unsigned int i = 0; i < mCurrentSPoints->GetNumberOfPoints(); i++)
  {
    sAddPoint(i, false);
  }
  sTableWidget->resizeColumnsToContents();
}
void fTumorPanel::SetCurrentTPoints(Landmarks* lm)
{
  tLoadButton->setEnabled(true);
  tSaveButton->setEnabled(true);
  tRemoveButton->setEnabled(true);
  tRemoveAllButton->setEnabled(true);
  mCurrentTPoints = lm;

  tTableWidget->clearContents();


  float x, y, z;
  // double val;
  bool bValid;
  tTableWidget->setRowCount(mCurrentTPoints->GetNumberOfPoints());
  for (unsigned int i = 0; i < mCurrentTPoints->GetNumberOfPoints(); i++)
  {
#ifdef USE_LANDMARKS_LPS_COORD
    x = mCurrentTPoints->mLandmarks[i].coordinates[0];
    y = mCurrentTPoints->mLandmarks[i].coordinates[1];
    z = mCurrentTPoints->mLandmarks[i].coordinates[2];
#endif
#ifdef USE_LANDMARKS_RAS_COORD
    x = -mCurrentTPoints->mLandmarks[i].coordinates[0];
    y = -mCurrentTPoints->mLandmarks[i].coordinates[1];
    z = mCurrentTPoints->mLandmarks[i].coordinates[2];
#endif
    // val = mCurrentTPoints->mLandmarks[i].pixel_value;
    bValid = mCurrentTPoints->mLandmarks[i].bValid;

    std::string tissue_label = labels_laconic[mCurrentTPoints->mLandmarks[i].id];


    if (bValid)
    {

      QTableWidgetItem* labelItem = new QTableWidgetItem(QString::fromStdString(tissue_label));
      //xItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      labelItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
      tTableWidget->setItem(i, 0, labelItem);

      QTableWidgetItem* xItem = new QTableWidgetItem(QString::number(x, 'f', 3));
      //xItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      xItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
      tTableWidget->setItem(i, 1, xItem);

      QTableWidgetItem* yItem = new QTableWidgetItem(QString::number(y, 'f', 3));
      //yItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      yItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
      tTableWidget->setItem(i, 2, yItem);

      QTableWidgetItem* zItem = new QTableWidgetItem(QString::number(z, 'f', 3));
      //zItem->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      zItem->setTextAlignment(Qt::AlignCenter | Qt::AlignVCenter);
      tTableWidget->setItem(i, 3, zItem);
    }
  }
  emit UpdateRenderWindows();
  //	
  //	for (unsigned int i = 0; i < tTableWidget->rowCount(); i++) 
  //{
  //	tAddPoint(i);
  //}
  tTableWidget->resizeColumnsToContents();
  UpdateCurrentPoints();
}

void fTumorPanel::sTableDoubleClicked(int row, int col)
{
  if (mCurrentSPoints == NULL)
    return;
  if (mCurrentSPoints->GetNumberOfPoints() == 0)
    return;

  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TUMOR_POINTS, row, col);

  if (row >= 0 && row < (int)mCurrentSPoints->GetNumberOfPoints())
  {
    if (mCurrentSPoints->mLandmarks[row].bValid)
    {
      double x, y, z;
      x = mCurrentSPoints->mLandmarks[row].coordinates[0];
      y = mCurrentSPoints->mLandmarks[row].coordinates[1];
      z = mCurrentSPoints->mLandmarks[row].coordinates[2];
      emit MoveSlicerCursor(x, y, z);
    }
  }
}
void fTumorPanel::tTableDoubleClicked(int row, int col)
{
  if (mCurrentTPoints == NULL)
    return;
  if (mCurrentTPoints->GetNumberOfPoints() == 0)
    return;

  mTissueType = mCurrentTPoints->mLandmarks[row].id;
  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TISSUE_POINTS, mTissueType, col);
  emit SetTissueTableIndex(row);

  switch (mTissueType)
  {
  case CSF:
    m_vectorOfSeedPointLabels[CSF].radioButton->setChecked(true);
    break;
  case GM:
    m_vectorOfSeedPointLabels[GM].radioButton->setChecked(true);
    break;
  case WM:
    m_vectorOfSeedPointLabels[WM].radioButton->setChecked(true);
    break;
  case VS:
    m_vectorOfSeedPointLabels[VS].radioButton->setChecked(true);
    break;
  case ED:
    m_vectorOfSeedPointLabels[ED].radioButton->setChecked(true);
    break;
  case NCR:
    m_vectorOfSeedPointLabels[NCR].radioButton->setChecked(true);
    break;
  case TU:
    m_vectorOfSeedPointLabels[TU].radioButton->setChecked(true);
    break;
  case NE:
    m_vectorOfSeedPointLabels[NE].radioButton->setChecked(true);
    break;
  case CB:
    m_vectorOfSeedPointLabels[CB].radioButton->setChecked(true);
    break;
  case VT:
    m_vectorOfSeedPointLabels[VT].radioButton->setChecked(true);
    break;
  case CAN:
    m_vectorOfSeedPointLabels[CAN].radioButton->setChecked(true);
    break;
  case CAE:
    m_vectorOfSeedPointLabels[CAE].radioButton->setChecked(true);
    break;
  case RTN:
    m_vectorOfSeedPointLabels[RTN].radioButton->setChecked(true);
    break;
  case RTE:
    m_vectorOfSeedPointLabels[RTE].radioButton->setChecked(true);
    break;
  default:
    break;
  }

  double x, y, z;
  x = mCurrentTPoints->mLandmarks[row].coordinates[0];
  y = mCurrentTPoints->mLandmarks[row].coordinates[1];
  z = mCurrentTPoints->mLandmarks[row].coordinates[2];
  emit MoveSlicerCursor(x, y, z);
}

void fTumorPanel::sTableFocused(bool bFocused)
{
  mTumorPointsSelected = true;
}

void fTumorPanel::tTableFocused(bool bFocused)
{
  mTumorPointsSelected = false;
}

void fTumorPanel::SetCurrentSelectedTissueType()
{
  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TISSUE_POINTS, mTissueType, 0);
}
void fTumorPanel::tRemoveAllPoints()
{
  m_typeRadBtnAllTissues->setChecked(true);
  int rowCount = mCurrentTPoints->GetNumberOfPoints();
  for (int rowIndex = rowCount - 1; rowIndex >= 0; rowIndex--)
  {
    mCurrentTPoints->RemoveLandmark(rowIndex);
    tTableWidget->removeRow(rowIndex);
  }
  emit UpdateRenderWindows();
  UpdateCurrentPoints();
}
void fTumorPanel::sRemoveAllPoints()
{
  m_typeRadBtnTumor->setChecked(true);
  int rowCount = mCurrentSPoints->GetNumberOfPoints();
  for (int rowIndex = rowCount - 1; rowIndex >= 0; rowIndex--)
  {
    mCurrentSPoints->RemoveLandmark(rowIndex);
    sTableWidget->removeRow(rowIndex);
  }
  emit UpdateRenderWindows();
}
void fTumorPanel::HighlightCurrentSelctedPoints(double x, double y, double z, double X, double Y, double Z, double value)
{
  if (mCurrentTPoints == NULL)
    return;

  int rowCount = mCurrentTPoints->GetNumberOfPoints();
  for (int rowIndex = 0; rowIndex<rowCount; rowIndex++)
  {
    float roundX = floorf(x * 1000) / 1000;
    float roundY = floorf(y * 1000) / 1000;
    float roundZ = floorf(z * 1000) / 1000;

    float roundCurrentX = floorf(mCurrentTPoints->mLandmarks[rowIndex].coordinates[0] * 1000) / 1000;
    float roundCurrentY = floorf(mCurrentTPoints->mLandmarks[rowIndex].coordinates[1] * 1000) / 1000;
    float roundCurrentZ = floorf(mCurrentTPoints->mLandmarks[rowIndex].coordinates[2] * 1000) / 1000;

    if (roundX == roundCurrentX && roundY == roundCurrentY && roundZ == roundCurrentZ)
    {
      tTableWidget->selectRow(rowIndex);
      QTableWidgetItem * item = tTableWidget->item(rowIndex, 0);
      tTableWidget->scrollToItem(item, QAbstractItemView::PositionAtTop);
      tTableWidget->setItemSelected(item, true);
      break;
    }
  }
}
void fTumorPanel::HandleKeyPressingEventTTable()
{
  if (mCurrentTPoints == NULL)
    return;

  if (mCurrentTPoints->GetNumberOfPoints() == 0)
    return;

  QList<QTableWidgetItem*> items = tTableWidget->selectedItems();
  if (items.empty())
    return;

  int row = items[0]->row();

  double x, y, z;
  x = mCurrentTPoints->mLandmarks[row].coordinates[0];
  y = mCurrentTPoints->mLandmarks[row].coordinates[1];
  z = mCurrentTPoints->mLandmarks[row].coordinates[2];
  emit MoveSlicerCursor(x, y, z);
}
void fTumorPanel::HandleDownKeyEventTTable()
{
  if (mCurrentTPoints == NULL)
    return;

  if (mCurrentTPoints->GetNumberOfPoints() == 0)
    return;

  QList<QTableWidgetItem*> items = tTableWidget->selectedItems();
  if (items.empty())
    return;

  unsigned int rowindex = items[0]->row() + 1;
  if (rowindex > mCurrentTPoints->GetNumberOfPoints())
    return;

  for (unsigned int i = 0; i < mCurrentTPoints->GetNumberOfPoints(); i++)
  {
    QTableWidgetItem * item3 = tTableWidget->item(i, 0);
    tTableWidget->setItemSelected(item3, false);
  }
  tTableWidget->selectRow(rowindex);
  QTableWidgetItem * item2 = tTableWidget->item(rowindex, 0);
  tTableWidget->scrollToItem(item2, QAbstractItemView::PositionAtTop);
  tTableWidget->setItemSelected(item2, true);

}
void fTumorPanel::HandleUpKeyEventTTable()
{
  if (mCurrentTPoints == NULL)
    return;

  if (mCurrentTPoints->GetNumberOfPoints() == 0)
    return;

  QList<QTableWidgetItem*> items = tTableWidget->selectedItems();
  if (items.empty())
    return;

  int rowindex = items[0]->row() - 1;
  if (rowindex < 0)
    return;
  for (unsigned int i = 0; i < mCurrentTPoints->GetNumberOfPoints(); i++)
  {
    QTableWidgetItem * item3 = tTableWidget->item(i, 0);
    tTableWidget->setItemSelected(item3, false);
  }
  tTableWidget->selectRow(rowindex);
  QTableWidgetItem * item2 = tTableWidget->item(rowindex, 0);
  tTableWidget->scrollToItem(item2, QAbstractItemView::PositionAtTop);
  tTableWidget->setItemSelected(item2, true);
}

void fTumorPanel::HandleDeleteKeyEventTTable()
{
  tRemoveSelectedPoints();
}


void fTumorPanel::HandleKeyPressingEventSTable()
{
  if (mCurrentSPoints == NULL)
    return;

  if (mCurrentSPoints->GetNumberOfPoints() == 0)
    return;

  QList<QTableWidgetItem*> items = sTableWidget->selectedItems();
  if (items.empty())
    return;

  int row = items[0]->row();

  double x, y, z;
  x = mCurrentSPoints->mLandmarks[row].coordinates[0];
  y = mCurrentSPoints->mLandmarks[row].coordinates[1];
  z = mCurrentSPoints->mLandmarks[row].coordinates[2];
  emit MoveSlicerCursor(x, y, z);
}
void fTumorPanel::HandleDownKeyEventSTable()
{
  if (mCurrentSPoints == NULL)
    return;

  if (mCurrentSPoints->GetNumberOfPoints() == 0)
    return;

  QList<QTableWidgetItem*> items = sTableWidget->selectedItems();
  if (items.empty())
    return;

  unsigned int rowindex = items[0]->row() + 1;
  if (rowindex > mCurrentSPoints->GetNumberOfPoints())
    return;



  for (unsigned int i = 0; i < mCurrentSPoints->GetNumberOfPoints(); i++)
  {
    QTableWidgetItem * item3 = sTableWidget->item(i, 0);
    sTableWidget->setItemSelected(item3, false);
  }
  sTableWidget->selectRow(rowindex);
  QTableWidgetItem * item2 = sTableWidget->item(rowindex, 0);
  sTableWidget->scrollToItem(item2, QAbstractItemView::PositionAtTop);
  sTableWidget->setItemSelected(item2, true);

}
void fTumorPanel::HandleUpKeyEventSTable()
{
  if (mCurrentSPoints == NULL)
    return;

  if (mCurrentSPoints->GetNumberOfPoints() == 0)
    return;

  QList<QTableWidgetItem*> items = sTableWidget->selectedItems();
  if (items.empty())
    return;

  int rowindex = items[0]->row() - 1;
  if (rowindex < 0)
    return;

  for (unsigned int i = 0; i < mCurrentSPoints->GetNumberOfPoints(); i++)
  {
    QTableWidgetItem * item3 = sTableWidget->item(i, 0);
    sTableWidget->setItemSelected(item3, false);
  }
  sTableWidget->selectRow(rowindex);
  QTableWidgetItem * item2 = sTableWidget->item(rowindex, 0);
  sTableWidget->scrollToItem(item2, QAbstractItemView::PositionAtTop);
  sTableWidget->setItemSelected(item2, true);
}

void fTumorPanel::HandleDeleteKeyEventSTable()
{
  sRemoveSelectedPoints();
}

void fTumorPanel::SetSeedType()
{
  if (mCurrentSPoints == NULL)
  {
    emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TUMOR_POINTS, 0, 0);
  }
  else
  {
    emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TUMOR_POINTS, mCurrentSPoints->GetNumberOfPoints(), 0);
  }

  for (size_t i = 0; i < NumberOfTissueTypes; i++)
  {
    m_vectorOfSeedPointLabels[i].radioButton->setEnabled(false);
  }

}
void fTumorPanel::RadioButtonLabels()
{
  if (m_saveType == SaveType::Generic)
  {
    for (size_t i = 0; i < NumberOfTissueTypes; i++)
    {
      m_vectorOfSeedPointLabels[i].text_label = genericLabels_laconic[i];
      m_vectorOfSeedPointLabels[i].text_toolTip = genericLabels_verbose[i];
    }
  }
  else
  {
    for (size_t i = 0; i < NumberOfTissueTypes; i++)
    {
      m_vectorOfSeedPointLabels[i].text_label = labels_laconic[i];
      m_vectorOfSeedPointLabels[i].text_toolTip = labels_verbose[i];
    }
  }
}

void fTumorPanel::SetTissueType(int inputTissueType)
{
  mTissueType = inputTissueType;
  SetTissueType();
}

void fTumorPanel::SetTissueType()
{
  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TISSUE_POINTS, mTissueType, 0);
  m_saveType = SaveType::allTissues;

  m_vectorOfSeedPointLabels[BG].radioButton->setEnabled(false);

  for (size_t i = 1; i < NumberOfTissueTypes; i++)
  {
    m_vectorOfSeedPointLabels[i].radioButton->setEnabled(true);
  }
  initializationFileName = "initializedPoints_all.txt";

  RadioButtonLabels();
}

void fTumorPanel::SetGLISTRTissueType()
{
  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TISSUE_POINTS, mTissueType, 0);

  m_vectorOfSeedPointLabels[CSF].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[GM].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[WM].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[VS].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[ED].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[NCR].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[TU].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[NE].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[CB].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[VT].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[CAN].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[CAE].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[RTN].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[RTE].radioButton->setEnabled(false);

  initializationFileName = "initializedPoints_glistr.txt";
  m_saveType = SaveType::GLISTR;
  RadioButtonLabels();
}

void fTumorPanel::SetPORTRPRETissueType()
{
  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TISSUE_POINTS, mTissueType, 0);

  m_vectorOfSeedPointLabels[CSF].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[GM].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[WM].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[VS].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[ED].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[NCR].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[TU].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[NE].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[CB].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[VT].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[CAN].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[CAE].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[RTN].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[RTE].radioButton->setEnabled(false);

  initializationFileName = "initializedPoints_portrPre.txt";
  m_saveType = SaveType::PORTRPre;
  RadioButtonLabels();
}

void fTumorPanel::SetPORTRPOSTTissueType()
{
  emit SetActiveLandmarksTypeSignal(LANDMARK_TYPE::TISSUE_POINTS, mTissueType, 0);

  m_vectorOfSeedPointLabels[CSF].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[GM].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[WM].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[VS].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[ED].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[NCR].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[TU].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[NE].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[CB].radioButton->setEnabled(false);
  m_vectorOfSeedPointLabels[VT].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[CAN].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[CAE].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[RTN].radioButton->setEnabled(true);
  m_vectorOfSeedPointLabels[RTE].radioButton->setEnabled(true);

  initializationFileName = "initializedPoints_portrPost.txt";
  m_saveType = SaveType::PORTRPost;
  RadioButtonLabels();
}

void fTumorPanel::UpdateCurrentPoints()
{
  if (mCurrentTPoints == NULL)
    return;

  // initialize counters to zero every time the update is run - done to ensure repition doesn't happen in count 
  for (size_t j = 0; j < NumberOfTissueTypes; j++)
  {
    m_vectorOfSeedPointLabels[j].counter = 0;
  }

  for (unsigned int i = 0; i < mCurrentTPoints->GetNumberOfPoints(); i++)
  {
    switch (mCurrentTPoints->mLandmarks[i].id)
    {
    case BG:
      m_vectorOfSeedPointLabels[BG].counter++;
      break;
    case CSF:
      m_vectorOfSeedPointLabels[CSF].counter++;
      break;
    case GM:
      m_vectorOfSeedPointLabels[GM].counter++;
      break;
    case WM:
      m_vectorOfSeedPointLabels[WM].counter++;
      break;
    case VS:
      m_vectorOfSeedPointLabels[VS].counter++;
      break;
    case ED:
      m_vectorOfSeedPointLabels[ED].counter++;
      break;
    case NCR:
      m_vectorOfSeedPointLabels[NCR].counter++;
      break;
    case TU:
      m_vectorOfSeedPointLabels[TU].counter++;
      break;
    case NE:
      m_vectorOfSeedPointLabels[NE].counter++;
      break;
    case CB:
      m_vectorOfSeedPointLabels[CB].counter++;
      break;
    case VT:
      m_vectorOfSeedPointLabels[VT].counter++;
      break;
    case CAN:
      m_vectorOfSeedPointLabels[CAN].counter++;
      break;
    case CAE:
      m_vectorOfSeedPointLabels[CAE].counter++;
      break;
    case RTN:
      m_vectorOfSeedPointLabels[RTN].counter++;
      break;
    case RTE:
      m_vectorOfSeedPointLabels[RTE].counter++;
      break;
    default:
      cbica::Logging(loggerFile, "Undefined mCurrentTPoints->mLandmarks");
      exit(EXIT_FAILURE);
    }
  }

  for (size_t i = 0; i < NumberOfTissueTypes; i++)
  {
    std::string text = m_vectorOfSeedPointLabels[i].text_label + "::" + std::to_string(m_vectorOfSeedPointLabels[i].counter);
    m_vectorOfSeedPointLabels[i].radioButton->setText(text.c_str());
    m_vectorOfSeedPointLabels[i].radioButton->setToolTip(m_vectorOfSeedPointLabels[i].text_toolTip.c_str());
  }

}

void fTumorPanel::helpClicked()
{
  emit helpClicked_Interaction("gs_seedpoints.html");
}