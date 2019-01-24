/**
\file CAPTk.cpp

\brief Main entry point for CAPTk
*/

#include <QApplication>
#include "fMainWindow.h"
//#include "CAPTk.h"

#include "vtkFileOutputWindow.h"
#include "itkFileOutputWindow.h"
//#include "vtkOutputWindow"

#include "cbicaCmdParser.h"
#include "cbicaUtilities.h"
#include "yaml-cpp/yaml.h"

///// debug
//#define _CRTDBG_MAP_ALLOC
//#include <stdlib.h>
//#include <crtdbg.h>
//#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )      
//#define new DBG_NEW   
///// debug

void setStyleSheet(const std::string &styleFileName = CAPTK_STYLESHEET_FILE)
{
  std::string styleFileFullPath = getCaPTkDataDir() + "../etc/" + CAPTK_STYLESHEET_FILE;

  QFile f(styleFileFullPath.c_str());
  if (!f.exists())
  {
#if WIN32
    styleFileFullPath = captk_currentApplicationPath + "../../data/" + CAPTK_STYLESHEET_FILE;
#else
    styleFileFullPath = captk_currentApplicationPath + "../data/" + CAPTK_STYLESHEET_FILE;
#endif
    if (!f.exists())
    {
      cbica::Logging(loggerFile, "Unable to set stylesheet, file not found");
    }
  }
  else
  {
    f.open(QFile::ReadOnly | QFile::Text);
    QTextStream ts(&f);
    qApp->setStyleSheet(ts.readAll());
    f.close();
  }
}

#ifdef _WIN32
// ensures no console pops up when launching the program
int CALLBACK WinMain(HINSTANCE hInstance, HINSTANCE prevInstance, LPSTR lpCmdLine, int nShowCmd)
{
  int argc = __argc;
  char **argv = __argv;
#else
int main(int argc, char** argv)
{
#endif

  QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());

  #if __APPLE__
  // this->
  QApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);
  #endif

  //! Support for High DPI monitors..works on windows but still some menu issues are seen
  //! Needs to be tested on Linux and Mac
  QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
  //QCoreApplication::setAttribute(Qt::AA_UseOpenGLES);

  QApplication app(argc, argv);

  //cbica::setEnvironmentVariable("QT_QPA_PLATFORM_PLUGIN_PATH", captk_currentApplicationPath + "/platforms");
  cbica::setEnvironmentVariable("QT_OPENGL", "software");
  
  // parse the command line
  auto parser = cbica::CmdParser(argc, argv, "CaPTk");
  parser.ignoreArgc1();

  parser.addOptionalParameter("i", "images", cbica::Parameter::FILE, "NIfTI or DICOM", "Input coregistered image(s) to load into CaPTk", "Multiple images are delineated using ','");
  parser.addOptionalParameter("m", "mask", cbica::Parameter::FILE, "NIfTI or DICOM", "Input mask [coregistered with image(s)] to load into CaPTk", "Accepts only one file");
  parser.addOptionalParameter("tu", "tumorPt", cbica::Parameter::FILE, ".txt", "Tumor Point file for the image(s) being loaded");
  parser.addOptionalParameter("ts", "tissuePt", cbica::Parameter::FILE, ".txt", "Tissue Point file for the image(s) being loaded");
  parser.addOptionalParameter("a", "advanced", cbica::Parameter::BOOLEAN, "none", "Advanced visualizer which does *not* consider", "origin information during loading");
  
  // parser.addOptionalParameter("de", "direction", cbica::Parameter::STRING, "", "Calls Directionality Estimator CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("eg", "egfrviii", cbica::Parameter::STRING, "", "Calls EGFRvIII PHI Calculator CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("fe", "feature", cbica::Parameter::STRING, "", "Calls Feature Extractor CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("ge", "geodesic", cbica::Parameter::STRING, "", "Calls Geodesic Segmentation CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("is", "imgsub", cbica::Parameter::STRING, "", "Calls Imaging SubType CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("ms", "molsub", cbica::Parameter::STRING, "", "Calls Molecular SubType CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("pa", "population", cbica::Parameter::STRING, "", "Calls Population Atlas CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("re", "recurrence", cbica::Parameter::STRING, "", "Calls Recurrence Estimator CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("ss", "sbrtSeg", cbica::Parameter::STRING, "", "Calls SBRT Segmentation CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("sa", "sbrtAna", cbica::Parameter::STRING, "", "Calls SBRT Analyze CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("su", "survival", cbica::Parameter::STRING, "", "Calls Survival Predictor CLI", "This needs to be passed first (before any images/masks)");
  // parser.addOptionalParameter("ws", "whites", cbica::Parameter::STRING, "", "Calls WhiteStripe CLI", "This needs to be passed first (before any images/masks)");
  
  parser.exampleUsage("-i C:/data/input1.nii.gz,C:/data/input2.nii.gz -m C:/data/inputMask.nii.gz -tu C:/data/init_seed.txt -ts C:/data/init_GLISTR.txt");

  std::string cmd_inputs, cmd_mask, cmd_tumor, cmd_tissue;
  
  if (parser.isPresent("i"))
  {
    parser.getParameterValue("i", cmd_inputs);
  }
  if (parser.isPresent("m"))
  {
    parser.getParameterValue("m", cmd_mask);
  }
  if (parser.isPresent("tu"))
  {
    parser.getParameterValue("tu", cmd_tumor);
  }
  if (parser.isPresent("ts"))
  {
    parser.getParameterValue("ts", cmd_tissue);
  }

  // check for CWL command coming in through the command line after "CaPTk"
  if (cmd_inputs.empty() && (argc > 1))
  {
    auto cwlFiles = cbica::getCWLFilesInApplicationDir();
    for (auto & file : cwlFiles)
    {
      // Check for filename without cwl extension
      if (cbica::getFilenameBase(file).find(std::string(argv[1])) != std::string::npos) 
      {
        // Get base command
        //std::ofstream selected_file;
        //selected_file.open(file.c_str());
        auto config = YAML::LoadFile(file);
        // Get all args passed to application
        std::string argv_complete;        
        for (size_t i = 2; i < argc; i++) // 2 because the argv[1] is always the "application"
        {
          argv_complete += " " + std::string(argv[i]);
        }
        // Pass them in
        return std::system((getApplicationPath(config["baseCommand"].as<std::string>()) + argv_complete).c_str());
      }
    }
  }

  ///// debug
  //HANDLE hLogFile;

  //hLogFile = CreateFile("MemoryLeaks.txt", GENERIC_WRITE, FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
  //
  //_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
  //_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  //_CrtSetReportFile(_CRT_WARN, hLogFile);
  //_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  //_CrtSetReportFile(_CRT_ERROR, hLogFile);
  //_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);
  //_CrtSetReportFile(_CRT_ASSERT, hLogFile);
  ///// debug

  // if (parser.isPresent("de"))
  // {
  //   int temp;
  //   parser.compareParameter("de", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("DirectionalityEstimate") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("eg"))
  // {
  //   int temp;
  //   parser.compareParameter("eg", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("EGFRvIIISurrogateIndex") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("fe"))
  // {
  //   int temp;
  //   parser.compareParameter("fe", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("FeatureExtraction") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("ge"))
  // {
  //   int temp;
  //   parser.compareParameter("ge", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("GeodesicSegmentation") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("is"))
  // {
  //   int temp;
  //   parser.compareParameter("is", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("ImagingSubtypePredictor") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("ms"))
  // {
  //   int temp;
  //   parser.compareParameter("ms", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("MolecularSubtypePredictor") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("pa"))
  // {
  //   int temp;
  //   parser.compareParameter("pa", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("PopulationAtlases") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("re"))
  // {
  //   int temp;
  //   parser.compareParameter("re", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("RecurrenceEstimator") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("ss"))
  // {
  //   int temp;
  //   parser.compareParameter("ss", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("SBRT_Lung_Segment") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("sa"))
  // {
  //   int temp;
  //   parser.compareParameter("sa", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("SBRT_Lung_Analyze") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("su"))
  // {
  //   int temp;
  //   parser.compareParameter("su", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("SurvivalPredictor") + argv_complete).c_str());
  // }
  // else if (parser.isPresent("ws"))
  // {
  //   int temp;
  //   parser.compareParameter("ws", temp);
  //   std::string argv_complete;
  //   for (size_t i = temp + 1; i < argc; i++)
  //   {
  //     argv_complete = argv_complete + " " + std::string(argv[i]);
  //   }
  //   return std::system((getApplicationPath("WhiteStripe") + argv_complete).c_str());
  // }

  //vtkOpenGLRenderWindow::SetGlobalMaximumNumberOfMultiSamples(0);

  //auto defaultFormat = QVTKOpenGLWidget::defaultFormat();
  //defaultFormat.setSamples(0);
  //QSurfaceFormat::setDefaultFormat(defaultFormat);

  //VTK_MODULE_INIT(vtkRenderingFreeType);

  //! redirect the vtkoutputwindow contents to file
  auto fileOutputWindow = vtkSmartPointer< vtkFileOutputWindow >::New();
  fileOutputWindow->SetFileName((loggerFolderBase + "vtk_errors.txt").c_str());
  vtkOutputWindow::SetInstance(fileOutputWindow);

  //! redirect the itk output window contents to file
  auto itkOutputWindow = itk::FileOutputWindow::New();
  itkOutputWindow->SetFileName((loggerFolderBase + "itk_errors.txt").c_str());
  itkOutputWindow->FlushOn();
  itk::OutputWindow::SetInstance(itkOutputWindow);

  fMainWindow window; // initialize main app
  if (parser.isPresent("a"))
  {
    window.EnableAdvancedVisualizer();
  }

  // get everything to QString
  std::vector< QString > inputFiles_QString;
  QString inputMask;

  if (!cmd_inputs.empty()) // the check is only for input images since images can be rendered with ROI but not the other way round
  {
    auto inputFiles = cbica::stringSplit(cmd_inputs, ",");
    auto inputMasks = cbica::stringSplit(cmd_mask, ",");

    for (size_t i = 0; i < inputFiles.size(); i++)
    {
      if (cbica::fileExists(inputFiles[i])) // check for file validity
      {
        inputFiles_QString.push_back(inputFiles[i].c_str());
      }
    }

    if (inputFiles_QString.size() != 0)
    {
      for (size_t i = 0; i < inputMasks.size(); i++)
      {
        if (cbica::fileExists(inputMasks[i])) // check for file validity
        {
          inputMask = inputMasks[i].c_str();
          break; // because only 1 annotated ROI can be drawn up at one time
        }
      }
    }
    else // at this point, give an error if there are no valid image files passed
    {
      std::cerr << "No valid images were provided.\n";
      return EXIT_FAILURE;
    }
  }

  cbica::createDir(loggerFolderBase);
  cbica::createDir(loggerFolder);
  cbica::createDir(captk_StuffFolderBase);
  cbica::createDir(captk_SampleDataFolder);
  cbica::createDir(captk_PretrainedFolder);

  setStyleSheet();
#ifndef _WIN32
  std::string old_locale = setlocale(LC_NUMERIC, NULL);
  setlocale(LC_NUMERIC, "POSIX");
#endif

  app.processEvents();

  if (inputFiles_QString.size() == 0)
  {
#ifndef CAPTK_PACKAGE_PROJECT
    std::string testData = "/data/AAAC0_flair_pp_shrunk.nii.gz";
    std::string testData_mask = "/data/AAAC0_flair_pp_shrunk_testTumor.nii.gz";

    // if test data is present, load it up otherwise forget about it
    std::string appDir = captk_currentApplicationPath;
    if (cbica::fileExists(appDir + testData))
    {
      inputFiles_QString.push_back((appDir + testData).c_str());
      inputMask.push_back((appDir + testData_mask).c_str());
    }
    else if (cbica::fileExists(appDir + "/.." + testData))
    {
      inputFiles_QString.push_back((appDir + "/.." + testData).c_str());
      inputMask.push_back((appDir + "/.." + testData_mask).c_str());
    }
    else if (cbica::fileExists(appDir + "/../.." + testData))
    {
      inputFiles_QString.push_back((appDir + "/../.." + testData).c_str());
      inputMask.push_back((appDir + "/../.." + testData_mask).c_str());
    }
#endif
  }
  else
    window.loadFromCommandLine(inputFiles_QString, inputMask.toStdString(), cmd_tumor, cmd_tissue/*, true*/); // at this point, inputFiles_QString will have at least 1 value



  //window.setFixedSize(QSize(2735, 1538)); // useful when doing video recording from SP4 and maintain 16x9 ratio
  window.showMaximized();
  window.show();

  // show the "about" screen in the first run
  if (!cbica::fileExists(tutorialScreen))
  {
      auto rec = QApplication::desktop()->screenGeometry();
      // std::cout << "Detected Size: " << rec.width() << "x" << rec.height() << "\n";
    window.about();
  }
  

  
#ifndef _WIN32
  setlocale(LC_NUMERIC, old_locale.c_str());
#endif

  return app.exec();
}

