#include "confetti.h"
#include "multinomialModel.h"
#include "itkImage.h"
#include <ctime>

//#include <opencv2/opencv.hpp>
// #include "opencv2/core.hpp"

//------------------Utility functions----------------------------------
bool icasecmp(const string& l, const string& r)
{
  return l.size() == r.size()
    && equal(l.cbegin(), l.cend(), r.cbegin(),
    [](string::value_type l1, string::value_type r1)
  { return toupper(l1) == toupper(r1); });
}
bool cmdOptionExists(const vector<string>& allArgs, const string& option)
{
  for (auto &arg : allArgs)
  {
    if (icasecmp(arg, option))
    {
      return true;
    }
  }
  return false;
}
static string getCmdOption(const vector<string>& allArgs, const std::string & option)
{
  for (size_t i = 0; i < allArgs.size() - 1;i++)
  {
    if (icasecmp(allArgs[i], option))
    {
      return allArgs[i + 1];
    }
  }
  return string();
}
inline bool fileExist(const std::string& name) {
  ifstream f(name.c_str());
  return f.good();
}
inline bool checkExtn(const std::string& name, const std::string& extn) {
  return  (name.substr(name.find_last_of(".") + 1) == extn);
}
inline bool isFilesEqual(const string& file1, const string& file2)
{
  ifstream in1(file1.c_str(), std::ifstream::in);
  ifstream in2(file2.c_str(), std::ifstream::in);
  if (!in1.good() || !in2.good())
  {
    return false;
  }
  ifstream::pos_type size1, size2;
  size1 = in1.seekg(0, ifstream::end).tellg();
  in1.seekg(0, ifstream::beg);
  size2 = in2.seekg(0, ifstream::end).tellg();
  in2.seekg(0, ifstream::beg);
  if (size1 != size2)
  {
    return false;
  }
  static const size_t BLOCKSIZE = 4096;
  size_t remaining = size1;
  while (remaining)
  {
    char buffer1[BLOCKSIZE], buffer2[BLOCKSIZE];
    size_t size = std::min(BLOCKSIZE, remaining);
    in1.read(buffer1, size);
    in2.read(buffer2, size);

    if (0 != memcmp(buffer1, buffer2, size))
    {
      return false;
    }
    remaining -= size;
  }
  return true;
}

//------------------Main functions----------------------------------
void version()
{
  cout << "ConfettiCore version 1.1" << endl;
}
inline void copyrightNotice()
{
  cout <<
    "\n==========================================================================\n" <<
    "Contact: software@cbica.upenn.edu\n\n" <<
    "Copyright (c) 2017 University of Pennsylvania. All rights reserved.\n" <<
    "See COPYING file or https://www.med.upenn.edu/cbica/captk/license.html" <<
    "\n==========================================================================\n";
}
void help()
{
  cout << "Help: 3  Individual steps of Confetti and their command line usage is given below" << endl;
  cout << "[Genrate Connectivity Signature]: signature -i <tdi_paths.csv> -f <fiberFile.Bfloat> -o <outputSignature.csv>" << endl;
  cout << "[Genrate Cluster]: cluster -s <signatures.csv> -k <200> -t <inTemplateDir> -o <clusterIDs.csv> -w <outTemplateDir>" << endl;
  cout << "[Extract Tracks]: extract -t <templateDir> -f <fiberFile.Bfloat>  -c <clusterIDs.csv> -o <outputDir>" << endl;
}
int runConfettiCmd(const vector<string>& allArgs)
{
  std::clock_t    start;
  start = std::clock();
  try
  {
    if (cmdOptionExists(allArgs, "signature"))
    {

      string tdiCsvFile = getCmdOption(allArgs, "-i");
      string signatureFileName = getCmdOption(allArgs, "-o");
      string fiberFile = getCmdOption(allArgs, "-f");
      if (!fileExist(tdiCsvFile))
      {
        cout << "Error: cannot find tdiCsvFile:" + tdiCsvFile << endl;
        return EXIT_FAILURE;
      }
      if (!fileExist(fiberFile))
      {
        cout << "Error: cannot find fiberFile:" + fiberFile << endl;
        return EXIT_FAILURE;
      }
      if (!checkExtn(signatureFileName, "csv"))
      {
        cout << "Error: require signature file as .csv:" + signatureFileName << endl;
        return EXIT_FAILURE;
      }
      Confetti obj;
      return obj.genConnectivitySig( tdiCsvFile, signatureFileName, fiberFile);
    }
    else if (cmdOptionExists(allArgs, "cluster"))
    {
      string signatureFile = getCmdOption(allArgs, "-s");
      string templateDir = getCmdOption(allArgs, "-t");
      string outFile = getCmdOption(allArgs, "-o");
      bool bUseTemplate = cmdOptionExists(allArgs,"-t");
      string outTemplateDir = getCmdOption(allArgs, "-w");

      int K = 200;
      if (cmdOptionExists(allArgs,"-k")) 
      {
        K = stoi(getCmdOption(allArgs, "-k"));
      }
      if (!fileExist(signatureFile))
      {
        cout << "Error: cannot find signatureFile:" + signatureFile << endl;
        return EXIT_FAILURE;
      }
      if (bUseTemplate)
      {
        if (!fileExist(templateDir + "/template_Alp.csv"))
        {
          cout << "Error: cannot find template:" + templateDir + "/template_Alp.csv" << endl;
          return EXIT_FAILURE;
        }
        if (!fileExist(templateDir + "/template_Pr.csv"))
        {
          cout << "Error: cannot find template:" + templateDir + "/template_Pr.csv" << endl;
          return EXIT_FAILURE;
        }
      }
      if (!outTemplateDir.empty())
      {
        if (!Utils::dirExists(outTemplateDir))
        {
          cout << "Error: output template directory doesn't exist:" + outTemplateDir << endl;
          return EXIT_FAILURE;
        }
      }

      if (!checkExtn(outFile, "csv"))
      {
        cout << "Error: require outFile as .csv:" + outFile << endl;
        return EXIT_FAILURE;
      }
      MNM mm;
      return mm.genAdaptiveCluster(signatureFile, templateDir, outFile, K, bUseTemplate, outTemplateDir);
    }
    else if (cmdOptionExists(allArgs, "extract"))
    {

      string outDIR = getCmdOption(allArgs, "-o");
      string templateDir = getCmdOption(allArgs, "-t");
      string inputLabelFile = getCmdOption(allArgs, "-c");
      string inputBfloatFile = getCmdOption(allArgs, "-f");
      if (!fileExist(inputLabelFile))
      {
        cout << "Error: cannot find inputLabelFile:" + inputLabelFile << endl;
        return EXIT_FAILURE;
      }
      if (!fileExist(inputBfloatFile))
      {
        cout << "Error: cannot find fiberFile:" + inputBfloatFile << endl;
        return EXIT_FAILURE;
      }

      string inputIDFile = templateDir + "/tract_ids.csv";
      string inputNamesFile = templateDir + "/tract_names.csv";
      if (!fileExist(inputIDFile))
      {
        cout << "Error: cannot find inputIDFile:" + inputIDFile << endl;
        return EXIT_FAILURE;
      }
      if (!fileExist(inputNamesFile))
      {
        cout << "Error: cannot find inputNamesFile:" + inputNamesFile << endl;
        return EXIT_FAILURE;
      }
      Confetti obj;
      return obj.extractTract(outDIR, templateDir, inputIDFile, inputNamesFile, inputLabelFile, inputBfloatFile);
    }
    else
    {
      cout << "No valid arguments provided!" << endl;
      help();
    }
    std::cout << "Time: " << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;
  }
  catch (...)
  {
    std::cout << "Error in runConfettiCmd" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
int runConfettiTest(const vector<string>& allArgs)
{
  string dataDir = "../Data/";
  if (!fileExist(dataDir + "input/87_labels.csv"))
  {
    dataDir = "Data/";
  }
  try
  {

    if (cmdOptionExists(allArgs, "Connectivity"))
    {
      string tdiCsvFile = dataDir + "input/87_labels.csv";
      string signatureFileName = dataDir + "output/signaturesSmall.csv";
      string signatureFileNameGT = dataDir + "output/signaturesSmallGT.csv";
      string fiberFile = dataDir + "input/fiberSmall.Bfloat";
      vector<string> testArgs;
      testArgs.push_back(allArgs[0]);
      testArgs.push_back("signature");
      testArgs.push_back("-i");
      testArgs.push_back(tdiCsvFile);
      testArgs.push_back("-f");
      testArgs.push_back(fiberFile);
      testArgs.push_back("-o");
      testArgs.push_back(signatureFileName);
      if (EXIT_SUCCESS == runConfettiCmd(testArgs))
      {
        if (isFilesEqual(signatureFileName, signatureFileNameGT))
        {
          return EXIT_SUCCESS;
        }
      }
      return EXIT_FAILURE;
    }
    else if (cmdOptionExists(allArgs, "Clustering"))
    {
      string signatureFile = dataDir + "output/signaturesSmall.csv";
      string templateDir = dataDir + "input/bundle_template_20150306";
      string outTemplateDir = dataDir + "output/outTemplate";
      string outFile = dataDir + "output/clusterIDsSmall.csv";
      string outFileGT = dataDir + "output/clusterIDsSmallGT.csv";
      vector<string> testArgs;
      testArgs.push_back(allArgs[0]);
      testArgs.push_back("cluster");
      testArgs.push_back("-s");
      testArgs.push_back(signatureFile);
      testArgs.push_back("-t");
      testArgs.push_back(templateDir);
      testArgs.push_back("-o");
      testArgs.push_back(outFile);
      testArgs.push_back("-w");
      testArgs.push_back(outTemplateDir);
      if (EXIT_SUCCESS == runConfettiCmd(testArgs))
      {
        if (isFilesEqual(outFile, outFileGT))
        {
          return EXIT_SUCCESS;
        }
      }
      return EXIT_FAILURE;
    }
    else if (cmdOptionExists(allArgs, "Extraction"))
    {
      string outDIR = dataDir + "output/tract";
      string templateDir = dataDir + "input/bundle_template_20150306";
      string inputLabelFile = dataDir + "output/clusterIDsSmall.csv";
      string fiberFile = dataDir + "input/fiberSmall.Bfloat";
      vector<string> testArgs;
      testArgs.push_back(allArgs[0]);
      testArgs.push_back("extract");
      testArgs.push_back("-o");
      testArgs.push_back(outDIR);
      testArgs.push_back("-t");
      testArgs.push_back(templateDir);
      testArgs.push_back("-c");
      testArgs.push_back(inputLabelFile);
      testArgs.push_back("-f");
      testArgs.push_back(fiberFile);

      return runConfettiCmd(testArgs);
    }

  }
  catch (...)
  {
    std::cout << "Error in runConfettiCmd" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
void debugCode()
{
  vector<string> allArgs;
  allArgs.push_back("Extraction");
  runConfettiTest(allArgs);
  exit(0);
}
int main(int argc, char *argv[]) {
  //string imageFile = "../Data/input/100001_t0.connectivity/seeds_to_ROI10.nii.gz";
  //string fiberFile = "../Data/input/fiberSmall.Bfloat";
  //Confetti().createTDI(fiberFile, imageFile);

  copyrightNotice();//TBD remove in apps if it is not necessary 
  std::vector<string> allArgs(argv, argv + argc);
  if (cmdOptionExists(allArgs, "-v") || cmdOptionExists(allArgs, "-version"))
  {
    version();
    return EXIT_SUCCESS;
  }
  if (cmdOptionExists(allArgs, "-h") || cmdOptionExists(allArgs, "-help"))
  {
    help();
    return EXIT_SUCCESS;
  }
  if (cmdOptionExists(allArgs, "-test"))
  {
    return runConfettiTest(allArgs);
  }
  else
  {
    return runConfettiCmd(allArgs);
  }
}
