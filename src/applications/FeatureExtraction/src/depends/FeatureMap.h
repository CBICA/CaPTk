/**
\file FeatureMap.h

This file extracts the feature parameters from the parameter file.

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software/license.html
*/
#pragma once

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>
#include <mutex>

#if defined(_WIN32)&& !defined(__GNUC__)
  #pragma warning( disable : 4503 )// Warnig related to decorated name length size
#endif

//! Splits the string
inline std::vector<std::string> splitStr(const std::string& str, const std::string& delim)
{
  std::vector<std::string> tokens;
  size_t prev = 0, pos = 0;
  do
  {
    pos = str.find(delim, prev);
    if (pos == std::string::npos) pos = str.length();
    std::string token = str.substr(prev, pos - prev);
    if (!token.empty()) tokens.push_back(token);
    prev = pos + delim.length();
  } while (pos < str.length() && prev < str.length());
  return tokens;
}

////! Removes the specified character from string
//inline void removeChar(std::string& str, char ch)
//{
//  str.erase(remove(str.begin(), str.end(), ch), str.end());
//}

/**
\brief This is a helper class used to parse the parameter file passed to FeatureExtraction
*/
class FeatureMap
{
public:
  using Container = std::map< std::string, std::map< std::string, std::string > >; //! helper alias for a big data structure

  FeatureMap(const std::string& fileName, std::string selected_features)
  {
    m_fileName = fileName;
    prepareFeatureMap(selected_features);
  }
  ~FeatureMap()
  {
  }
  std::string getLastError()
  {
    return m_lastError;
  }
  Container getFeatureMap()
  {
    return m_featureParams;
  }
  void setFeatureMap(Container featureMap)
  {
    m_featureParams = featureMap;
  }
private:

  Container m_featureParams;
  
  void prepareFeatureMap(std::string selected_features)
  {
    std::vector< std::string > header;
    std::vector< std::vector< std::string > > matrix;
    readCsv(m_fileName, header, matrix);

    for (auto row : matrix)
    {
      if (!selected_features.empty())
      {
        if (selected_features == row[0])
        {
          std::map<std::string, std::string> props;
          for (size_t i = 2; i < row.size(); i++)
          {
            props[header[i]] = row[i];
            if (header[i] == "Default")
            {
              props["Value"] = row[i];
            }
          }

          m_featureParams[row[1]] = props;
        }
      }
    }
   
  }
  std::string trim(const std::string &s)
  {
    auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c){return std::isspace(c); });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c){return std::isspace(c); }).base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
  }
  bool readCsv(const std::string& fileName, std::vector<std::string>& header, std::vector<std::vector<std::string> >&   matrix)
  {
    std::mutex file_mutex;
    const std::lock_guard< std::mutex > lock(file_mutex); // this will release the lock when the function closes
    std::ifstream file(fileName.c_str());
    if (!file)
    {
      m_lastError = "Could not read: " + fileName;
      return false;
    }
    while (file)
    {
      std::string line;
      std::getline(file, line);
      // fix line ending problems 
      std::remove_copy(line.begin(), line.end(), line.begin(), '\r');
      std::stringstream lineStream(line);
      std::vector<std::string> row;
      std::string cell;
      while (getline(lineStream, cell, ','))
      {
        row.push_back(trim(cell));
      }
      if (!row.empty() && header.empty())
      {
        header = row;
      }
      else if (!row.empty())
      {
        row.resize(header.size());
        if (row[0].empty() && !matrix.empty())
        {
          row[0] = matrix.back()[0]; //Copy previous feature family name if empty
        }
        matrix.push_back(row);
      }
    }
    if (header.empty() || matrix.empty())
    {
      m_lastError = "Empty file: " + fileName;
      return false;
    }
    return true;
  }

  std::string m_lastError;
  std::string m_fileName;

 // std::map<std::string, std::vector<std::map<std::string, std::string> > > m_featureMap;
};

