#pragma once

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>

#if defined(_WIN32)&& !defined(__GNUC__)
  #pragma warning( disable : 4503 )// Warnig related to decorated name length size
#endif

////TBD add documentation once class is frozen
//inline std::vector<std::string> splitStr(const std::string& str, const std::string& delim)
//{
//  std::vector<std::string> tokens;
//  size_t prev = 0, pos = 0;
//  do
//  {
//    pos = str.find(delim, prev);
//    if (pos == std::string::npos) pos = str.length();
//    std::string token = str.substr(prev, pos - prev);
//    if (!token.empty()) tokens.push_back(token);
//    prev = pos + delim.length();
//  } while (pos < str.length() && prev < str.length());
//  return tokens;
//}
//inline void removeChar(std::string& str, char ch)
//{
//  str.erase(remove(str.begin(), str.end(), ch), str.end());
//}
class FeatureParser
{
public:
  FeatureParser(const std::string& fileName)
  {
    m_fileName = fileName;
    m_featureMap = prepareFeatureMap();
  }
  ~FeatureParser()
  {
  }
  std::string getLastError()
  {
    return m_lastError;
  }
  std::map<std::string, std::vector<std::map<std::string, std::string> > > getFeatureMap()
  {
    return m_featureMap;
  }
  void setFeatureMap(std::map<std::string, std::vector<std::map<std::string, std::string> > > featureMap)
  {
    m_featureMap = featureMap;
  }
private:
  std::map<std::string, std::vector<std::map<std::string, std::string> > > prepareFeatureMap()
  {
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > matrix;
    readCsv(m_fileName, header, matrix);
    std::map< std::string, std::vector<std::map<std::string, std::string> > > features;
    for (auto row : matrix)
    {
      std::map<std::string, std::string> props;
      for (size_t i = 1; i < row.size(); i++)
      {
        props[header[i]] = row[i];
        if (header[i] == "Default")
        {
          props["Value"] = row[i];
        }
      }
      features[row[0]].push_back(props);
    }
    return features;
  }
  std::string trim(const std::string &s)
  {
    auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c){return std::isspace(c); });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c){return std::isspace(c); }).base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
  }
  bool readCsv(const std::string& fileName, std::vector<std::string>& header, std::vector<std::vector<std::string> >&   matrix)
  {
    std::ifstream file(fileName.c_str());
    if (!file)
    {
      m_lastError = "Could not read: " + fileName;
      return false;
    }
    while (file)
    {
      std::string line;
      getline(file, line);
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
      else if(!row.empty())
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
  std::map<std::string, std::vector<std::map<std::string, std::string> > > m_featureMap;
};

