/*=========================================================================

  Program:   ALFABIS fast medical image registration programs
  Language:  C++
  Website:   github.com/pyushkevich/greedy
  Copyright (c) Paul Yushkevich, University of Pennsylvania. All rights reserved.

  This program is part of ALFABIS: Adaptive Large-Scale Framework for
  Automatic Biomedical Image Segmentation.

  ALFABIS development is funded by the NIH grant R01 EB017255.

  ALFABIS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ALFABIS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ALFABIS.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/
#ifndef __GreedyException_h_
#define __GreedyException_h_

#include <cstdio>
#include <cstdarg>

/**
 * A simple exception class with string formatting
 */
class GreedyException : public std::exception
{
public:

  GreedyException(const char *format, ...)
    {
    buffer = new char[4096];
    va_list args;
    va_start (args, format);
    //vsprintf (buffer,format, args);
    va_end (args);
    }

  virtual const char* what() const throw() { return buffer; }

  virtual ~GreedyException() throw() { delete buffer; }

private:

  char *buffer;

};


#endif
