/*=========================================================================

  Program:   ITK General
  Module:    $RCSfile: StackSlices.cxx,v $
  Language:  C++
  Date:      $Date: 2010/04/04 11:23:11 $
  Version:   $Revision: 1.22 $

  Author: Yuanjie Zheng (zheng.vision@gmail.com)
  Institution: PICSL

  This software is distributed WITHOUT ANY WARRANTY; without even 
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
  PURPOSE.

=========================================================================*/
#ifndef __AUXILIARYFUNCTIONS_h_
#define __AUXILIARYFUNCTIONS_h_

#include <float.h>//DBL_MAX


bool isinf(double x)
{
	if(x>DBL_MAX)
		return true;
	else
		return false;
}

#endif // __AUXILIARYFUNCTIONS_h_
