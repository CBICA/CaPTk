/**
\file  fBottomImageInfoTip.h

\brief Declaration of fBottomImageInfoTip class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#ifndef _fBottomImageInfoTip_h_
#define _fBottomImageInfoTip_h_


//#include "CAPTk.h"
#include "ui_fBottomImageInfoTip.h"

/**
\class fBottomImageInfoTip

\brief This class controls the elements in the bottom panel of the UI
*/
class fBottomImageInfoTip : public QWidget, private Ui::fBottomImageInfoTip
{
  Q_OBJECT

public:
  //!Constructor
  fBottomImageInfoTip(QWidget* parent = 0);

  //!Destructor
  ~fBottomImageInfoTip() {}

  /**
  \brief Display the name of currently displayed image
  \param test Image name
  */
  void setFileName(QString text);

  /**
  \brief Display the origin of currently displayed image
  \param x x-coordinate of the origin
  \param y y-coordinate of the origin
  \param z z-coordinate of the origin
  */
  void setOrigin(double x, double y, double z);

  /**
  \brief Display the spacing of currently displayed image
  \param x Spacing in x-direction
  \param y Spacing in y-direction
  \param z Spacing in z-direction
  */
  void setSpacing(double x, double y, double z);

  /**
  \brief Display the size of currently displayed image
  \param x Size in x-direction
  \param y Size in y-direction
  \param z Size in z-direction
  */
  void setSizePixel(double x, double y, double z);

  /**
  \brief Display the information of currently selected voxel
  \param visibility Whether information needs to be displayed or not
  \param x World position x
  \param y World position y
  \param z World position z
  \param X Pixel position x
  \param Y Pixel position y
  \param Z Pixel position z
  \param value Intensity value of the voxel
  */
  void setCurrentInfo(int visibility, double x, double y, double z, double X, double Y, double Z, double value);

  public slots:
  /**
  \brief Move the cursor to given pixel position 
  */
  void pixelPosButtonClicked();

  /**
 \brief Display Z slice position 
 \param z slize position
 */
  void setZSlicePosition(int zslice);

  /**
\brief Display intensity value at cursor position
\param value intensity value
*/
  void setIntensityValue(double value);

signals:
  void MoveSlicerCursor(double, double, double, int);
};


#endif
