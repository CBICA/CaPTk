/**
\file  fBottomImageInfoTip.cpp

\brief Implementation of fBottomImageInfoTip class

https://www.med.upenn.edu/sbia/software/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/sbia/software-agreement.html

*/

#include "fBottomImageInfoTip.h"

fBottomImageInfoTip::fBottomImageInfoTip(QWidget * parent) :QWidget(parent)
{
  setupUi(this);
  connect(pixelPosButton, SIGNAL(clicked()), this, SLOT(pixelPosButtonClicked()));
}

void fBottomImageInfoTip::setFileName(QString text)
{
  imageLabel->setText(text);
}

void fBottomImageInfoTip::setSizePixel(double x, double y, double z)
{
  QString text = "";
  text += QString::number((int)x) + "  ";
  text += QString::number((int)y) + "  ";
  text += QString::number((int)z) + "  ";
  text += "  ";
  sizePixelLabel->setText(text);
}

void fBottomImageInfoTip::setOrigin(double x, double y, double z)
{
  QString text = "";
  text += QString::number(x, 'f', 3) + "  ";
  text += QString::number(y, 'f', 3) + "  ";
  text += QString::number(z, 'f', 3) + "  ";
  text += "  ";
  originLabel->setText(text);
}

void fBottomImageInfoTip::setSpacing(double x, double y, double z)
{
  QString text = "";
  text += QString::number(x, 'f', 3) + "  ";
  text += QString::number(y, 'f', 3) + "  ";
  text += QString::number(z, 'f', 3) + "  ";
  text += "  ";
  spacingLabel->setText(text);
}

void fBottomImageInfoTip::setCurrentInfo(int visibility, double x, double y, double z, double X, double Y, double Z, double value)
{
	//! we always want the info to be visible, refer issue #88
	visibility = 1;
	//
	if (visibility)
	{
		pixelPosX->setText(QString::number((int)X));
		pixelPosY->setText(QString::number((int)Y));
		pixelPosZ->setText(QString::number((int)Z));
		valueLabel->setText(QString::number(value));
	}
	else
	{
		pixelPosX->setText(QString(""));
		pixelPosY->setText(QString(""));
		pixelPosZ->setText(QString(""));
		valueLabel->setText(QString(""));
	}
}

void fBottomImageInfoTip::setZSlicePosition(int zslice)
{
  pixelPosZ->setText(QString::number(zslice));
}

void fBottomImageInfoTip::setIntensityValue(double value)
{
  valueLabel->setText(QString::number(value));
}

void fBottomImageInfoTip::pixelPosButtonClicked()
{
  double X, Y, Z;
  X = pixelPosX->text().toDouble();
  Y = pixelPosY->text().toDouble();
  Z = pixelPosZ->text().toDouble();
  emit MoveSlicerCursor(X, Y, Z, 1);
}

