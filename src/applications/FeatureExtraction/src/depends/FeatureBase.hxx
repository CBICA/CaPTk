/**
\file FeatureBase.hxx

This file contains the definitions for FeatureBase.

https://www.med.upenn.edu/cbica/captk/ <br>
software@cbica.upenn.edu

Copyright (c) 2018 University of Pennsylvania. All rights reserved. <br>
See COPYING file or https://www.med.upenn.edu/cbica/software-agreement.html
*/
#pragma once

#include "FeatureBase.h"

template< class TImage >
void FeatureBase< TImage >::SetInputImage(const typename TImage::Pointer image)
{
  m_inputImage = image;
  m_algorithmDone = false;
}

template< class TImage >
void FeatureBase< TImage >::SetInputMask(const typename TImage::Pointer image)
{
  m_Mask = image;
  TConstIteratorType maskIterator(image, image->GetBufferedRegion());
  for (maskIterator.GoToBegin(); !maskIterator.IsAtEnd(); ++maskIterator)
  {
    if (maskIterator.Get() > 0)
    {
      m_nonZeroIndeces.push_back(maskIterator.GetIndex());
    }
  }
  m_algorithmDone = false;
}

template< class TImage >
void FeatureBase< TImage >::SetNonZeroIndeces(std::vector< typename TImage::IndexType > &nonZeroIndeces)
{
  m_nonZeroIndeces = nonZeroIndeces;

  if (!m_inputImage.IsNull())
  {
    m_Mask = cbica::CreateImage< TImage >(m_inputImage);

    TIteratorType maskIt(m_Mask, m_Mask->GetBufferedRegion());
    for (size_t i = 0; i < m_nonZeroIndeces.size(); i++)
    {
      maskIt.SetIndex(m_nonZeroIndeces[i]);
      maskIt.Set(1);
    }
  }

  m_algorithmDone = false;
}