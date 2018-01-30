/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkSignalIntensityToS0ImageFilter.hxx$
  Language:  C++
  Date:      $Date: 2012/03/7 $
  Version:   $Revision: 1.0 $

  =========================================================================*/
#ifndef _itkSignalIntensityToS0ImageFilter_hxx
#define _itkSignalIntensityToS0ImageFilter_hxx

#include "itkSignalIntensityToS0ImageFilter.h"

#include "SignalComputationUtils.h"
#include "ITKDataConversionUtils.hxx"

namespace itk
{

  template <class TInputImage, class TOutputImage>
  SignalIntensityToS0ImageFilter<TInputImage, TOutputImage>::SignalIntensityToS0ImageFilter() : m_S0GradThresh(15.0f)
  { }

  template <class TInputImage, class TOutputImage>
  void SignalIntensityToS0ImageFilter<TInputImage, TOutputImage>
#if ITK_VERSION_MAJOR < 4
    ::ThreadedGenerateData(const typename Superclass::OutputImageRegionType & outputRegionForThread,
                           int itkNotUsed(threadId) )
#else
    ::ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread,
                           ThreadIdType itkNotUsed(threadId))
#endif
  {
    const InputImageType* inputVectorVolume = this->GetInput();
    OutputImageType* outputVolume = this->GetOutput();

    InputImageConstIterType inputVectorVolumeIter(inputVectorVolume, outputRegionForThread);
    OutputImageIterType outputVolumeIter(outputVolume, outputRegionForThread);

    while (!inputVectorVolumeIter.IsAtEnd())
    {
      typedef typename InputImageType::InternalPixelType InPixelType;
      typedef typename InternalVectorVoxelType::ComponentType OutPixelType;
      InternalVectorVoxelType vectorVoxel = ITKUtils::convertVectorType<InPixelType, OutPixelType>(inputVectorVolumeIter.Get());
      const OutputPixelType s0Value = computeS0(inputVectorVolume->GetNumberOfComponentsPerPixel(), vectorVoxel.GetDataPointer());
      outputVolumeIter.Set(s0Value);

      ++outputVolumeIter;
      ++inputVectorVolumeIter;
    }
  }

  template <class TInputImage, class TOutput>
  typename SignalIntensityToS0ImageFilter<TInputImage, TOutput>::OutputPixelType
  SignalIntensityToS0ImageFilter<TInputImage, TOutput>::computeS0(int signalSize, const float* signal)
  {
    int batIndex;
    try {
      batIndex = m_batEstimator->getBATIndex(signalSize, signal);
    }
    catch (...) {
      return 0;
    }

    double s0 = 0;
    float* signalGradient = new float[signalSize];
    compute_gradient(signalSize, signal, signalGradient);

    int count = 0;
    double sum = 0;
    for (int i = 0; i < batIndex; i++)
    {
      sum += signal[i];
      if (signalGradient[i] < m_S0GradThresh) {
        s0 += signal[i];
        count++;
      }
    }
    if (batIndex > 0) {
      if (count) {
        s0 /= count;
      }
      else {
        s0 = sum / (batIndex);
      }
    }
    else {
      s0 = signal[0];
    }

    delete[] signalGradient;
    return static_cast<OutputPixelType>(s0);
  }


  /** Standard "PrintSelf" method */
  template <class TInputImage, class TOutput>
  void SignalIntensityToS0ImageFilter<TInputImage, TOutput>::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);
    os << indent << "S0GradThresh: " << m_S0GradThresh << std::endl;
  }

} // end namespace itk

#endif
