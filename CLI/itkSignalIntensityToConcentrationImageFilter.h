/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkSignalIntensityToConcentrationImageFilter.h $
  Language:  C++
  Date:      $Date: 2012/03/07 $
  Version:   $Revision: 0.0 $

  =========================================================================*/
#ifndef __itkSignalIntensityToConcentrationImageFilter_h
#define __itkSignalIntensityToConcentrationImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkSignalIntensityToS0ImageFilter.h"
#include "itkImageFileWriter.h"


namespace itk
{
  /** \class SignalIntensityToConcentrationImageFilter
   * \brief Convert from signal intensities to concentrations.
   *
   * This converts an VectorImage of signal intensities into a
   * VectorImage of concentration values. Typical use if for the output
   * image's pixel component type to be floating point.
   *
   * An second input, specifying the location of the arterial input
   * function, allows for the calculation to be adjusted for blood
   * versus tissue.
   *
   * \note
   * This work is part of the National Alliance for Medical Image Computing
   * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
   * for Medical Research, Grant U54 EB005149.
   *
   */
  template <class TInputImage, class TMaskImage, class TOutputImage>
  class SignalIntensityToConcentrationImageFilter : public ImageToImageFilter < TInputImage, TOutputImage >
  {
  public:
    /** Standard class typedefs. */
    typedef SignalIntensityToConcentrationImageFilter           Self;
    typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
    typedef SmartPointer<Self>                                  Pointer;
    typedef SmartPointer<const Self>                            ConstPointer;

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage                                   InputImageType;
    typedef typename InputImageType::Pointer              InputImagePointerType;
    typedef typename InputImageType::ConstPointer         InputImageConstPointerType;
    typedef typename InputImageType::PixelType            InputPixelType;
    typedef itk::ImageRegionConstIterator<InputImageType> InputImageConstIterType;

    typedef TMaskImage                                   InputMaskType;
    typedef itk::ImageRegionConstIterator<InputMaskType> InputMaskConstIterType;

    typedef TOutputImage                              OutputImageType;
    typedef typename OutputImageType::Pointer         OutputImagePointerType;
    typedef typename OutputImageType::ConstPointer    OutputImageConstPointerType;
    typedef typename OutputImageType::PixelType       OutputPixelType;
    typedef typename OutputImageType::RegionType      OutputImageRegionType;
    typedef itk::ImageRegionIterator<OutputImageType> OutputIterType;

  private:
    typedef itk::Image<float, TInputImage::ImageDimension> InternalVolumeType;
    typedef typename InternalVolumeType::Pointer           InternalVolumePointerType;
    typedef itk::ImageRegionIterator<InternalVolumeType>   InternalVolumeIterType;
    typedef typename InternalVolumeType::RegionType        InternalVolumeRegionType;
    typedef typename InternalVolumeType::SizeType          InternalVolumeSizeType;

    typedef itk::VectorImage<float, TInputImage::ImageDimension>    InternalVectorVolumeType;
    typedef typename InternalVectorVolumeType::Pointer              InternalVectorVolumePointerType;
    typedef itk::ImageRegionConstIterator<InternalVectorVolumeType> InternalVectorVolumeConstIterType;
    typedef typename InternalVectorVolumeType::RegionType           InternalVectorVolumeRegionType;
    typedef typename InternalVectorVolumeType::SizeType             InternalVectorVolumeSizeType;

    typedef itk::VariableLengthVector<float> InternalVectorVoxelType;

  public:
    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SignalIntensityToConcentrationImageFilter, ImageToImageFilter);

    /** Parameter getters/setters. */
    itkGetMacro(T1PreBlood, float);
    itkSetMacro(T1PreBlood, float);
    itkGetMacro(T1PreTissue, float);
    itkSetMacro(T1PreTissue, float);
    itkGetMacro(TR, float);
    itkSetMacro(TR, float);
    itkGetMacro(FA, float);
    itkSetMacro(FA, float);
    itkGetMacro(RGD_relaxivity, float);
    itkSetMacro(RGD_relaxivity, float);
    itkGetMacro(S0GradThresh, float);
    itkSetMacro(S0GradThresh, float);

    void SetBatEstimator(const BolusArrivalTime::BolusArrivalTimeEstimator* batEstimator)
    {
      m_batEstimator = batEstimator;
      this->Modified();
    }

    // Set a mask image for specifying the location of the arterial
    // input function. This is interpretted as a binary image with
    // nonzero values only at the arterial input function locations.
    void SetAIFMask(InputMaskType* aifMaskVolume)
    {
      this->SetNthInput(1, const_cast<InputMaskType*>(aifMaskVolume));
    }

    // Get the mask image assigned as the arterial input function
    const InputMaskType* GetAIFMask() const
    {
      return dynamic_cast<const InputMaskType*>(this->ProcessObject::GetInput(1));
    }

    // Set a mask image for specifying the location of voxels for model fit.
    void SetROIMask(InputMaskType* roiMaskVolume)
    {
      this->SetNthInput(2, const_cast<InputMaskType*>(roiMaskVolume));
    }

    // Get the mask image specifying the location of voxels for model fit.
    const InputMaskType* GetROIMask() const
    {
      return dynamic_cast<const InputMaskType*>(this->ProcessObject::GetInput(2));
    }

    // Set a T1 Map image for T1
    void SetT1Map(InputMaskType* T1MapVolume)
    {
      this->SetNthInput(3, const_cast<InputMaskType*>(T1MapVolume));
    }

    // Get the mask image specifying the location of voxels for model fit.
    const InputMaskType* GetT1Map() const
    {
      return dynamic_cast<const InputMaskType*>(this->ProcessObject::GetInput(3));
    }

  protected:
    SignalIntensityToConcentrationImageFilter();

    virtual ~SignalIntensityToConcentrationImageFilter()
    { }

    void GenerateData();
    OutputImageType* GetAllocatedOutputVolume(const InputImageType* inputVectorVolume);
    InternalVolumePointerType GetS0Image(const InputImageType* inputVectorVolume);
    void convertSignalToConcentration(unsigned int signalSize,
                                      const float* SignalIntensityCurve,
                                      const float T1Pre, float TR, float FA,
                                      float* concentration,
                                      float RGd_relaxivity,
                                      float s0,
                                      float S0GradThresh);

    void PrintSelf(std::ostream& os, Indent indent) const;

  private:
    SignalIntensityToConcentrationImageFilter(const Self &); // not implemented on purpose
    void operator=(const Self &); // not implemented on purpose

    float m_T1PreBlood;
    float m_T1PreTissue;
    float m_TR;
    float m_FA;
    float m_RGD_relaxivity;
    float m_S0GradThresh;
    const BolusArrivalTime::BolusArrivalTimeEstimator* m_batEstimator;

    //------------------------------------------------------------------------------------------------------------------------
    //! Private internal helper class to handle getting the correct T1Pre value.
    //! Use it like an Iterator to walk through the voxel positions.
    class T1PreValueMapper {
    public:
      //! Instantiate the Mapper by providing ROI mask, AIF mask, and/or T1 Map (all of which are optional and my be NULL if not available).
      //! Also provide default constant Tissue and Blood value (these are required inputs).
      T1PreValueMapper(const InputMaskType* roiMask, const InputMaskType* aifMask, const InputMaskType* t1Map, float t1PreTissue, float t1PreBlood);
      virtual ~T1PreValueMapper();

      //! Returns the T1Pre value for the current voxel position, based on the availability and validity of ROI/AIF mask and T1 Map at this position.
      float Get();
      void GoToBegin();
      T1PreValueMapper& operator++();

    private:
      InputMaskConstIterType* getNewConstMaskIterOrNull(const InputMaskType* inMask);

      InputMaskConstIterType* roiMaskVolumeIter;
      InputMaskConstIterType* aifMaskVolumeIter;
      InputMaskConstIterType* T1MapVolumeIter;
      float m_T1PreTissue;
      float m_T1PreBlood;

    }; // end T1PreValueIterator class
    //------------------------------------------------------------------------------------------------------------------------

  };

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSignalIntensityToConcentrationImageFilter.hxx"
#endif

#endif
