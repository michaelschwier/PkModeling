/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConcentrationToQuantitativeImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2012/05/01 14:28:51 $
  Version:   $Revision: 0.0 $
  =========================================================================*/
#ifndef __itkConcentrationToQuantitativeImageFilter_h
#define __itkConcentrationToQuantitativeImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"
#include "PkSolver.h"
#include <string>
#include "AIF/ArterialInputFunction.h"

namespace itk
{
  /** \class ConcentrationToQuantitativeImageFilter
   * \brief Calculates quantitative imaging parameters from concentration curves.
   *
   * This filter computes Pk modeling quantitative images from
   * concentration curves. The input volume is a vector image of
   * concentration curves represented in floating point.  The output is
   * a series of floating point images of quantitative parameters.
   *
   * An second input, specifying the location of the arterial input
   * function, allows for the calculation to be adjusted for blood
   * verses tissue.
   *
   * \note
   * This work is part of the National Alliance for Medical Image Computing
   * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
   * for Medical Research, Grant U54 EB005149.
   *
   */

  template <class TInputImage, class TMaskImage, class TOutputImage>
  class ITK_EXPORT ConcentrationToQuantitativeImageFilter : public ImageToImageFilter < TInputImage, TOutputImage >
  {
  public:
    typedef ConcentrationToQuantitativeImageFilter Self;

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage                             VectorVolumeType;
    typedef typename VectorVolumeType::Pointer      VectorVolumePointerType;
    typedef typename VectorVolumeType::ConstPointer VectorVolumeConstPointerType;
    typedef typename VectorVolumeType::PixelType    VectorVolumePixelType;
    typedef typename VectorVolumeType::RegionType   VectorVolumeRegionType;
    typedef typename VectorVolumeType::SizeType     VectorVolumeSizeType;
    typedef itk::ImageRegionConstIterator<VectorVolumeType> VectorVolumeConstIterType;
    typedef itk::ImageRegionIterator<VectorVolumeType> VectorVolumeIterType;

    typedef TMaskImage                            MaskVolumeType;
    typedef typename MaskVolumeType::Pointer      MaskVolumePointerType;
    typedef typename MaskVolumeType::ConstPointer MaskVolumeConstPointerType;
    typedef typename MaskVolumeType::PixelType    MaskVolumePixelType;
    typedef typename MaskVolumeType::RegionType   MaskVolumeRegionType;
    typedef typename MaskVolumeType::SizeType     MaskVolumeSizeType;
    typedef itk::ImageRegionConstIterator<MaskVolumeType> MaskVolumeConstIterType;

    typedef TOutputImage                                    OutputVolumeType;
    typedef typename OutputVolumeType::Pointer              OutputVolumePointerType;
    typedef typename OutputVolumeType::ConstPointer         OutputVolumeConstPointerType;
    typedef typename OutputVolumeType::PixelType            OutputVolumePixelType;
    typedef typename OutputVolumeType::RegionType           OutputVolumeRegionType;
    typedef typename OutputVolumeType::IndexType            OutputVolumeIndexType;
    typedef itk::ImageRegionIterator<OutputVolumeType>      OutputVolumeIterType;
    typedef itk::ImageRegionConstIterator<OutputVolumeType> OutputVolumeConstIterType;

    /** Methods to implement smart pointers and work with the itk object factory
     **/
    typedef ProcessObject              ProcessObjectSuperclass;
    typedef SmartPointer< Self >       Pointer;
    //typedef SmartPointer< const Self > ConstPointer;
    typedef itk::DataObject::Pointer   DataObjectPointer;

    typedef ProcessObject::DataObjectPointerArraySizeType DataObjectPointerArraySizeType;
    using ProcessObjectSuperclass::MakeOutput;
    virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType idx);

    typedef itk::VariableLengthVector<float>          VectorVoxelType;

    /** Standard class typedefs. */
    typedef ImageToImageFilter<VectorVolumeType, OutputVolumeType> Superclass;
    typedef SmartPointer<const Self>                  ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(ConcentrationToQuantitativeImageFilter, ImageToImageFilter);

    /** ImageDimension enumeration */
    itkStaticConstMacro(VectorVolumeDimension, unsigned int,
      VectorVolumeType::ImageDimension);
    itkStaticConstMacro(MaskVolumeDimension, unsigned int,
      MaskVolumeType::ImageDimension);
    itkStaticConstMacro(OutputVolumeDimension, unsigned int,
      OutputVolumeType::ImageDimension);

    /** Set and get the parameters to control the calculation of
    quantified valued */
    itkGetMacro(T1Pre, float);
    itkSetMacro(T1Pre, float);
    itkGetMacro(TR, float);
    itkSetMacro(TR, float);
    itkGetMacro(FA, float);
    itkSetMacro(FA, float);
    itkGetMacro(RGD_relaxivity, float);
    itkSetMacro(RGD_relaxivity, float);
    itkGetMacro(S0GradThresh, float);
    itkSetMacro(S0GradThresh, float);
    itkGetMacro(fTol, float);
    itkSetMacro(fTol, float);
    itkGetMacro(gTol, float);
    itkSetMacro(gTol, float);
    itkGetMacro(xTol, float);
    itkSetMacro(xTol, float);
    itkGetMacro(epsilon, float);
    itkSetMacro(epsilon, float);
    itkGetMacro(maxIter, int);
    itkSetMacro(maxIter, int);
    itkGetMacro(hematocrit, float);
    itkSetMacro(hematocrit, float);
    itkGetMacro(AUCTimeInterval, float);
    itkSetMacro(AUCTimeInterval, float);
    itkGetMacro(ModelType, int);
    itkSetMacro(ModelType, int);

    void SetBatEstimator(const BolusArrivalTime::BolusArrivalTimeEstimator* batEstimator)
    {
      m_batEstimator = batEstimator;
    }

    const BolusArrivalTime::BolusArrivalTimeEstimator* GetBatEstimator() const
    {
      return this->m_batEstimator;
    }

    void SetTiming(const std::vector<float>& inputTiming);
    const std::vector<float>& GetTiming();

    /// Set a mask to specify where the model fit is be calculated from the
    /// input concentration image.
    void SetROIMask(const MaskVolumeType* volume);

    /// Get the mask that specifies from where the model fit is calculated
    const TMaskImage* GetROIMask() const;


    /// Set the AIF
    void SetAIF(const ArterialInputFunction* aif);

    /// Get the quantitative output images
    TOutputImage* GetKTransOutput();
    TOutputImage* GetVEOutput();
    TOutputImage* GetFPVOutput();
    TOutputImage* GetMaxSlopeOutput();
    TOutputImage* GetAUCOutput();
    TOutputImage* GetRSquaredOutput();
    TOutputImage* GetBATOutput();

    VectorVolumeType* GetFittedDataOutput();

    TOutputImage* GetOptimizerDiagnosticsOutput();

  protected:
    ConcentrationToQuantitativeImageFilter();
    ~ConcentrationToQuantitativeImageFilter(){
    }
    void PrintSelf(std::ostream& os, Indent indent) const;

    void BeforeThreadedGenerateData();

#if ITK_VERSION_MAJOR < 4
    void ThreadedGenerateData( const OutputVolumeRegionType& outputRegionForThread, int threadId );

#else
    void ThreadedGenerateData(const OutputVolumeRegionType& outputRegionForThread,
      ThreadIdType threadId);

#endif


  private:
    ConcentrationToQuantitativeImageFilter(const Self &); // purposely not implemented
    void operator=(const Self &); // purposely not implemented

    float  m_T1Pre;
    float  m_TR;
    float  m_FA;
    float  m_RGD_relaxivity;
    float  m_S0GradThresh;
    float  m_fTol;
    float  m_gTol;
    float  m_xTol;
    float  m_epsilon;
    int    m_maxIter;
    float  m_hematocrit;
    float  m_AUCTimeInterval;
    int    m_AIFBATIndex;
    int    m_ModelType;
    const BolusArrivalTime::BolusArrivalTimeEstimator* m_batEstimator;
    const ArterialInputFunction* m_aif;

    std::vector<float> m_Timing;

    // variables to cache information to share between threads
    std::vector<float> m_AIF;
    float  m_aifAUC;
  };

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConcentrationToQuantitativeImageFilter.hxx"
#endif

#endif
