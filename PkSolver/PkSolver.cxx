/*=auto=========================================================================

  Portions (c) Copyright 2009 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: pk_solver.cxx,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.13 $

  =========================================================================auto=*/

#include <vnl/algo/vnl_convolve.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include "PkSolver.h"
#include "SignalComputationUtils.h"
#include "itkTimeProbesCollectorBase.h"
#include <string>

namespace itk
{
  using namespace SignalUtils;
  //
  // Support routines/classes used internally in the PkSolver
  //
  static itk::TimeProbesCollectorBase probe;

  // TODO This is a very bad itermediate hack during step-wise refactoring. Will refactor PkSolver into a proper class and make this a member. 
  const BolusArrivalTime::BolusArrivalTimeEstimator* m_batEstimator;

  //
  // Implementation of the PkSolver API
  //
  //

  unsigned pk_solver(int signalSize, const float* timeAxis,
    const float* PixelConcentrationCurve,
    const float* BloodConcentrationCurve,
    float& Ktrans, float& Ve, float& Fpv,
    float fTol, float gTol, float xTol,
    float epsilon, int maxIter,
    float hematocrit,
    itk::LevenbergMarquardtOptimizer* optimizer,
    LMCostFunction* costFunction,
    int modelType,
    const BolusArrivalTime::BolusArrivalTimeEstimator* batEstimator)
  {
    //std::cout << "in pk solver" << std::endl;
    // probe.Start("pk_solver");
    // Note the unit: timeAxis should be in minutes!! This could be related to the following parameters!!
    // fTol      =  1e-4;  // Function value tolerance
    // gTol      =  1e-4;  // Gradient magnitude tolerance
    // xTol      =  1e-5;  // Search space tolerance
    // epsilon   =  1e-9;    // Step
    // maxIter   =   200;  // Maximum number of iterations
    //std::cerr << "In pkSolver!" << std::endl;

    m_batEstimator = batEstimator;

    // Levenberg Marquardt optimizer

    //////////////
    LMCostFunction::ParametersType initialValue;
    if (modelType == itk::LMCostFunction::TOFTS_2_PARAMETER)
    {
      initialValue = LMCostFunction::ParametersType(2); ///...
    }
    else
    {
      initialValue = LMCostFunction::ParametersType(3);
      initialValue[2] = 0.1;     //f_pv //...
    }
    initialValue[0] = 0.1;     //Ktrans //...
    initialValue[1] = 0.5;     //ve //...

    costFunction->SetNumberOfValues(signalSize);


    costFunction->SetCb(BloodConcentrationCurve, signalSize); //BloodConcentrationCurve
    costFunction->SetCv(PixelConcentrationCurve, signalSize); //Signal Y
    costFunction->SetTime(timeAxis, signalSize); //Signal X
    costFunction->SetHematocrit(hematocrit);
    costFunction->GetValue(initialValue); //...
    costFunction->SetModelType(modelType);

    optimizer->UseCostFunctionGradientOff();

    try
    {
      optimizer->SetCostFunction(costFunction);
    }
    catch (itk::ExceptionObject & e)
    {
      std::cout << "Exception thrown ! " << std::endl;
      std::cout << "An error ocurred during Optimization" << std::endl;
      std::cout << e << std::endl;
      return false;
    }

    itk::LevenbergMarquardtOptimizer::InternalOptimizerType * vnlOptimizer = optimizer->GetOptimizer();//...

    vnlOptimizer->set_f_tolerance(fTol); //...
    vnlOptimizer->set_g_tolerance(gTol); //...
    vnlOptimizer->set_x_tolerance(xTol); //...
    vnlOptimizer->set_epsilon_function(epsilon); //...
    vnlOptimizer->set_max_function_evals(maxIter); //...

    // We start not so far from the solution

    optimizer->SetInitialPosition(initialValue); //...

    try
    {
      //  probe.Start("optimizer");
      optimizer->StartOptimization();
      //   probe.Stop("optimizer");
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception thrown ! " << std::endl;
      std::cerr << "An error ocurred during Optimization" << std::endl;
      std::cerr << "Location    = " << e.GetLocation() << std::endl;
      std::cerr << "Description = " << e.GetDescription() << std::endl;
      return false;
    }
    //vnlOptimizer->diagnose_outcome();
    //std::cerr << "after optimizer!" << std::endl;
    itk::LevenbergMarquardtOptimizer::ParametersType finalPosition;
    finalPosition = optimizer->GetCurrentPosition();
    //std::cerr << finalPosition[0] << ", " << finalPosition[1] << ", " << finalPosition[2] << std::endl;


    //Solution: remove the scale of 100
    Ktrans = finalPosition[0];
    Ve = finalPosition[1];
    if (modelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
    {
      Fpv = finalPosition[2];
    }


    std::string diagnosticsCode = optimizer->GetStopConditionDescription();
    unsigned errorCode;
    for (errorCode = 0; errorCode < NumOptimizerDiagnosticCodes; errorCode++)
    {
      if (diagnosticsCode.find(OptimizerDiagnosticStrings[errorCode]) != std::string::npos)
        break;
    }

    // "Project" back onto the feasible set.  Should really be done as a
    // constraint in the optimization.
    if (Ve < 0)
    {
      Ve = 0;
      errorCode |= VE_CLAMPED;
    }
    if (Ve > 1)
    {
      Ve = 1;
      errorCode |= VE_CLAMPED;
    }
    if (Ktrans < 0)
    {
      Ktrans = 0;
      errorCode |= KTRANS_CLAMPED;
    }
    if (Ktrans > 5)
    {
      Ktrans = 5;
      errorCode |= KTRANS_CLAMPED;
    }

    //if((Fpv>1)||(Fpv<0)) Fpv = 0;
    //  probe.Stop("pk_solver");
    return errorCode;
  }

  void pk_report()
  {
    probe.Report();
  }

  void pk_clear()
  {
    probe.Clear();
  }

  float area_under_curve(int signalSize,
    const float* timeAxis,
    const float* concentration,
    int BATIndex,
    float aucTimeInterval)
  {
    float auc = 0.0f;
    if (BATIndex >= signalSize) return auc;
    //std::cerr << std::endl << "BATIndex:"<<BATIndex << std::endl;

    int lastIndex = BATIndex;
    //std::cerr << std::endl << "timeAxis[0]:"<< timeAxis[0] << std::endl;
    float targetTime = timeAxis[BATIndex] + aucTimeInterval;
    //std::cerr << std::endl << "targetTime:"<< targetTime << std::endl;
    float tempTime = timeAxis[BATIndex + 1];
    //std::cerr << std::endl << "tempTime:"<< tempTime << std::endl;

    //find the last index
    while ((tempTime < targetTime) && (lastIndex < (signalSize - 2)))
    {
      lastIndex += 1;
      tempTime = timeAxis[lastIndex + 1];
      //std::cerr << std::endl << "tempTime"<<tempTime << std::endl;
    }

    if ((lastIndex - BATIndex) == 0) return auc = aucTimeInterval*concentration[BATIndex];

    //extract time and concentration
    float * concentrationValues = new float[lastIndex - BATIndex + 2]();
    float * timeValues = new float[lastIndex - BATIndex + 2]();

    //find the extra time and concentration value for auc
    float y1, y2, x1, x2, slope, b, targetX, targetY;
    y2 = concentration[lastIndex + 1];
    y1 = concentration[lastIndex];
    x2 = timeAxis[lastIndex + 1];
    x1 = timeAxis[lastIndex];
    slope = (y2 - y1) / (x2 - x1);
    b = y1 - slope*x1;
    targetX = timeAxis[BATIndex] + aucTimeInterval;
    targetY = slope*targetX + b;
    if (targetX > timeAxis[signalSize - 1])
    {
      targetX = timeAxis[lastIndex + 1];
      targetY = concentration[lastIndex + 1];
    }
    concentrationValues[lastIndex - BATIndex + 1] = targetY; //put the extra value at the end
    timeValues[lastIndex - BATIndex + 1] = targetX; //put the extra time value at the end

    //	printf("lastIndex is %f\n", (float)lastIndex);
    for (int i = 0; i < (lastIndex - BATIndex + 1); ++i)
    {
      concentrationValues[i] = concentration[i + BATIndex];
      timeValues[i] = timeAxis[i + BATIndex];
      //printf("lastIndex is %f,%f\n", (float) concentrationValues[i],(float)timeValues[i]);
    }

    //get auc
    auc = intergrate(concentrationValues, timeValues, (lastIndex - BATIndex + 2));

    delete[] concentrationValues;
    delete[] timeValues;
    return auc;
  }

  float intergrate(float* yValues, float * xValues, int size)
  {
    float area = 0.0f;
    for (int i = 1; i < size; ++i)
    {
      area += (xValues[i] - xValues[i - 1])*(yValues[i] + yValues[i - 1]) / 2;
      //	std::cerr << std::endl << "area:" << area<<","<<yValues[i] << std::endl;
    }
    return area;
  }

  void compute_derivative(int signalSize,
    const float* SignalY,
    float* YDeriv)
  {
    YDeriv[0] = (float)((-3.0*SignalY[0] + 4.0*SignalY[1] - SignalY[2]) / 2.0);
    YDeriv[signalSize - 1] = (float)((3.0*SignalY[signalSize - 1] - 4.0*SignalY[signalSize - 2] + SignalY[signalSize - 3]) / 2.0);
    for (int i = 1; i < signalSize - 1; i++)
    {
      YDeriv[i] = (float)((SignalY[i + 1] - SignalY[i - 1]) / 2.0);
    }
  }

  void compute_derivative_forward(int signalSize,
    const float* SignalY,
    float* YDeriv)
  {
    YDeriv[signalSize - 1] = (float)(SignalY[signalSize - 1] - SignalY[signalSize - 2]);
    for (int i = 0; i < signalSize - 1; i++)
    {
      YDeriv[i] = (float)(SignalY[i + 1] - SignalY[i]);
    }
  }

  void compute_derivative_backward(int signalSize,
    const float* SignalY,
    float* YDeriv)
  {
    YDeriv[0] = (float)(-SignalY[0] + SignalY[1]);
    for (int i = 1; i < signalSize; i++)
    {
      YDeriv[i] = (float)(SignalY[i] - SignalY[i - 1]);
    }
  }

  void compute_gradient_old(int signalSize, const float* SignalY, float* SignalGradient)
  {
    typedef itk::Image<float, 1>   ImageType;
    typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType>  FilterType;

    ImageType::Pointer SignalImg1D = ImageType::New();
    ImageType::SizeType imgSize;
    imgSize[0] = signalSize;
    SignalImg1D->SetRegions(imgSize);
    SignalImg1D->Allocate();

    typedef itk::ImageRegionIterator< ImageType > IteratorType;
    IteratorType inputIt(SignalImg1D, SignalImg1D->GetLargestPossibleRegion());
    inputIt.GoToBegin();
    for (unsigned int i = 0; !inputIt.IsAtEnd(); i++, ++inputIt)
    {
      inputIt.Set(SignalY[i]);
    }

    FilterType::Pointer gfilter = FilterType::New();
    gfilter->SetInput(SignalImg1D);
    gfilter->Update();

    ImageType::Pointer GradientImg1D = gfilter->GetOutput();

    IteratorType outputIt(GradientImg1D, GradientImg1D->GetLargestPossibleRegion());
    outputIt.GoToBegin();
    for (unsigned int i = 0; !outputIt.IsAtEnd(); i++, ++outputIt) {
      SignalGradient[i] = outputIt.Get();
    }
  }

  void compute_gradient(int signalSize, const float* SignalY, float* SignalGradient)
  {
    compute_derivative(signalSize, SignalY, SignalGradient);
    for (int i = 0; i < signalSize; i++)
    {
      SignalGradient[i] = sqrt(SignalGradient[i] * SignalGradient[i]);
    }
  }

  void compute_gradient_forward(int signalSize, const float* SignalY, float* SignalGradient)
  {
    compute_derivative_forward(signalSize, SignalY, SignalGradient);
    for (int i = 0; i < signalSize; i++)
    {
      SignalGradient[i] = sqrt(SignalGradient[i] * SignalGradient[i]);
    }
  }

  void compute_gradient_backward(int signalSize, const float* SignalY, float* SignalGradient)
  {
    compute_derivative_backward(signalSize, SignalY, SignalGradient);
    for (int i = 0; i < signalSize; i++)
    {
      SignalGradient[i] = sqrt(SignalGradient[i] * SignalGradient[i]);
    }
  }

}; // end of namespace
