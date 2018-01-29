#include "ArterialInputFunctionPopulation.h"

#include "SignalComputationUtils.h"
#include <math.h>


ArterialInputFunctionPopulation::ArterialInputFunctionPopulation(const std::vector<float>& referenceSignalTime,
                                                                 const float bolusArrivalTimeFraction)
                                                                 : m_referenceSignalTime(referenceSignalTime),
                                                                   m_bolusArrivalTimeFraction(bolusArrivalTimeFraction)
{
  m_aif = computeAIF();
}

std::vector<float> ArterialInputFunctionPopulation::getSignalValues() const
{
  return m_aif;
}

unsigned int ArterialInputFunctionPopulation::getSignalSize() const
{
  return m_aif.size();
}

std::vector<float> ArterialInputFunctionPopulation::computeAIF() const
{
  std::vector<float> AIF;

  // Make a high resolution timing vector as input to the AIF construction.
  std::vector<float> aif_time(m_referenceSignalTime.size() * 10);
  const float final_time_point = m_referenceSignalTime[m_referenceSignalTime.size() - 1];
  const float resolution = final_time_point / (aif_time.size() - 1);
  for (std::size_t j = 0; j < aif_time.size(); ++j)
  {
    aif_time[j] = resolution * j;
  }

  std::size_t bolus_arrival_time_idx = aif_time.size() * m_bolusArrivalTimeFraction;

  std::size_t n = aif_time.size();
  AIF.resize(n);

  std::size_t numTimePoints = n - bolus_arrival_time_idx;
  std::vector<float> timeSinceBolus(numTimePoints);


  // These time points "start" when the bolus arrives.
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    timeSinceBolus[j] = aif_time[bolus_arrival_time_idx + j] - aif_time[bolus_arrival_time_idx];
  }

  // Parker
  // defining parameters
  const double a1(0.809);
  const double a2(0.330);
  const double T1(0.17406);
  const double T2(0.365);
  const double sigma1(0.0563);
  const double sigma2(0.132);
  const double alpha(1.050);
  const double beta(0.1685);
  const double s(38.078);
  const double tau(0.483);


  // term0=alpha*exp(-beta*t)./(1+exp(-s*(t-tau)));
  // Here the assumption is that time is in minutes, so must scale accordingly.
  // see Parker.
  std::vector<double> term0(numTimePoints);
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    term0[j] = alpha * exp(-beta*timeSinceBolus[j] / 60.0)
      / (1 + exp(-s * (timeSinceBolus[j] / 60.0 - tau)));
  }

  const double A1 = a1 / (sigma1 * pow((2 * M_PI), 0.5));

  // B1=exp(-(t-T1).^2./(2.*sigma1^2));
  double numerator, denominator;
  std::vector<double> B1(numTimePoints);
  denominator = 2.0 * pow(sigma1, 2.0);
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    numerator = -1 * pow(-(timeSinceBolus[j] / 60.0 - T1), 2.0);
    B1[j] = exp(numerator / denominator);
  }

  // term1=A1.*B1;
  std::vector<double> term1(numTimePoints);
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    term1[j] = A1 * B1[j];
  }

  // A2=a2/(sigma2*((2*pi)^0.5));
  const double A2 = a2 / (sigma2 * pow(2 * M_PI, 0.5));

  //B2=exp(-(t-T2).^2./(2.*sigma2^2));
  std::vector<double> B2(numTimePoints);
  denominator = 2.0 * pow(sigma2, 2.0);
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    numerator = -1 * pow(-(timeSinceBolus[j] / 60.0 - T2), 2.0);
    B2[j] = exp(numerator / denominator);
  }

  // term2=A2.*B2;
  std::vector<double> term2(numTimePoints);
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    term2[j] = A2 * B2[j];
  }

  // aifPost=term0+term1+term2;
  std::vector<double> aifPost(numTimePoints);
  for (std::size_t j = 0; j < numTimePoints; ++j)
  {
    aifPost[j] = term0[j] + term1[j] + term2[j];
  }

  // Initialize values before bolus arrival.
  for (std::size_t j = 0; j < bolus_arrival_time_idx; ++j)
  {
    AIF[j] = 0;
  }

  // Shift the data to take into account the bolus arrival time.
  for (std::size_t j = bolus_arrival_time_idx; j < AIF.size(); ++j)
  {
    AIF[j] = aifPost[j - bolus_arrival_time_idx];
  }

  // Resample back to signal (sequence) time.
  std::vector<float> rAIF = SignalUtils::resampleSignal(aif_time, AIF, m_referenceSignalTime);

  return rAIF;
}
