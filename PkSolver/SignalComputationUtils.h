#ifndef __SignalComputationUtils_h
#define __SignalComputationUtils_h

#include <vector>
#include <math.h>

// work around compile error on Windows.
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif


namespace SignalUtils {
  int getMaxPosition(int signalSize, const float* signal);
  int getMaxPositionInRange(int start, int stop, const float* signal);
  std::vector<float> resampleSignal(std::vector<float> signalTime, std::vector<float> signal, std::vector<float> referenceTime);
}

#endif
