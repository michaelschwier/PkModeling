#ifndef __ITKDataConversionUtils_hxx
#define __ITKDataConversionUtils_hxx

#include "itkVectorImage.h"

namespace ITKUtils {

  template<class TInputVoxelType, class TReturnVoxelType>
  itk::VariableLengthVector<TReturnVoxelType> convertVectorType(const itk::VariableLengthVector<TInputVoxelType>& inputVector)
  {
    itk::VariableLengthVector<TReturnVoxelType> outVector;
    outVector.SetSize(inputVector.GetSize());
    outVector.Fill(0);
    outVector += inputVector; // shorthand for a copy/cast
    return outVector;
  }

}

#endif
