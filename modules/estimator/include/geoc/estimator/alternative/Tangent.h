#ifndef GEOC_ESTIMATOR_TANGENT_H
#define GEOC_ESTIMATOR_TANGENT_H

#include <DGtal/helpers/StdDefs.h>

#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"

namespace GEOC
{
    namespace Estimator
    {
        namespace Alternative
        {
            template<typename IteratorType>
            struct MDSSLambdaTangent
            {
                typedef DGtal::ArithmeticalDSSComputer< IteratorType, int, 4 > SegmentComputer;
                typedef DGtal::TangentFromDSSEstimator<SegmentComputer> SCEstimator;

                typedef DGtal::PointVector<2,double> TangentVector;

                MDSSLambdaTangent(IteratorType itb,
                                  IteratorType ite,
                                  std::vector< TangentVector >& estimations);


            }
        }
    }
}

#include "Tangent.hpp"

#endif //GEOC_TANGENT_H
