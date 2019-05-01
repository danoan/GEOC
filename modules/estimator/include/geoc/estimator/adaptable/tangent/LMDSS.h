#ifndef GEOC_ESTIMATORS_CURVATURE_LMDSS_H
#define GEOC_ESTIMATORS_CURVATURE_LMDSS_H

#include <DGtal/kernel/PointVector.h>
#include <DGtal/geometry/curves/ArithmeticalDSSComputer.h>
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>

namespace GEOC
{
    namespace Estimator
    {
        namespace Alternative
        {
            template<typename IteratorType>
            struct LMDSSTangent
            {
                typedef DGtal::ArithmeticalDSSComputer< IteratorType, int, 4 > SegmentComputer;
                typedef DGtal::TangentFromDSSEstimator<SegmentComputer> SCEstimator;

                typedef DGtal::PointVector<2,double> TangentVector;

                LMDSSTangent(IteratorType itb,
                             IteratorType ite,
                             std::vector< TangentVector >& estimations,
                             double h)
                {}


            };
        }
    }
}

#endif //GEOC_ESTIMATORS_CURVATURE_LMDSS_H
