#ifndef GEOC_ESTIMATOR_STANDARD_CURVATURE_H
#define GEOC_ESTIMATOR_STANDARD_CURVATURE_H

#include "DGtal/geometry/curves/StabbingCircleComputer.h"
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"

namespace GEOC
{
    namespace Estimator
    {
        namespace Standard
        {
            template<typename IteratorType>
            struct MDCACurvature
            {
                typedef DGtal::StabbingCircleComputer<IteratorType> SegmentComputer;
                typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;

                MDCACurvature(IteratorType itb,
                              IteratorType ite,
                              std::vector<double>& estimations);
            };
        }
    }
}

#include "Curvature.hpp"

#endif //GEOC_ESTIMATOR_STANDARD_CURVATURE_H
