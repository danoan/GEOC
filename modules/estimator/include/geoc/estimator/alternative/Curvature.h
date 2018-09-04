#ifndef GEOC_ESTIMATOR_ALTERNATIVE_CURVATURE_H
#define GEOC_ESTIMATOR_ALTERNATIVE_CURVATURE_H

#include "DGtal/geometry/curves/StabbingCircleComputer.h"
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include <geoc/estimator/base/DCALambdaEstimator.h>
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"

#include "DGtal/geometry/curves/estimation/LambdaMST2D.h"

namespace GEOC
{
    namespace Estimator
    {
        namespace Alternative
        {
            template<typename IteratorType>
            struct DCALambdaCurvature
            {
                typedef DGtal::StabbingCircleComputer<IteratorType> SegmentComputer;
                typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;
                typedef DGtal::functors::Lambda64Function LambdaFunction;

                typedef GEOC::Estimator::Base::DCALambdaEstimator LambdaMSE;

                DCALambdaCurvature(IteratorType itb,
                                   IteratorType ite,
                                   std::vector<double>& estimations);
            };

            //Highest value DCA estimator
            template<typename IteratorType>
            struct HDCACurvature
            {
                typedef DGtal::StabbingCircleComputer<IteratorType> SegmentComputer;
                typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;

                HDCACurvature(IteratorType itb,
                              IteratorType ite,
                              std::vector<double>& estimations);
            };
        }
    }
}

#endif //GEOC_ESTIMATOR_ALTERNATIVE_
