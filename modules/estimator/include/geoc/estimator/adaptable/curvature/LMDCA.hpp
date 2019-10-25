#ifndef GEOC_ESTIMATOR_CURVATURE_LMDCA_H
#define GEOC_ESTIMATOR_CURVATURE_LMDCA_H

#include <DGtal/geometry/curves/StabbingCircleComputer.h>
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include <DGtal/geometry/curves/estimation/FunctorsLambdaMST.h>
#include <geoc/estimator/adaptable/base/LMDCAEstimator.h>

namespace GEOC
{
    namespace Estimator
    {
        namespace Alternative
        {
            template<typename IteratorType>
            struct LMDCACurvature
            {
                typedef DGtal::StabbingCircleComputer<IteratorType> SegmentComputer;
                typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;
                typedef DGtal::functors::Lambda64Function LambdaFunction;

                LMDCACurvature(IteratorType itb,
                                   IteratorType ite,
                                   std::vector<double>& estimations,
                               double h,
                               void* data)
                {
                    SegmentComputer sc;
                    SCEstimator f;

                    LMDCAEstimator<SegmentComputer,SCEstimator, LambdaFunction> LMDCACurvatureEstimator(sc,f);

                    LMDCACurvatureEstimator.init(h,itb,ite);
                    LMDCACurvatureEstimator.eval(itb,ite,std::back_inserter(estimations));
                }
            };
        }
    }
}
#endif //GEOC_ESTIMATOR_CURVATURE_LMDCA_H
