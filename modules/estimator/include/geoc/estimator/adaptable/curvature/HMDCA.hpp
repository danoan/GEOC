#ifndef GEOC_ESTIMATOR_CURVATURE_HMDCA_H
#define GEOC_ESTIMATOR_CURVATURE_HMDCA_H

#include <DGtal/geometry/curves/StabbingCircleComputer.h>
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include <geoc/estimator/adaptable/base/HMDCAEstimator.h>

namespace GEOC
{
    namespace Estimator
    {
        namespace Alternative
        {
            //Highest value DCA estimator
            template<typename IteratorType>
            struct HMDCACurvature
            {
                typedef DGtal::StabbingCircleComputer<IteratorType> SegmentComputer;
                typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;

                HMDCACurvature(IteratorType itb,
                              IteratorType ite,
                              std::vector<double>& estimations,
                               double h,
                               void* data)
                {
                    SegmentComputer sc;
                    SCEstimator f;

                    HMDCAEstimator<SegmentComputer,SCEstimator> PessimistMDCACurvature(sc,f);

                    PessimistMDCACurvature.init(h,itb,ite);

                    PessimistMDCACurvature.eval(itb,ite,std::back_inserter(estimations));
                }
            };
        }
    }
}
#endif //GEOC_ESTIMATOR_CURVATURE_HMDCA_H
