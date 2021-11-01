#ifndef GEOC_ESTIMATOR_CURVATURE_MDCA_H
#define GEOC_ESTIMATOR_CURVATURE_MDCA_H

#include <DGtal/geometry/curves/StabbingCircleComputer.h>
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include <DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h>

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
                              std::vector<double>& estimations,
                              double h,
                              void* data)
                {
                    SegmentComputer sc;
                    SCEstimator f;

                    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDCACurvatureEstimator(sc,f);

                    MCMDCACurvatureEstimator.init(itb,ite);
                    MCMDCACurvatureEstimator.eval(itb,ite,std::back_inserter(estimations),h);
                }
            };
        }
    }
}

#endif //GEOC_ESTIMATOR_CURVATURE_MDCA_H
