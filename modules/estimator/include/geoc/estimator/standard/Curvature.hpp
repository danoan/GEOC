#include "geoc/estimator/standard/Curvature.h"

using namespace GEOC::Estimator::Standard;

template<typename IteratorType>
MDCACurvature::MDCACurvature(IteratorType itb,
                             IteratorType ite,
                             std::vector<double>& estimations)
{
    SegmentComputer sc;
    SCEstimator f;

    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDCACurvatureEstimator(sc,f);

    MCMDCACurvatureEstimator.init(1.0,itb,ite);
    MCMDCACurvatureEstimator.eval(itb,ite,std::back_inserter(estimations));
})

