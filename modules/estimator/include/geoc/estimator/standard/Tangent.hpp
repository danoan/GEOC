#include "geoc/estimator/standard/Tangent.h"

using namespace GEOC::Estimator::Standard;

template<typename IteratorType>
MDSSTangent<IteratorType>::MDSSTangent(IteratorType itb,
                                       IteratorType ite,
                                       std::vector< TangentVector >& estimations,
                                       double h)
{
    SegmentComputer sc;
    SCEstimator f;

    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDSSTangentEstimator(sc,f);

    MCMDSSTangentEstimator.init(h,itb,ite);
    MCMDSSTangentEstimator.eval(itb,ite,std::back_inserter(estimations));
}