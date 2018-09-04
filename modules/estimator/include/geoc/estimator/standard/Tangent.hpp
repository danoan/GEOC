#include "geoc/estimator/standard/Tangent.h"

using namespace GEOC::Estimator::Standard;

template<typename IteratorType>
MDSSTangent<IteratorType>::MDSSTangent(IteratorType itb,
                                       IteratorType ite,
                                       std::vector< TangentVector >& estimations)
{
    SegmentComputer sc;
    SCEstimator f;

    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDSSTangentEstimator(sc,f);

    MCMDSSTangentEstimator.init(1.0,itb,ite);
    MCMDSSTangentEstimator.eval(itb,ite,std::back_inserter(estimations));
}