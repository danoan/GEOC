#include "geoc/estimator/alternative/Curvature.h"

using namespace GEOC::Estimator::Alternative;

template<typename IteratorType>
DCALambdaCurvature<IteratorType>::DCALambdaCurvature(IteratorType itb,
                                                     IteratorType ite,
                                                     std::vector<double> &estimations)
{
    SegmentComputer sc;
    SCEstimator f;

    DCALambdaEstimator<SegmentComputer,SCEstimator, LambdaFunction> DCALambdaCurvatureEstimator(sc,f);

    DCALambdaCurvatureEstimator.init(1.0,itb,ite);
    DCALambdaCurvatureEstimator.eval(itb,ite,std::back_inserter(estimations));
}

template<typename IteratorType>
HDCACurvature<IteratorType>::HDCACurvature(IteratorType itb,
                                           IteratorType ite,
                                           std::vector<double> &estimations)
{
    SegmentComputer sc;
    SCEstimator f;

    HDCASegmentComputer,SCEstimator> PessimistMDCACurvature(sc,f);

    PessimistMDCACurvature.init(1.0,itb,ite);

    PessimistMDCACurvature.eval(itb,ite,std::back_inserter(estimations));
}