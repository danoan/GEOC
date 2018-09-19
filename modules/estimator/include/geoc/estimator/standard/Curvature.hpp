#include "geoc/estimator/standard/Curvature.h"

using namespace GEOC::Estimator::Standard;

template<typename IteratorType>
MDCACurvature<IteratorType>::MDCACurvature(IteratorType itb,
                                           IteratorType ite,
                                           std::vector<double>& estimations,
                                           double h)
{
    SegmentComputer sc;
    SCEstimator f;

    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDCACurvatureEstimator(sc,f);

    MCMDCACurvatureEstimator.init(h,itb,ite);
    MCMDCACurvatureEstimator.eval(itb,ite,std::back_inserter(estimations));
}

template<typename IteratorType>
IICurvature<IteratorType>::IICurvature(IteratorType itb,
                                       IteratorType ite,
                                       std::vector<double>& estimations,
                                       double h,
                                       bool ccw)
{
    BoundingBox bb;
    DIPaCUS::Properties::CurveBoundingBox<IteratorType>(bb,itb,ite);

    Domain domain(bb.lb + DGtal::Z2i::Point(-2,-2),bb.ub+ DGtal::Z2i::Point(2,2));
    DigitalSet ds(domain);
    KSpace KImage;
    KImage.init(domain.lowerBound(),domain.upperBound(),true);



    DIPaCUS::Misc::CompactSetFromClosedCurve<IteratorType>(ds,itb,ite,ccw);

    double re_convolution_kernel = 3.0; // Euclidean radius of the convolution kernel. Set by user.


    MyIICurvatureFunctor curvatureFunctor; /// Functor used to convert volume -> curvature
    curvatureFunctor.init( h, re_convolution_kernel ); // Initialisation for a grid step and a given Euclidean radius of convolution kernel

    MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
    curvatureEstimator.attach( KImage, ds ); /// Setting a KSpace and a predicate on the object to evaluate
    curvatureEstimator.setParams( re_convolution_kernel/h ); /// Setting the digital radius of the convolution kernel
    curvatureEstimator.init( h, itb, ite); /// Initialisation for a given h

    std::back_insert_iterator< std::vector< double > > resultsIt( estimations ); /// output iterator for results of Integral Invariant curvature computation
    curvatureEstimator.eval( itb, ite, resultsIt ); /// Computation
}

