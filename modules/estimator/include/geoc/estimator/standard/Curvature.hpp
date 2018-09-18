#include "geoc/estimator/standard/Curvature.h"

using namespace GEOC::Estimator::Standard;

template<typename IteratorType>
MDCACurvature<IteratorType>::MDCACurvature(IteratorType itb,
                                           IteratorType ite,
                                           std::vector<double>& estimations)
{
    SegmentComputer sc;
    SCEstimator f;

    DGtal::MostCenteredMaximalSegmentEstimator<SegmentComputer,SCEstimator> MCMDCACurvatureEstimator(sc,f);

    MCMDCACurvatureEstimator.init(1.0,itb,ite);
    MCMDCACurvatureEstimator.eval(itb,ite,std::back_inserter(estimations));
}

template<typename IteratorType>
IICurvature<IteratorType>::IICurvature(IteratorType itb,
                                       IteratorType ite,
                                       std::vector<double>& estimations,
                                       bool ccw)
{
    typedef DGtal::Z2i::Domain Domain;

    typename DIPaCUS::Properties::CurveBoundingBox<IteratorType>::BoundingBox bb;
    DIPaCUS::Properties::CurveBoundingBox<IteratorType>(bb,itb,ite);

    Domain domain(bb.lb + DGtal::Z2i::Point(-2,-2),bb.ub+ DGtal::Z2i::Point(2,2));
    DigitalSet ds(domain);
    ds.clear();

    DIPaCUS::Misc::CompactSetFromClosedCurve<IteratorType>(ds,itb,ite,ccw);

    KSpace KImage;
    KImage.init(bb.lb,bb.ub,true);




    typedef DGtal::LightImplicitDigitalSurface< KSpace, DigitalSet > LightImplicitDigSurface;
    typedef DGtal::DigitalSurface< LightImplicitDigSurface > MyDigitalSurface;

    DGtal::SurfelAdjacency<KSpace::dimension> SAdj( true );
    KSpace::Surfel bel = DGtal::Surfaces<KSpace>::findABel( KImage, ds, 100000 );
    LightImplicitDigSurface LightImplDigSurf( KImage, ds, SAdj, bel );
    MyDigitalSurface digSurf( LightImplDigSurf );

    typedef DGtal::DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef DGtal::GraphVisitorRange< Visitor > VisitorRange;
    typedef VisitorRange::ConstIterator SurfelConstIterator;

    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();







    double h = 1.0;

    double re_convolution_kernel = 3.0; // Euclidean radius of the convolution kernel. Set by user.

    typedef DGtal::functors::IICurvatureFunctor<DGtal::Z2i::Space> MyIICurvatureFunctor;
    typedef DGtal::IntegralInvariantVolumeEstimator< DGtal::Z2i::KSpace, DGtal::Z2i::DigitalSet, MyIICurvatureFunctor > MyIICurvatureEstimator;
    typedef MyIICurvatureFunctor::Value Value;

    MyIICurvatureFunctor curvatureFunctor; /// Functor used to convert volume -> curvature
    curvatureFunctor.init( h, re_convolution_kernel ); // Initialisation for a grid step and a given Euclidean radius of convolution kernel

    MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
    curvatureEstimator.attach( KImage, ds ); /// Setting a KSpace and a predicate on the object to evaluate
    curvatureEstimator.setParams( re_convolution_kernel/h ); /// Setting the digital radius of the convolution kernel
    curvatureEstimator.init( h, abegin, aend); /// Initialisation for a given h

    std::vector< Value > results;
    std::back_insert_iterator< std::vector< Value > > resultsIt( results ); /// output iterator for results of Integral Invariant curvature computation
    curvatureEstimator.eval( abegin, aend, resultsIt ); /// Computation
}

