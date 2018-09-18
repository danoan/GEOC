#include <DIPaCUS/derivates/Misc.h>
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
                                       std::vector<double>& estimations)
{
    typedef DGtal::SetOfSurfels<KSpace> DigitalSurfelContainer;
    typedef DGtal::DigitalSurface< DigitalSurfelContainer > MyDigitalSurface;


    auto it = itb;
    KSpace::Point minKP,maxKP,currP;
    minKP = it->preCell().coordinates;
    maxKP = minKP;
    do{
        currP = it->preCell().coordinates;

        minKP[0] = currP[0] < minKP[0]?currP[0]:minKP[0];
        minKP[1] = currP[1] < minKP[1]?currP[1]:minKP[1];

        maxKP[0] = currP[0] > minKP[0]?currP[0]:maxKP[0];
        maxKP[1] = currP[1] > minKP[1]?currP[1]:maxKP[1];

        ++it;
    }while(it!=ite);

    minKP/=2;
    maxKP/=2;

    KSpace KImage;
    DGtal::Z2i::Domain domain(minKP,maxKP);
    KImage.init(minKP,maxKP,true);

    DigitalSet ds(domain);
    DigitalSet boundary(domain);
    KSpace::SCell interiorPixel = KImage.sIndirectIncident(*itb, KImage.sOrthDir(*itb));

    it = itb;
    do{
        boundary.insert( KImage.sCoords( KImage.sIndirectIncident(*it,KImage.sOrthDir(*it)) ) );
        ++it;
    }while(it!=ite);

    DIPaCUS::Misc::FillInterior(ds,KImage.sCoords(interiorPixel),boundary);


    DigitalSurfelContainer dsc(KImage,DGtal::SurfelAdjacency<KSpace::dimension>( true ));
    KSpace::SCell surfelModel = KImage.sCell( KSpace::Point(1,1), true );
    for(auto it=ds.begin();it!=ds.end();++it)
    {
        dsc.surfelSet().insert( KImage.sCell(*it,surfelModel) );
    }

    MyDigitalSurface digSurf(dsc);

    /// Construction of the shape + digitalization
    double h = 1.0;




    typedef DGtal::DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef DGtal::GraphVisitorRange< Visitor > VisitorRange;
    typedef VisitorRange::ConstIterator SurfelConstIterator;

    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    /// Integral Invariant stuff
    //! [IntegralInvariantUsage]
    double re_convolution_kernel = 3.0; // Euclidean radius of the convolution kernel. Set by user.

    typedef DGtal::functors::IICurvatureFunctor<DGtal::Z2i::Space> MyIICurvatureFunctor;
    typedef DGtal::IntegralInvariantVolumeEstimator< DGtal::Z2i::KSpace, DGtal::Z2i::DigitalSet, MyIICurvatureFunctor > MyIICurvatureEstimator;
    typedef MyIICurvatureFunctor::Value Value;

    MyIICurvatureFunctor curvatureFunctor; /// Functor used to convert volume -> curvature
    curvatureFunctor.init( h, re_convolution_kernel ); // Initialisation for a grid step and a given Euclidean radius of convolution kernel

    MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
    curvatureEstimator.attach( KImage, ds ); /// Setting a KSpace and a predicate on the object to evaluate
    curvatureEstimator.setParams( re_convolution_kernel/h ); /// Setting the digital radius of the convolution kernel
    curvatureEstimator.init( h, abegin, aend ); /// Initialisation for a given h

    std::vector< Value > results;
    std::back_insert_iterator< std::vector< Value > > resultsIt( results ); /// output iterator for results of Integral Invariant curvature computation
    curvatureEstimator.eval( abegin, aend, resultsIt ); /// Computation
}

