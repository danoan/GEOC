#include <DGtal/shapes/parametric/Ball2D.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <geoc/adapter/GridCurve.h>
#include <geoc/estimator/standard/Curvature.h>

#include "DGtal/io/boards/Board2D.h"

template<class TPredicate>
struct WithDomain
{
    typedef TPredicate Predicate;

    WithDomain(TPredicate predicate,
               DGtal::Z2i::Domain domain):predicate(predicate),
                                          domain(domain){}

    TPredicate predicate;
    DGtal::Z2i::Domain domain;
};

WithDomain<DGtal::Z2i::DigitalSet> dd(double radius)
{
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Space Space;

    typedef DGtal::Ball2D<Space> Ball2D;
    typedef DGtal::Z2i::Curve Curve;
    typedef Curve::SCell SCell;
    typedef Curve::Point Point;

    Ball2D ball(0,0,radius);
    DGtal::GaussDigitizer<Space,Ball2D> gd;
    gd.attach(ball);

    DGtal::Z2i::Domain domain(ball.getLowerBound() + Point(-2,-2),
                              ball.getUpperBound() + Point(2,2));

    gd.init(domain.lowerBound(),
            domain.upperBound(),
            1);

    KSpace KImage;
    KImage.init(domain.lowerBound(),
                domain.upperBound(),
                true);

    std::vector<Point> vectorOfPoint;
    DGtal::SurfelAdjacency<KSpace::dimension> SAdj(true);
    SCell bel = DGtal::Surfaces<KSpace>::findABel(KImage,gd);
    DGtal::Surfaces<KSpace>::track2DBoundaryPoints(vectorOfPoint,KImage,SAdj,gd,bel);

    Curve c(KImage);
    c.initFromVector(vectorOfPoint);

    DGtal::Z2i::DigitalSet ds(domain);
    DIPaCUS::Misc::CompactSetFromClosedCurve<Curve::ConstIterator>(ds,c.begin(),c.end());

    WithDomain<DGtal::Z2i::DigitalSet> wd(ds,ds.domain());

    return wd;
}

template<class TShape>
struct GDFactory
{
    typedef DGtal::GaussDigitizer< DGtal::Z2i::Space, TShape > MyGaussDigitizer;



    WithDomain<MyGaussDigitizer> operator()(TShape shape, double h)
    {
        MyGaussDigitizer digShape;
        digShape.attach( shape );
        digShape.init( shape.getLowerBound() + DGtal::Z2i::Point::diagonal(-1), shape.getUpperBound() + DGtal::Z2i::Point::diagonal(1), h );

        WithDomain<MyGaussDigitizer> wd(digShape,digShape.getDomain());

        return wd;
    };
};


template<class TWithDomain>
int fromDGtal(TWithDomain wd)
{
    typedef typename TWithDomain::Predicate TPointPredicate;
    using namespace DGtal;

    TPointPredicate digShape = wd.predicate;
    Z2i::Domain domainShape = wd.domain;

    /// Construction of the shape + digitalization
    double h = 1.0;

    typedef LightImplicitDigitalSurface< Z2i::KSpace, TPointPredicate > LightImplicitDigSurface;
    typedef DigitalSurface< LightImplicitDigSurface > MyDigitalSurface;



    Z2i::KSpace KSpaceShape;
    bool space_ok = KSpaceShape.init( domainShape.lowerBound(), domainShape.upperBound(), true );
    if ( !space_ok )
    {
        trace.error() << "Error in the Khamisky space construction." << std::endl;
    }

    SurfelAdjacency<Z2i::KSpace::dimension> SAdj( true );
    Z2i::KSpace::Surfel bel = Surfaces<Z2i::KSpace>::findABel( KSpaceShape, digShape, 100000 );
    LightImplicitDigSurface LightImplDigSurf( KSpaceShape, digShape, SAdj, bel );
    MyDigitalSurface digSurf( LightImplDigSurf );

    typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef GraphVisitorRange< Visitor > VisitorRange;
    typedef typename VisitorRange::ConstIterator SurfelConstIterator;

    VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
    SurfelConstIterator abegin = range.begin();
    SurfelConstIterator aend = range.end();

    /// Integral Invariant stuff
    //! [IntegralInvariantUsage]
    double re_convolution_kernel = 3.0; // Euclidean radius of the convolution kernel. Set by user.

    typedef functors::IICurvatureFunctor<Z2i::Space> MyIICurvatureFunctor;
    typedef IntegralInvariantVolumeEstimator< Z2i::KSpace, TPointPredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
    typedef MyIICurvatureFunctor::Value Value;

    MyIICurvatureFunctor curvatureFunctor; /// Functor used to convert volume -> curvature
    curvatureFunctor.init( h, re_convolution_kernel ); // Initialisation for a grid step and a given Euclidean radius of convolution kernel

    MyIICurvatureEstimator curvatureEstimator( curvatureFunctor );
    curvatureEstimator.attach( KSpaceShape, digShape ); /// Setting a KSpace and a predicate on the object to evaluate
    curvatureEstimator.setParams( re_convolution_kernel/h ); /// Setting the digital radius of the convolution kernel
    curvatureEstimator.init( h, abegin, aend ); /// Initialisation for a given h

    std::vector< Value > results;
    std::back_insert_iterator< std::vector< Value > > resultsIt( results ); /// output iterator for results of Integral Invariant curvature computation
    curvatureEstimator.eval( abegin, aend, resultsIt ); /// Computation
    //! [IntegralInvariantUsage]

    /// Drawing results
    Value min = std::numeric_limits < Value >::max();
    Value max = std::numeric_limits < Value >::min();
    for ( unsigned int i = 0; i < results.size(); ++i )
    {
        if ( results[ i ] < min )
        {
            min = results[ i ];
        }
        else if ( results[ i ] > max )
        {
            max = results[ i ];
        }
    }
    Board2D board;
    VisitorRange range2( new Visitor( digSurf, *digSurf.begin() ) );
    abegin = range2.begin();

    typedef GradientColorMap< Value > Gradient;
    Gradient cmap_grad( min, max );
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );

    board << SetMode( (*abegin).className(), "Paving" );
    std::string specificStyle = (*abegin).className() + "/Paving";
    for ( unsigned int i = 0; i < results.size(); ++i )
    {
        Z2i::KSpace::SCell currentCell = KSpaceShape.sIndirectIncident( *abegin, *KSpaceShape.sOrthDirs( *abegin ) ); // We apply the color to the inner spel (more visible than surfel)
        board << CustomStyle( specificStyle, new CustomColors( Color::Black, cmap_grad( results[ i ] )))
              << currentCell;
        ++abegin;
    }
    board.saveSVG ( "example-integralinvariant2D.svg" );

}


void old();

int main()
{
    typedef DGtal::Ball2D<DGtal::Z2i::Space> Ball2D;

    Ball2D ball(0,0,10);
    double h = 1.0;

    GDFactory<Ball2D> GDF;
    fromDGtal( GDF(ball,h) );


    WithDomain<DGtal::Z2i::DigitalSet> ds = dd(10);
    fromDGtal(ds);

    old();


    return 0;
}

void old()
{
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Space Space;

    typedef DGtal::Ball2D<Space> Ball2D;
    typedef DGtal::Z2i::Curve Curve;
    typedef Curve::SCell SCell;
    typedef Curve::Point Point;

    Ball2D ball(0,0,10);
    DGtal::GaussDigitizer<Space,Ball2D> gd;
    gd.attach(ball);

    DGtal::Z2i::Domain domain(ball.getLowerBound() + Point(-2,-2),
                              ball.getUpperBound() + Point(2,2));

    gd.init(domain.lowerBound(),
            domain.upperBound(),
            1);

    KSpace KImage;
    KImage.init(domain.lowerBound(),
                domain.upperBound(),
                true);

    std::vector<Point> vectorOfPoint;
    DGtal::SurfelAdjacency<KSpace::dimension> SAdj(true);
    SCell bel = DGtal::Surfaces<KSpace>::findABel(KImage,gd);
    DGtal::Surfaces<KSpace>::track2DBoundaryPoints(vectorOfPoint,KImage,SAdj,gd,bel);

    Curve c(KImage);
    c.initFromVector(vectorOfPoint);


    typedef GEOC::Adapter::GridCurve::SymmetricCurvature< GEOC::Estimator::Standard::MDCACurvature,true > ClosedSymmetricCurvature;
    typedef GEOC::Adapter::GridCurve::IdentityRangeCurvature< GEOC::Estimator::Standard::IICurvature,false > OpenIICurvature;


    std::vector<double> estimationsMDCA;
    ClosedSymmetricCurvature CSC(c.begin(),c.end(),KImage,estimationsMDCA);

    std::vector<double> estimationsII;
    OpenIICurvature CIIC(c.begin(),c.end(),KImage,estimationsII);

    double mdca,ii;
    assert(estimationsII.size()==estimationsMDCA.size());
    for(int i=0;i<estimationsII.size();++i)
    {
        mdca = estimationsMDCA[i];
        ii = estimationsII[i];
        std::cout << "MDCA: " << mdca << "\tII: " << ii << "\tDIFF: " << fabs(mdca-ii) << std::endl;
    }

}