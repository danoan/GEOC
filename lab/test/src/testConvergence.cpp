#include "TestConvergence.h"

using namespace GEOC::Test;

TestConvergence::CurveAndDS TestConvergence::prepareDigitalObjects(double radius, double h)
{
    Ball2D ball(0,0,radius);
    MyGaussDigitizer gd;
    gd.attach(ball);

    Domain unscaledDomain(ball.getLowerBound() + Point(-2,-2),
                          ball.getUpperBound() + Point(2,2));

    gd.init(unscaledDomain.lowerBound(),
            unscaledDomain.upperBound(),
            h);

    Domain domain = gd.getDomain();

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

    return CurveAndDS(c,ds);
}

TestConvergence::TestConvergence(double radius, double h)
{
    if(verbose) std::cout << "Test Convergence (radius=" << radius << ",h=" << h << ")" << std::endl;

    CurveAndDS cd = prepareDigitalObjects(radius,h);
    Curve& c = cd.curve;
    const Domain& domain = cd.digitalSet.domain();

    KSpace KImage;
    KImage.init(domain.lowerBound(),domain.upperBound(),true);

    typedef GEOC::Adapter::GridCurve::SymmetricCurvature< GEOC::Estimator::Standard::MDCACurvature,true > ClosedSymmetricCurvature;
    typedef GEOC::Adapter::GridCurve::IdentityRangeCurvature< GEOC::Estimator::Standard::IICurvature,false > OpenIICurvature;


    std::vector<double> estimationsMDCA;
    ClosedSymmetricCurvature CSC(c.begin(),c.end(),KImage,estimationsMDCA,h);

    std::vector<double> estimationsII;
    OpenIICurvature CIIC(c.begin(),c.end(),KImage,estimationsII,h);


    double mdca,ii;
    assert(estimationsII.size()==estimationsMDCA.size());

    double maxDiff=0;
    if(verbose)
    {
        for(int i=0;i<estimationsII.size();++i)
        {
            mdca = estimationsMDCA[i];
            ii = estimationsII[i];
            maxDiff = fabs(mdca-ii)>maxDiff?fabs(mdca-ii):maxDiff;
        }

        std::cout << "Max Diff: " << maxDiff << std::endl;
    }
}
