#ifndef GEOC_TESTCONVERGENCE_H
#define GEOC_TESTCONVERGENCE_H

#include <DGtal/helpers/StdDefs.h>
#include <DGtal/shapes/parametric/Ball2D.h>
#include <DGtal/shapes/GaussDigitizer.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>

#include <DIPaCUS/derivates/Misc.h>
#include <geoc/estimator/standard/Curvature.h>
#include <geoc/adapter/GridCurve.h>

namespace GEOC
{
    namespace Test
    {
        extern bool verbose;
        extern bool displayOutput;

        class TestConvergence
        {
        public:
            typedef DGtal::Z2i::DigitalSet DigitalSet;
            typedef DGtal::Z2i::Domain Domain;
            typedef DGtal::Z2i::KSpace KSpace;
            typedef DGtal::Z2i::Space Space;
            typedef DGtal::Z2i::Curve Curve;
            typedef Curve::SCell SCell;
            typedef KSpace::Point Point;

            typedef DGtal::Ball2D<Space> Ball2D;
            typedef DGtal::GaussDigitizer<Space, Ball2D > MyGaussDigitizer;

        private:
            struct CurveAndDS
            {
                CurveAndDS(Curve& c, DigitalSet& ds):curve(c),digitalSet(ds){}

                Curve curve;
                DigitalSet digitalSet;
            };

        public:
            //Compare values from MDCA and II curvature estimators.
            TestConvergence(double radius, double h);

        private:
            CurveAndDS prepareDigitalObjects(double radius, double h);

        };
    }
}

#endif //GEOC_TESTCONVERGENCE_H
