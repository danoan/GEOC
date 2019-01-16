#ifndef GEOC_ESTIMATOR_CURVATURE_II_H
#define GEOC_ESTIMATOR_CURVATURE_II_H

#include <DGtal/helpers/StdDefs.h>
#include <DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h>
#include <DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h>
#include <DIPaCUS/components/properties/boundingBox.h>
#include <DIPaCUS/components/Properties.h>
#include <DIPaCUS/derivates/Misc.h>

namespace GEOC
{
    namespace Estimator
    {
        namespace Standard
        {
            template<typename IteratorType>
            struct IICurvature
            {
                typedef DGtal::Z2i::KSpace KSpace;
                typedef DGtal::Z2i::DigitalSet DigitalSet;
                typedef DGtal::Z2i::Domain Domain;
                typedef DGtal::functors::IICurvatureFunctor<DGtal::Z2i::Space> MyIICurvatureFunctor;
                typedef DGtal::IntegralInvariantVolumeEstimator< KSpace, DigitalSet, MyIICurvatureFunctor > MyIICurvatureEstimator;

                typedef typename DIPaCUS::Properties::BoundingBox BoundingBox;

                //It is expected a counter-clockwise curve.
                IICurvature(IteratorType itb,
                            IteratorType ite,
                            std::vector<double>& estimations,
                            double h,
                            bool ccw=true)
                {
                    BoundingBox bb;
                    DIPaCUS::Properties::curveBoundingBox<IteratorType>(bb,itb,ite);

                    //The DGtal class that computes the II estimator needs the origin in the domain
                    if(bb.lb(0) >0 ) bb.lb(0) = 0;
                    if(bb.lb(1) >0 ) bb.lb(1) = 0;

                    if(bb.ub(0)< 0) bb.ub(0) = 0;
                    if(bb.ub(1)< 0) bb.ub(1) = 0;

                    Domain domain(bb.lb + DGtal::Z2i::Point(-2,-2),bb.ub+ DGtal::Z2i::Point(2,2));
                    DigitalSet ds(domain);

                    DIPaCUS::Misc::CompactSetFromClosedCurve<IteratorType>(ds,itb,ite,ccw);

                    KSpace KImage;
                    KImage.init(ds.domain().lowerBound(),ds.domain().upperBound(),true);

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
            };
        }
    }
}
#endif //GEOC_ESTIMATOR_CURVATURE_II_H
