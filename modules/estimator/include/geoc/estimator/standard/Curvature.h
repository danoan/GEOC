#ifndef GEOC_ESTIMATOR_STANDARD_CURVATURE_H
#define GEOC_ESTIMATOR_STANDARD_CURVATURE_H


#include <DIPaCUS/derivates/Misc.h>
#include <DIPaCUS/components/Properties.h>

#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/geometry/curves/StabbingCircleComputer.h"
#include <DGtal/geometry/curves/estimation/SegmentComputerEstimators.h>
#include "DGtal/geometry/curves/estimation/MostCenteredMaximalSegmentEstimator.h"

#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/ExplicitDigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"

namespace GEOC
{
    namespace Estimator
    {
        namespace Standard
        {
            template<typename IteratorType>
            struct MDCACurvature
            {
                typedef DGtal::StabbingCircleComputer<IteratorType> SegmentComputer;
                typedef DGtal::CurvatureFromDCAEstimator<SegmentComputer, false> SCEstimator;

                MDCACurvature(IteratorType itb,
                              IteratorType ite,
                              std::vector<double>& estimations,
                              double h);
            };

            template<typename IteratorType>
            struct IICurvature
            {
                typedef DGtal::Z2i::KSpace KSpace;
                typedef DGtal::Z2i::DigitalSet DigitalSet;
                typedef DGtal::Z2i::Domain Domain;
                typedef DGtal::functors::IICurvatureFunctor<DGtal::Z2i::Space> MyIICurvatureFunctor;
                typedef DGtal::IntegralInvariantVolumeEstimator< KSpace, DigitalSet, MyIICurvatureFunctor > MyIICurvatureEstimator;

                typedef typename DIPaCUS::Properties::CurveBoundingBox<IteratorType>::BoundingBox BoundingBox;

                //It is expected a counter-clockwise curve.
                IICurvature(IteratorType itb,
                            IteratorType ite,
                            std::vector<double>& estimations,
                            double h,
                            bool ccw=true);
            };
        }
    }
}

#include "Curvature.hpp"

#endif //GEOC_ESTIMATOR_STANDARD_CURVATURE_H
