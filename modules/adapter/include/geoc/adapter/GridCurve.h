#ifndef GEOC_ADAPTER_GRIDCURVE_H
#define GEOC_ADAPTER_GRIDCURVE_H

#include <DGtal/helpers/StdDefs.h>
#include <geoc/adapter/base/GeneralAdapter.h>

namespace GEOC
{
    namespace Adapter
    {
        namespace GridCurve
        {
            typedef DGtal::Z2i::Curve Curve;
            typedef Curve::ConstIterator IteratorType;

            template<template<typename> class TEstimator, bool closedCurve >
            class IdentityRangeCurvature: public GeneralAdapter::IdentityRangeCurvature<IteratorType,TEstimator,closedCurve>
            {
            private:
               typedef GeneralAdapter::IdentityRangeCurvature<IteratorType,TEstimator,closedCurve> BaseClass;
            public:
                typedef typename BaseClass::EstimationValue EstimationValue;
            public:
                IdentityRangeCurvature(IteratorType begin,
                                       IteratorType end,
                                       std::vector<EstimationValue>& estimations,
                                       double h):BaseClass(begin,end,estimations,h)
                {}

            };

            template<template<typename> class TEstimator, bool closedCurve>
            class SymmetricCurvature:public GeneralAdapter::SymmetricCurvature<IteratorType,TEstimator,closedCurve>
            {
            private:
                typedef GeneralAdapter::SymmetricCurvature<IteratorType,TEstimator,closedCurve> BaseClass;
            public:
                typedef typename BaseClass::EstimationValue EstimationValue;
            public:
                SymmetricCurvature(IteratorType begin,
                                   IteratorType end,
                                   const KSpace& KImage,
                                   std::vector<EstimationValue>& estimations,
                                   double h):BaseClass(begin,end,KImage,estimations,h)
                {}
            };

            template<template<typename> class TEstimator, bool closedCurve>
            class SymmetricTangent:public GeneralAdapter::SymmetricTangent<IteratorType,TEstimator,closedCurve>
            {
            private:
                typedef GeneralAdapter::SymmetricTangent<IteratorType,TEstimator,closedCurve> BaseClass;
            public:
                typedef typename BaseClass::EstimationValue EstimationValue;
            public:
                SymmetricTangent(IteratorType begin,
                                 IteratorType end,
                                 const KSpace& KImage,
                                 std::vector<EstimationValue>& estimations,
                                 double h):BaseClass(begin,end,KImage,estimations,h)
                {}
            };
        }
    }
}

#endif //GEOC_ADAPTER_GRIDCURVE_H
