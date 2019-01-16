#ifndef GEOC_ADAPTER_GLUEDCURVE_H
#define GEOC_ADAPTER_GLUEDCURVE_H

#include <geoc/adapter/base/GeneralAdapter.h>
#include "gcurve/GluedCurve.h"

namespace GEOC
{
    namespace Adapter
    {
        namespace GluedCurve
        {
            typedef GCurve::GluedCurve::MyGluedLinelsIterator IteratorType;

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

#endif //GEOC_ADAPTER_GLUEDCURVE_H
