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
            using Symmetric = GeneralAdapter::Symmetric<IteratorType,TEstimator,closedCurve>;

            template<template<typename> class TEstimator, bool closedCurve>
            using Identity = GeneralAdapter::Identity<IteratorType,TEstimator,closedCurve>;
            
        }
    }
}

#endif //GEOC_ADAPTER_GLUEDCURVE_H
