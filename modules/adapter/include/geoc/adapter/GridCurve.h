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

            template<template<typename> class TEstimator, bool closedCurve>
            using Symmetric = GeneralAdapter::Symmetric<IteratorType,TEstimator,closedCurve>;

            template<template<typename> class TEstimator, bool closedCurve>
            using Identity = GeneralAdapter::Identity<IteratorType,TEstimator,closedCurve>;
        }
    }
}

#endif //GEOC_ADAPTER_GRIDCURVE_H
