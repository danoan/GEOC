#ifndef GEOC_ADAPTER_GRIDCURVE_H
#define GEOC_ADAPTER_GRIDCURVE_H

#include <DGtal/helpers/StdDefs.h>

#include "DIPaCUS/derivates/Misc.h"

namespace GEOC
{
    namespace Adapter
    {
        namespace GridCurve
        {
            template<typename TEstimator>
            class SymmetricCurvature
            {
            public:
                typedef DGtal::Z2i::Curve Curve;
                typedef DGtal::Z2i::KSpace KSpace;

                typedef DGtal::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

                typedef DGtal::ConstRangeAdapter< Curve::ConstIterator,
                        AdapterFunctor,
                        AdapterFunctor::Output > SCellToPointRangeAdapter;
            public:
                SymmetricCurvature(Curve::ConstIterator begin,
                                   Curve::ConstIterator end,
                                   const KSpace& KImage,
                                   std::vector<double>& estimations,
                                   bool closedCurve);

            };

            template<typename TEstimator>
            class SymmetricTangent
            {
            public:
                typedef DGtal::Z2i::Curve Curve;
                typedef DGtal::Z2i::KSpace KSpace;

                typedef DGtal::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

                typedef DGtal::ConstRangeAdapter< Curve::ConstIterator,
                        AdapterFunctor,
                        AdapterFunctor::Output > SCellToPointRangeAdapter;

                typedef DGtal::PointVector<2,double> TangentVector;

            public:
                SymmetricTangent(Curve::ConstIterator begin,
                                 Curve::ConstIterator end,
                                 const KSpace& KImage,
                                 std::vector< TangentVector >& estimationsTangent,
                                 bool closedCurve);
            };
        }
    }
}

#include "GridCurve.hpp"

#endif //GEOC_ADAPTER_GRIDCURVE_H
