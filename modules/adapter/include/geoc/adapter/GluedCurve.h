#ifndef GEOC_ADAPTER_GLUEDCURVE_H
#define GEOC_ADAPTER_GLUEDCURVE_H

#include "gcurve/GluedLinelsIterator.h"
#include "DIPaCUS/derivates/Misc.h"

namespace GEOC
{
    namespace Adapter
    {
        namespace GluedCurve
        {
            template<typename TEstimator>
            class SymmetricCurvature
            {
            public:
                typedef DGtal::Z2i::Curve Curve;

                typedef DGtal::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

                typedef DGtal::ConstRangeAdapter< GCurve::GluedLinelsIterator,
                        AdapterFunctor,
                        AdapterFunctor::Output > GluedRangeAdapter;


                typedef ConstRangeAdapter< Curve::ConstIterator,
                        AdapterFunctor,
                        AdapterFunctor::Output > CurveRangeAdapter;

            public:
                SymmetricCurvature(GCurve::GluedLinelsIterator begin,
                                   GCurve::GluedLinelsIterator end,
                                   const KSpace& KImage,
                                   std::vector<double>& estimations);
            };

            template<typename TEstimator>
            class SymmetricTangent
            {
            public:
                typedef DGtal::Z2i::Curve Curve;

                typedef DGtal::functors::SCellToIncidentPoints<KSpace> AdapterFunctor;

                typedef DGtal::ConstRangeAdapter< GCurve::GluedLinelsIterator,
                        AdapterFunctor,
                        AdapterFunctor::Output > GluedRangeAdapter;


                typedef ConstRangeAdapter< Curve::ConstIterator,
                        AdapterFunctor,
                        AdapterFunctor::Output > CurveRangeAdapter;

                typedef DGtal::PointVector<2,double> TangentVector;

                SymmetricTangent(GCurve::GluedLinelsIterator begin,
                                 GCurve::GluedLinelsIterator end,
                                 const KSpace& KImage,
                                 std::vector< TangentVector >& estimations);
            };

            template<typename TTangentAdapter>
            class ProjectedLength
            {
                typedef typename TTangentAdapter::TangentVector TangentVector;

                ProjectedLength(GCurve::GluedLinelsIterator begin,
                                GCurve::GluedLinelsIterator end,
                                const KSpace& KImage,
                                std::vector< double >& estimations);
            };


            template<typename TTangentAdapter>
            class SinCosLength
            {
                typedef typename TTangentAdapter::TangentVector TangentVector;

                SinCosLength(GCurve::GluedLinelsIterator begin,
                             GCurve::GluedLinelsIterator end,
                             const KSpace &KImage,
                             std::vector<double> &estimations);
            };

        }
    }
}

#endif //GEOC_ADAPTER_GLUEDCURVE_H
