#ifndef GEOC_API_GRIDCURVE_TANGENT_H
#define GEOC_API_GRIDCURVE_TANGENT_H

#include <geoc/estimator/adaptable/Tangent.h>
#include <geoc/adapter/GridCurve.h>

namespace GEOC
{
    namespace API
    {
        namespace GridCurve
        {
            namespace Tangent
            {
                namespace EstimationAlgorithms
                {
                    template<typename T>
                    using ALG_MDSS = GEOC::Estimator::Standard::MDSSTangent<T>;

                    template<typename T>
                    using ALD_LMDSS = GEOC::Estimator::Alternative::LMDSSTangent<T>;
                }


                typedef DGtal::Z2i::KSpace KSpace;
                typedef DGtal::Z2i::Curve Curve;
                typedef Curve::ConstIterator CurveIterator;

                typedef DGtal::PointVector<2,double> TangentVector;
                typedef std::vector<TangentVector> EstimationsVector;

                template< template<typename> class TAlgorithm>
                void symmetricOpen(const KSpace& KImage,
                                   CurveIterator begin,
                                   CurveIterator end,
                                   EstimationsVector& ev,
                                   double h=1.0)
                {
                    GEOC::Adapter::GridCurve::Symmetric< TAlgorithm,false>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                void symmetricClosed(const KSpace& KImage,
                                     CurveIterator begin,
                                     CurveIterator end,
                                     EstimationsVector& ev,
                                     double h=1.0)
                {
                    GEOC::Adapter::GridCurve::Symmetric< TAlgorithm,true>(begin,end,KImage,ev,h);
                }


            };
        }
    }
}

#endif