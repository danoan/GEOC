#ifndef GEOC_API_GRIDCURVE_CURVATURE_H
#define GEOC_API_GRIDCURVE_CURVATURE_H

#include <DGtal/helpers/StdDefs.h>

#include <geoc/adapter/GridCurve.h>
#include <geoc/estimator/adaptable/Curvature.h>

namespace GEOC
{
    namespace API
    {
        namespace GridCurve
        {
            namespace Curvature
            {
                namespace EstimationAlgorithms
                {
                    template<typename T>
                    using ALG_MDCA = GEOC::Estimator::Standard::MDCACurvature<T>;

                    template<typename T>
                    using ALG_II = GEOC::Estimator::Standard::IICurvature<T>;

                    template<typename T>
                    using ALG_LMDCA = GEOC::Estimator::Alternative::LMDCACurvature<T>;

                    template<typename T>
                    using ALG_HMDCA = GEOC::Estimator::Alternative::HMDCACurvature<T>;
                }

                typedef DGtal::Z2i::Curve Curve;
                typedef Curve::ConstIterator CurveIterator;

                typedef std::vector<double> EstimationsVector;

                template< template<typename> class TAlgorithm>
                void symmetricOpen(const KSpace& KImage,
                                   CurveIterator begin,
                                   CurveIterator end,
                                   EstimationsVector& ev,
                                   double h,
                                   void* data)
                {
                    GEOC::Adapter::GridCurve::Symmetric< TAlgorithm,false>(begin,end,KImage,ev,h,data);
                }

                template< template<typename> class TAlgorithm>
                void symmetricClosed(const KSpace& KImage,
                                     CurveIterator begin,
                                     CurveIterator end,
                                     EstimationsVector& ev,
                                     double h,
                                     void* data)
                {
                    GEOC::Adapter::GridCurve::Symmetric< TAlgorithm,true>(begin,end,KImage,ev,h,data);
                }

                template< template<typename> class TAlgorithm>
                void identityOpen(const KSpace& KImage,
                                  CurveIterator begin,
                                  CurveIterator end,
                                  EstimationsVector& ev,
                                  double h,
                                  void* data)
                {
                    GEOC::Adapter::GridCurve::Identity< TAlgorithm,false>(begin,end,KImage,ev,h,data);
                }

                template< template<typename> class TAlgorithm>
                void identityClosed(const KSpace& KImage,
                                    CurveIterator begin,
                                    CurveIterator end,
                                    EstimationsVector& ev,
                                    double h,
                                    void* data)
                {
                    GEOC::Adapter::GridCurve::Identity< TAlgorithm,true>(begin,end,KImage,ev,h,data);
                }

            };
        }
    }
}

#endif //GEOC_API_GRIDCURVE_CURVATURE_H
