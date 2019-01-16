#ifndef GEOC_API_GLUEDCURVE_LENGTH_H
#define GEOC_API_GLUEDCURVE_LENGTH_H

#include <geoc/estimator/non-adaptable/length/SinCos.hpp>
#include <geoc/estimator/non-adaptable/length/Projected.hpp>
#include <geoc/adapter/GluedCurve.h>

namespace GEOC
{
    namespace API
    {
        namespace GridCurve
        {
            namespace Length
            {
                namespace EstimationAlgorithms
                {
                    template<typename TTangentAdapter>
                    using ALG_PROJECTED = GEOC::Estimator::Standard::ProjectedLength<TTangentAdapter>;

                    template<typename TTangentAdapter>
                    using ALG_SINCOS = GEOC::Estimator::Alternative::SinCosLength<TTangentAdapter>;
                }

                typedef DGtal::Z2i::KSpace KSpace;
                typedef GEOC::Adapter::GluedCurve::IteratorType CurveIterator;
                typedef std::vector<double> EstimationsVector;


                template<typename TAlgorithm>
                void mdssOpen(const KSpace& KImage,
                              CurveIterator begin,
                              CurveIterator end,
                              EstimationsVector& ev,
                              double h=1.0)
                {
                    typedef GEOC::Adapter::GeneralAdapter::Symmetric<CurveIterator,
                            GEOC::Estimator::Standard::MDSSTangent,
                            false> TangentAdapter;
                    TAlgorithm<TangentAdapter>(begin,end,KImage,ev,h);
                }

                template<typename TAlgorithm>
                void mdssClosed(const KSpace& KImage,
                                CurveIterator begin,
                                CurveIterator end,
                                EstimationsVector& ev,
                                double h=1.0)
                {
                    typedef GEOC::Adapter::GeneralAdapter::Symmetric<CurveIterator,
                            GEOC::Estimator::Standard::MDSSTangent,
                            true> TangentAdapter;
                    TAlgorithm<TangentAdapter>(begin,end,KImage,ev,h);
                }

                template<typename TAlgorithm>
                void lmdssOpen(const KSpace& KImage,
                               CurveIterator begin,
                               CurveIterator end,
                               EstimationsVector& ev,
                               double h=1.0)
                {
                    typedef GEOC::Adapter::GeneralAdapter::Symmetric<CurveIterator,
                            GEOC::Estimator::Alternative::LMDSSTangent,
                            false> TangentAdapter;
                    TAlgorithm<TangentAdapter>(begin,end,KImage,ev,h);
                }

                template<typename TAlgorithm>
                void lmdssClosed(const KSpace& KImage,
                                 CurveIterator begin,
                                 CurveIterator end,
                                 EstimationsVector& ev,
                                 double h=1.0)
                {
                    typedef GEOC::Adapter::GeneralAdapter::Symmetric<CurveIterator,
                            GEOC::Estimator::Alternative::LMDSSTangent,
                            true> TangentAdapter;
                    TAlgorithm<TangentAdapter>(begin,end,KImage,ev,h);
                }

            }
        }
    }
}

#endif //GEOC_API_GLUEDCURVE_LENGTH_H