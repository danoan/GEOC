#ifndef GEOC_API_GLUEDCURVE_TANGENT_H
#define GEOC_API_GLUEDCURVE_TANGENT_H

#include <geoc/estimator/adaptable/Tangent.h>
#include <geoc/adapter/GluedCurve.h>

namespace GEOC
{
    namespace API
    {
        namespace GluedCurve
        {
            class Tangent
            {
            public:

                template<typename T>
                using ALG_MDSS = GEOC::Estimator::Standard::MDSSTangent<T>;

                template<typename T>
                using ALD_LMDSS = GEOC::Estimator::Alternative::LMDSSTangent<T>;

            public:
                typedef DGtal::Z2i::KSpace KSpace;
                typedef GEOC::Adapter::GluedCurve::IteratorType CurveIterator;

                typedef DGtal::PointVector<2,double> TangentVector;
                typedef std::vector<TangentVector> EstimationsVector;


            private:
                Tangent();
                Tangent(const Tangent&)=delete;
                Tangent& operator=(const Tangent&)=delete;

            public:

                template< template<typename> class TAlgorithm>
                static void symmetricOpen(const KSpace& KImage,
                                          CurveIterator begin,
                                          CurveIterator end,
                                          EstimationsVector& ev,
                                          double h=1.0)
                {
                    GEOC::Adapter::GluedCurve::Symmetric< TAlgorithm,false>(begin,end,KImage,ev,h);
                }

                template< template<typename> class TAlgorithm>
                static void symmetricClosed(const KSpace& KImage,
                                            CurveIterator begin,
                                            CurveIterator end,
                                            EstimationsVector& ev,
                                            double h=1.0)
                {
                    GEOC::Adapter::GluedCurve::Symmetric< TAlgorithm,true>(begin,end,KImage,ev,h);
                }


            };
        }
    }
}
#endif //GEOC_API_GLUEDCURVE_TANGENT_H