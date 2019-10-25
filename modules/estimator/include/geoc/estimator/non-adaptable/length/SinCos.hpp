#ifndef GEOC_SINCOS_H
#define GEOC_SINCOS_H

#include <vector>
#include <DGtal/helpers/StdDefs.h>

namespace GEOC
{
    namespace Estimator
    {
        namespace Alternative
        {
            template<typename TTangentAdapter>
            class SinCosLength
            {
            public:
                typedef DGtal::Z2i::KSpace KSpace;

                typedef typename TTangentAdapter::MyIterator MyIterator;
                typedef typename TTangentAdapter::EstimationValue MyTangentValue;

                typedef double EstimationValue;

                typedef TTangentAdapter MyTangentAdapter;

            public:
                SinCosLength(MyIterator begin,
                             MyIterator end,
                             const KSpace &KImage,
                             std::vector<EstimationValue> &estimations,
                             double h,
                             void* data)
                {
                    std::vector< MyTangentValue > tangentEstimations;
                    MyTangentAdapter(begin,end,KImage,tangentEstimations,h,data);

                    auto it = begin;
                    int i = 0;
                    do
                    {
                        //1.0/(cos+sin) Length estimation
                        estimations.push_back( h*1.0/( fabs( tangentEstimations[i][0] ) + fabs( tangentEstimations[i][1] ) ) );
                        ++i;
                        ++it;
                    }while(it!=end);

                };
            };

        }
    }
}
#endif //GEOC_SINCOS_H
