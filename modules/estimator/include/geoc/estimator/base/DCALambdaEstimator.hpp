#include "geoc/estimator/base/DCALambdaEstimator.h"

using namespace GEOC::Estimator::Base;


// ----------------------- Standard services ------------------------------

// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction>
inline
DCALambdaEstimator<SegmentComputer, SCEstimator,LambdaFunction>
::DCALambdaEstimator() {}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction>
inline
DCALambdaEstimator<SegmentComputer, SCEstimator, LambdaFunction>
::DCALambdaEstimator(const SegmentComputer &aSegmentComputer,
                                const SCEstimator &aSCEstimator)
        : myH(0), mySC(aSegmentComputer), mySCEstimator(aSCEstimator) {}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction>
inline
void
DCALambdaEstimator<SegmentComputer, SCEstimator, LambdaFunction>
::init(const double h, const ConstIterator &itb, const ConstIterator &ite) {

    myH = h;
    myBegin = itb;
    myEnd = ite;
    if (this->isValid())
        mySCEstimator.init(myH, myBegin, myEnd);
}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction>
inline
bool
DCALambdaEstimator<SegmentComputer, SCEstimator, LambdaFunction>::isValid() const {
    return ((myH > 0) && (DGtal::isNotEmpty(myBegin, myEnd)));
}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction>
template<typename OutputIterator>
inline
OutputIterator
DCALambdaEstimator<SegmentComputer, SCEstimator, LambdaFunction>
::eval(const ConstIterator &itb, const ConstIterator &ite,
       OutputIterator result)
{

    Segmentation seg(myBegin, myEnd, mySC);
    seg.setSubRange(itb, ite);
    if ((myBegin != itb) || (myEnd != ite)) { //if subrange
        seg.setMode("MostCentered++");
    } else {//whole range
        seg.setMode("MostCentered");
    }

    int intervalLength = DGtal::rangeSize(itb,ite);
    std::vector<std::pair<double, double> > outValues(intervalLength);


    for(auto it=outValues.begin();it!=outValues.end();++it)
    {
        it->first = 0;
        it->second = 0;
    }

    if (this->isValid())
    {

        SegmentIterator segIt = seg.begin();
        SegmentIterator segItEnd = seg.end();


        do {
            mySCEstimator.attach(*segIt);

            int pos = 0;
            int segLength = DGtal::rangeSize(segIt->begin(), segIt->end());

            auto segMemberIt = segIt->begin();
            do
            {
                int baseIndex = DGtal::rangeSize(itb, segMemberIt)%intervalLength;
                double lbda = lambdaFunction((double) pos / segLength);
                outValues[baseIndex].first += lbda;
                outValues[baseIndex].second += lbda * mySCEstimator.eval(segMemberIt);

                ++pos;
                ++segMemberIt;
            }while(segMemberIt!=segIt->end());

            ++segIt;
        } while (segIt != segItEnd);

    }

    for (int i=0;i<outValues.size();++i)
    {
        if(outValues[i].first==0)
        {
            result = 0;
        }else
        {
            double x = outValues[i].second / outValues[i].first;
            result = x;
        }
    }

    return result;

}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator, typename LambdaFunction>
inline
typename DCALambdaEstimator<SegmentComputer, SCEstimator, LambdaFunction>::Quantity
DCALambdaEstimator<SegmentComputer, SCEstimator, LambdaFunction>
::eval(const ConstIterator &it) {

    throw DGtal::InputException("Not implemented");

    if (this->isValid()) {
        if (DGtal::isNotEmpty(it, myEnd)) {
            int pos;
            int segmentLength;
            lambdaFunction((double) pos / segmentLength);

            return 0;
        } else {
            std::cerr
                    << "[DGtal::LambdaMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& it)]"
                    << " ERROR. Iterator is invalid (==myEnd)." << std::endl;
            throw DGtal::InputException();
            return Quantity();
        }

    } else {
        std::cerr
                << "[DGtal::LambdaMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& it)]"
                << " ERROR. Object is not initialized." << std::endl;
        throw DGtal::InputException();
        return Quantity();
    }
}