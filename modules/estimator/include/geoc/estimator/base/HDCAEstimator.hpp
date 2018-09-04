#include "HDCAEstimator.h"

using namespace GEOC::Estimator::Base;

// ----------------------- Standard services ------------------------------

// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator>
inline
HDCAEstimator<SegmentComputer, SCEstimator>
::HDCAEstimator() {}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator>
inline
HDCAEstimator<SegmentComputer, SCEstimator>
::HDCAEstimator(const SegmentComputer &aSegmentComputer,
                                   const SCEstimator &aSCEstimator)
        : myH(0), mySC(aSegmentComputer), mySCEstimator(aSCEstimator) {}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator>
inline
void
HDCAEstimator<SegmentComputer, SCEstimator>
::init(const double h, const ConstIterator &itb, const ConstIterator &ite) {

    myH = h;
    myBegin = itb;
    myEnd = ite;
    if (this->isValid())
        mySCEstimator.init(myH, myBegin, myEnd);
}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator>
inline
bool
HDCAEstimator<SegmentComputer, SCEstimator>::isValid() const {
    return ((myH > 0) && (DGtal::isNotEmpty(myBegin, myEnd)));
}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator>
template<typename OutputIterator>
inline
OutputIterator
HDCAEstimator<SegmentComputer, SCEstimator>
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
    std::vector<double> outValues(intervalLength);


    for(auto it=outValues.begin();it!=outValues.end();++it)
    {
        *it= 0;
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
                double v = mySCEstimator.eval(segMemberIt);
                outValues[baseIndex] = fabs(v)>fabs(outValues[baseIndex])?fabs(v):outValues[baseIndex];

                ++pos;
                ++segMemberIt;
            }while(segMemberIt!=segIt->end());

            ++segIt;
        } while (segIt != segItEnd);

    }

    for (int i=0;i<outValues.size();++i)
    {
        result = outValues[i];
    }

    return result;

}


// ------------------------------------------------------------------------
template<typename SegmentComputer, typename SCEstimator>
inline
typename HDCAEstimator<SegmentComputer, SCEstimator>::Quantity
PessimistMaximalSegmentEstimator<SegmentComputer, SCEstimator>
::eval(const ConstIterator &it) {

    throw DGtal::InputException("Not implemented");

    if (this->isValid()) {
        if (DGtal::isNotEmpty(it, myEnd)) {


            return 0;
        } else {
            std::cerr
                    << "[DGtal::PessimistMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& it)]"
                    << " ERROR. Iterator is invalid (==myEnd)." << std::endl;
            throw DGtal::InputException();
            return Quantity();
        }

    } else {
        std::cerr
                << "[DGtal::PessimistMaximalSegmentEstimator<SegmentComputer,SCEstimator>::eval(const ConstIterator& it)]"
                << " ERROR. Object is not initialized." << std::endl;
        throw DGtal::InputException();
        return Quantity();
    }
}
