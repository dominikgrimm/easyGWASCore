#include "CEasyGWAS/utils/CMisc.h"
#include <algorithm>
#include <vector>

ArgSort::ArgSort(VectorXd const& v) {
    __vpointer = v;
}

bool ArgSort::operator() (int64 const& x, int64 const& y) {
    return __vpointer(x) < __vpointer(y);
}

VectorXd ArgSort::getIndices() {
    VectorXd ids = VectorXd::LinSpaced(__vpointer.rows(),0,__vpointer.rows()-1);
    std::sort(ids.data(),ids.data()+ids.rows(), *this);
    return ids;
}
