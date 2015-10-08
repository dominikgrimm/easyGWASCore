#ifndef CMISC
#define CMISC


#include "CEasyGWAS/globals.h"

/*
*Helper to perform an argsort operation on eigen Vector
*/

class ArgSort {
    private:
        VectorXd __vpointer;

    public:
        bool operator() (int64 const&, int64 const&);
        ArgSort(VectorXd const&);

        VectorXd getIndices();
};


#endif
