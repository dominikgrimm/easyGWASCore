%module CEasyGWAS

/*
*Include all headers
*/
%{
#define SWIG_FILE_WITH_INIT
#define SWIG
#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/types.h"
#include "CEasyGWAS/regression/CRegression.h"
#include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
#include "CEasyGWAS/gwas/CGWASData.h"
#include "CEasyGWAS/gwas/CMetaGWAS.h"
#include "CEasyGWAS/meta/CMetaAnalysis.h"
#include "CEasyGWAS/meta/CEffectSize.h"
#include "CEasyGWAS/kernel/CKernels.h"
#include "CEasyGWAS/stats/CChi2.h"
#include "CEasyGWAS/stats/CGaussian.h"
#include "CEasyGWAS/stats/CBeta.h"
#include "CEasyGWAS/stats/CFisherF.h"
#include "CEasyGWAS/stats/CStudentT.h"
#include "CEasyGWAS/stats/CGamma.h"
/*
*/
%}

%include "types.i"

/*
*Include std lib specific swig interfaces
*/
%include "exception.i"
%include "std_string.i"
%include "std_map.i"
%include "std_vector.i"

namespace std {
        %template(StringVector) vector<string>;
        %template(CharVector) vector<char>;
        %template(CharCharVector) vector< vector<char> >;
        %template(Uint64Vector) vector<uint64>;
        %template(UintVector) vector<uint>;
        %template(Float64Vector) vector<float64>;
}

/*
*Add external numpy interfaces to access python specific numpy types
*/
%include "External/numpy.i"


%init %{
        import_array();
%}

/*
*Include Eigen lib to convert between numpy and eigen types
*/
%include "eigen.i"


/*
*Include internal interfaces
*/
%include "regression/CRegression.i"
%include "gwas/CSingleTraitGWAS.i"
%include "meta/CMetaAnalysis.i"
%include "meta/CEffectSize.i"
%include "gwas/CGWASData.i"
%include "gwas/CMetaGWAS.i"
%include "kernel/CKernels.i"
%include "stats/CChi2.i"
%include "stats/CGaussian.i"
%include "stats/CBeta.i"
%include "stats/CFisherF.i"
%include "stats/CStudentT.i"
%include "stats/CGamma.i"
/*
*Ignore some global stuff
*/
