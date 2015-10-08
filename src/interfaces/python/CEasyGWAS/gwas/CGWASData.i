%ignore GWASData;
%ignore GWASResults;
%ignore CGWASDataException;

%ignore CGWASDataHelper::encodeHomozygousData(GWASData*);
%ignore CGWASDataHelper::encodeHeterozygousData(GWASData*);
%ignore CGWASDataHelper::encodeHeterozygousData(GWASData*,uint const&);
%ignore CGWASDataHelper::filterSNPsByMAF;
%ignore CGWASDataHelper::filterUniqueSNPs;
%ignore CGWASDataHelper::createSNPHash;
%ignore CGWASDataHelper::removeSamples4MissingData;

%include "CEasyGWAS/gwas/CGWASData.h"
