namespace CSingleTraitGWAS {

%ignore LinearRegression::LinearRegression(VectorXd const&,MatrixXd const&);
%ignore LinearRegression::LinearRegression(VectorXd const&,MatrixXd const&, MatrixXd const&);
%ignore LinearRegression::LRTEST;        
%ignore LinearRegression::FTEST; 
%ignore LinearRegression::setTestStatistic;

%ignore LogisticRegression::LogisticRegression(VectorXd const&,MatrixXd const&);
%ignore LogisticRegression::LogisticRegression(VectorXd const&,MatrixXd const&, MatrixXd const&);

%ignore EMMAX::EMMAX(VectorXd const&,MatrixXd const&, MatrixXd const&);
%ignore EMMAX::EMMAX(VectorXd const&,MatrixXd const&, MatrixXd const&,MatrixXd const&);
       
%ignore FaSTLMM::FaSTLMM(VectorXd const&,MatrixXd const&, MatrixXd const&);
%ignore FaSTLMM::FaSTLMM(VectorXd const&,MatrixXd const&, MatrixXd const&,MatrixXd const&);
}


%include "CEasyGWAS/gwas/CSingleTraitGWAS.h"
