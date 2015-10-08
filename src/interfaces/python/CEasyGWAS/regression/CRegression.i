
%ignore CLinearRegression::CLinearRegression(bool const&);
%ignore CLinearRegression::fit(bool const&);
%ignore CLinearRegression::fit(VectorXd const&, MatrixXd const&);
%ignore CLinearRegression::fit(VectorXd const&, MatrixXd const&, bool const&);

%ignore CLogisticRegression::CLogisticRegression(bool const&);
%ignore CLogisticRegression::fit(bool const&);
%ignore CLogisticRegression::fit(VectorXd const&, MatrixXd const&);
%ignore CLogisticRegression::fit(VectorXd const&, MatrixXd const&, bool const&);


#%ignore CLinearMixedRegression::CLinearMixedRegression(bool const&);
%ignore CLinearMixedRegression::CLinearMixedRegression(VectorXd const&, MatrixXd const&, MatrixXd const&);
%ignore CLinearMixedRegression::CLinearMixedRegression(VectorXd const&, MatrixXd const&, MatrixXd const&, bool const&);
%ignore CLinearMixedRegression::fit(VectorXd const&, MatrixXd const&, MatrixXd const&);
%ignore CLinearMixedRegression::predict(VectorXd*, MatrixXd const&, MatrixXd const&);


%include "CEasyGWAS/regression/CRegression.h"
