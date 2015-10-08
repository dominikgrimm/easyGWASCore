#ifndef CRegressionClass 
#define CRegressionClass

#include <string>
#include <Eigen/Dense>

#include "CEasyGWAS/globals.h"
#include "CEasyGWAS/optimizer/CRootFinding.h"

/*
*CRegression Exception Class
*/
class CRegressionException {
	private:
		std::string __error_msg;
	public:
		CRegressionException(std::string const& error_msg) : __error_msg(error_msg) {
			std::cout << RED << "CRegression Exception: " << error_msg << BLACK << "\n";
		}

		std::string what() {
			return __error_msg;
		}
};


/*
* CRegression Parent class
* Contains General Methods and Getter for different types of Regression Models
*/
class CRegression {
	
	protected:
		std::string _model_info;
		bool _intercept;
		MatrixXd _x;
		VectorXd _y;
		VectorXd _betas;
		VectorXd _yhat;
		VectorXd _residuals;
		uint64 _n_samples;
		uint64 _n_features;
		float64 _loglikelihood;
		uint _rank;

		void _checkDimensions(uint64*,uint64*) throw (CRegressionException);
	public:
		CRegression();
		CRegression(bool const&);
		~CRegression() {};
		
		virtual void fit() throw (CRegressionException) {};
		virtual void fit(bool const&) throw (CRegressionException) {};
		//virtual void fit(VectorXd const&, MatrixXd const&) throw (CRegressionException) {};
		//virtual void fit(VectorXd const&, MatrixXd const&, bool const&) throw (CRegressionException) {};

		//Several Getter Methods
		virtual MatrixXd getCovarianceBetas();
		virtual VectorXd getStdBetas();
		VectorXd getBetas() const;
		VectorXd getResiduals() const;
		virtual VectorXd getYHat();
		virtual float64 getRSquared() const;
		virtual float64 getAdjustedRSquared() const;
		virtual float64 getLogLikelihood() const;
		
		void setX(MatrixXd const&);
		void setY(VectorXd const&);
		
		float64 getAIC() const;
		float64 getAICc() const;
		float64 getBIC() const;

		uint getDF() const; //get degree of freedom
		
		void print();

};

/*
* CLinearRegression class: Inherits variables and methods from the general CRegression class
* This class fits a linear regression model assuming Gaussian noise using maximum likelihood estimates
*/
class CLinearRegression: public CRegression {
	
	private:
 		float64 __variance;

	protected:
		void _estimateBetas();
		void _estimateLogLikelihood();

	public:
		CLinearRegression();
		CLinearRegression(bool const&);
		~CLinearRegression() {};

		MatrixXd getCovarianceBetas();
		VectorXd getStdBetas();
		float64 getMSE();
		float64 getRMSE();
		
		void fit() throw (CRegressionException);
		void fit(bool const&) throw (CRegressionException);
		void fit(VectorXd const&, MatrixXd const&) throw (CRegressionException);
		void fit(VectorXd const&, MatrixXd const&, bool const&) throw (CRegressionException);
		void predict(VectorXd*,MatrixXd const&) throw (CRegressionException);
		
};

/*
* CLogisticRegression class: Inherits variables and methods from the general CRegression class
* This class fits a logistic regression model assuming Gaussian noise using maximum likelihood estimates
*/
class CLogisticRegression: public CRegression {
	
	private:
		float64 __epsilon;
		uint __iterations;
		MatrixXd __F;

	protected:
		void _estimateBetas();
		void _estimateLogLikelihood();

	public:
		CLogisticRegression();
		CLogisticRegression(bool const&);
		CLogisticRegression(bool const&,float64 const&);
		~CLogisticRegression() {};
		
		uint getIterations();
		MatrixXd getCovarianceBetas();
		VectorXd getStdBetas();
		VectorXd getYHat();

		void fit() throw (CRegressionException);
		void fit(bool const&) throw (CRegressionException);
		void fit(VectorXd const&, MatrixXd const&) throw (CRegressionException);
		void fit(VectorXd const&, MatrixXd const&, bool const&) throw (CRegressionException);
		void predict(VectorXd*,MatrixXd const&) throw (CRegressionException);
		

};

/*
* CLinearMixedRegression class: Inherits variables and methods from the general CRegression class
* This class fits a linear mixed regression model assuming Gaussian noise using maximum likelihood estimates
*/
class CLinearMixedRegression: public CRegression {

	private:
 		float64 __logDelta;
		bool __transform;
		bool __REML;
		bool __logDeltaSet;
		bool __use_brent;
		uint __nInterval;
		float64 __logDeltaMin;
		float64 __logDeltaMax;
		float64 __sigma;

		MatrixXd __U;
		MatrixXd __Ux;
		MatrixXd __Uy;
		MatrixXd __S;
		MatrixXd __K;
		MatrixXd __covarianceBetas;
		VectorXd __logDeltas;

		float64 __optimizeDelta();
		void __evaluateNLL(float64 const&,VectorXd*,float64*,float64*);
		

	public:
		CLinearMixedRegression();
		CLinearMixedRegression(bool const&);
		CLinearMixedRegression(MatrixXd const&, MatrixXd const&, MatrixXd const&);
		CLinearMixedRegression(MatrixXd const&, MatrixXd const&, MatrixXd const&, bool const&);
		~CLinearMixedRegression() {};

		VectorXd getYHat();
		MatrixXd getCovarianceBetas();
		VectorXd getStdBetas();
		float64 getMSE();
		float64 getRMSE();
		float64 getLogDelta();
		float64 getLogSigma();

		void fit() throw (CRegressionException);
		void fit(VectorXd const&, MatrixXd const&, MatrixXd const&) throw (CRegressionException);
		void predict(VectorXd*,MatrixXd const&, MatrixXd const&) throw (CRegressionException);
        VectorXd predict(MatrixXd const&,MatrixXd const&) throw (CRegressionException);
		
		void setK(MatrixXd const&);
		void setInterval(uint const&);
		void setLogDeltaMin(float64 const&);
		void setLogDeltaMax(float64 const&);
		void setLogDelta(float64 const&);
		void setREML(bool const&);
		void setBrent(bool const&);
};
		
//Function class template to evaluate brents method
class BrentFunction: public CBrentFunction {
	private:
		MatrixXd __Ux;
		MatrixXd __Uy;
		MatrixXd __S;
		bool __REML;
	public:	
		BrentFunction(MatrixXd const&, MatrixXd const&, MatrixXd const&, bool const&);
		float64 evaluate(float64 const&);
};

#endif //CRegressionClass
