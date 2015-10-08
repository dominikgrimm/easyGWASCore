%include "CEasyGWAS/types.h"

%{
        //#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
        #define SWIG_FILE_WITH_INIT
        #include "Eigen/Core"
%}

/*
*Convert Numpy Vector of size n x 1 into a eigen vector of size n x 1
*/
//INPUT
%typemap(in,fragment="NumPy_Fragments") VectorXd const & (VectorXd vec) {
        if(array_numdims($input)>2) {
                PyErr_SetString(PyExc_ValueError,"This is not a vector. Vector needs size n x 1 or n!");
        }
        int isNewObject = 0;
        PyArrayObject* input_array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &isNewObject);
        
        int rows = array_size(input_array,0);
        int cols = 0;
        if(array_numdims($input)<2) {
                cols = 1;
        } else {
                cols = array_size(input_array,1);
        }
        
        $1 = &vec;
        (*$1) = Eigen::Map<VectorXdNumpy>((float64*)array_data(input_array),rows,cols).cast<float64>();

        if(isNewObject) {
                Py_DECREF(input_array);
        }
}

/*
*Convert Eigen Vector of size n x 1 into a numpy vector of size n x 1
*/
//Output
%typemap(out,fragment="NumPy_Fragments") VectorXd {
        PyObject** out = &$result;
        VectorXd* in = &$1;
        npy_intp dims[2] = {in->rows(),in->cols()};
        *out = PyArray_SimpleNew(2,dims,NPY_DOUBLE);
        float64* data = (float64*)array_data(*out);
        //float64* data = (float64*)PyArray_DATA(*out);
        Eigen::Map<VectorXdNumpy>(data,dims[0],dims[1]) = (*in);
        //float64* data = (float64*)array_data(out);
        //Eigen::Map<VectorXdNumpy>(data,dims[0],dims[1]) = (*in);
        
        /*PyObject** out = &$result;
        VectorXd* in = &$1;
        npy_intp dims[2] = {in->rows(),in->cols()};
        *out = PyArray_SimpleNew(2,dims,PyArray_DOUBLE);
        //float64* data = (float64*)PyArray_DATA(*out);
        //Eigen::Map<VectorXdNumpy>(data,dims[0],dims[1]) = (*in);
        float64* data = (float64*)array_data(out);
        Eigen::Map<VectorXdNumpy>(data,dims[0],dims[1]) = (*in);
        */
}

/*
*Convert Numpy Matrix of size n x m into a eigen matrix of size n x m
*/
//INPUT
%typemap(in,fragment="NumPy_Fragments") MatrixXd const & (MatrixXd mat) {
        if(array_numdims($input)>2) {
                PyErr_SetString(PyExc_ValueError,"Only 1 or 2 dimensional arrays supported!");
        }
        int isNewObject = 0;
        PyArrayObject* input_array = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, &isNewObject);
        
        int rows = array_size(input_array,0);
        int cols = 0;
        if(array_numdims($input)<2) {
                cols = 1;
        } else {
                cols = array_size(input_array,1);
        }
        
        $1 = &mat;
        (*$1) = Eigen::Map<MatrixXdNumpy>((float64*)array_data(input_array),rows,cols).cast<float64>();

        if(isNewObject) {
                Py_DECREF(input_array);
        }
}

/*
*Convert Eigen Matrix of size n x m into a numpy matrix of size n x m
*/
//Output
%typemap(out,fragment="NumPy_Fragments") MatrixXd {
        PyObject** out = &$result;
        MatrixXd* in = &$1;
        npy_intp dims[2] = {in->rows(),in->cols()};
        *out = PyArray_SimpleNew(2,dims,NPY_DOUBLE);
        // *out = PyArray_SimpleNew(2,dims,PyArray_DOUBLE);
        //float64* data = (float64*)PyArray_DATA(*out);
        float64* data = (float64*)array_data(*out);
        Eigen::Map<MatrixXdNumpy>(data,dims[0],dims[1]) = (*in);
        //float64* data = (float64*)array_data(out);
        //Eigen::Map<VectorXdNumpy>(data,dims[0],dims[1]) = (*in);
        //float64* data = (float64*)PyArray_DATA(*out);
        //for(int i=0; i!= dims[0]; i++)
        //        for(int j=0;j!=dims[1];j++)
        //                data[i*dims[1]+j] = in->coeff(i,j);
}

/*
*Typechecks to support different notations and types
*/
/*
%typecheck(SWIG_TYPECHECK_FLOAT_ARRAY) MatrixXd, MatrixXd *, MatrixXd*, MatrixXd&, MatrixXd &,
                                       const MatrixXd, const MatrixXd&, const MatrixXd &,
                                       MatrixXd const, MatrixXd const&, MatrixXd const &,
                                       VectorXd, VectorXd *, VectorXd*, VectorXd&, VectorXd &,
                                       const VectorXd, const VectorXd&, const VectorXd &,
                                       VectorXd const, VectorXd const&, VectorXd const & {
        $1 = (array_type($input) == NPY_DOUBLE) || (array_type($input) == NPY_FLOAT);
}
*/
