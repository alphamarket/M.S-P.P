#ifndef NNCONF_HPP
#define NNCONF_HPP

#include "inc/stdafx.hpp"

#include <stdexcept>

#include "Eigen/Core"

typedef scalar input;
typedef scalar output;
typedef scalar weight;
typedef scalar (*activation_func_t)(scalar);

#define e_raise(s) throw std::runtime_error(s);

inline Eigen::VectorXf pow(const Eigen::VectorXf& v, scalar p) {
    Eigen::VectorXf o(v);
    for(int i = 0; i < v.rows(); i++)
        o(i) = pow(o(i), p);
    return o;
}

template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
inline void replace_nan(Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols>& m, scalar p) {
    for(size_t i = 0; i < m.rows(); i++)
        for(size_t j = 0; j < m.cols(); j++)
            if(isnan(m(i, j)))
                m(i, j) = p;
}

struct learning_configs {
    scalar
        _learning_rate,
        _momentum_rate,
        _regularization_rate;
    /**
     * @brief learning_rate construct a learning rate
     * @param r the rate for updating biases and layers' weights
     */
    learning_configs(scalar r) : learning_configs(r, r / 2, r / 4) { }
    /**
     * @brief learning_rate construct a learning rate
     * @param b the rate for updating biases' weights
     * @param l the rate for updating layers' weights
     */
    learning_configs(scalar lr, scalar mr, scalar rr) {
        this->_learning_rate = lr;
        this->_momentum_rate = mr;
        this->_regularization_rate = rr;
    }
    friend ostream& operator<<(ostream& os, const learning_configs& lc) {
        os
            <<"Learning rate: "<<lc._learning_rate<<endl
            <<"Momentum rate: "<<lc._momentum_rate<<endl
            <<"Regularization rate: "<<lc._regularization_rate<<endl;
        return os;
    }
};

#endif // NNCONF_HPP

