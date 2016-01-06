#ifndef LOAD_MAT_HPP
#define LOAD_MAT_HPP

#include "inc/stdafx.hpp"

#include <matio.h>

#include "Eigen/Core"

vector<Eigen::VectorXf*> load_mnist_mat(string path, string type) {
    mat_t* matfp = Mat_Open(path.c_str(),MAT_ACC_RDONLY);
    if(matfp == nullptr) throw runtime_error("Error opening MAT file!");

    matvar_t* matvar = Mat_VarReadNextInfo(matfp);

    if(matvar == nullptr) throw runtime_error("No variable named found in `"+path+"` file");

    size_t n_rows = matvar->dims[0], n_cols = matvar->dims[1];

    double* data_ptr = new double[n_rows*n_cols];
    int start[2] = {0, 0};
    int stride[2] = {1, 1};
    int edge[2] = {(int)n_rows, (int)n_cols};
    Mat_VarReadData(matfp, matvar, data_ptr, start, stride, edge);
    vector<Eigen::VectorXf*> images;
    if(type == "images") {
        for(size_t j=0; j<n_cols; j++) {
            Eigen::VectorXf* v = new Eigen::VectorXf(n_rows);
            for(size_t i=0; i<n_rows; i++)
                (*v)[i] = data_ptr[j*n_rows+i];
            images.push_back(v);
        }
    } else if(type == "labels") {
        for(size_t i=0; i<n_rows; i++) {
            Eigen::VectorXf* v = new Eigen::VectorXf(10);
            v->setZero();
            (*v)[data_ptr[i]] = 1;
            images.push_back(v);
        }
    }
    else { throw runtime_error("Unknown type `"+type+"` specified!"); }
    delete[] data_ptr;
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return images;
}


#endif // LOAD_MAT_HPP

