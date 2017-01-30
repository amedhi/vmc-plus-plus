/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-16 22:12:57
* Last Modified by:   amedhi
* Last Modified time: 2017-01-30 16:11:27
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <Eigen/Core>

#ifndef MATRIX_H
#define MATRIX_H

using Vector3i = Eigen::Vector3i;
using Vector3d = Eigen::Vector3d;
using Diagonal3d = Eigen::DiagonalMatrix<double,3>;
using RealMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;

#endif
