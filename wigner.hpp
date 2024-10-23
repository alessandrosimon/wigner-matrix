#pragma once

#include <iostream>
#include <complex>
#include <vector>
#include <memory>
#include <cassert>
#include <cstdio> 

#define LAPACK_COMPLEX_CUSTOM
#define lapack_complex_double std::complex<double>
#define lapack_complex_float std::complex<float>
#include <lapacke.h>

#include <Eigen/Dense>


using Eigen::MatrixXcd;
using Eigen::VectorXcd;


class WignerD
{
public:
  WignerD(int lmax) : _lmax{lmax}
  {
    // cache eigenvector for later
    for (int l=0; l<=lmax; l++)
      {
	this->jEigenvecs.push_back(this->getEigenvecs(l));
      }
  }

  VectorXcd eval_col(int l, int n, double phi, double theta, double chi)
  {
    // evaluate one column of the wigner matrix. n is the physical index
    // i.e. going from -l to +l

    assert(l <= this->_lmax && "l is larger than lmax");
    assert(std::abs(n) <= l && "|m| must not be larger than l");

    // compute exp factor for phi and chi
    VectorXcd factor = (VectorXcd::LinSpaced(2*l + 1, -l, l)* (-1.0)* phi*std::complex<double>(0., 1.0)).array().exp();
    factor *= std::exp(std::complex(0., 1.)* (-1.0)* chi*static_cast<double>(n));
 
    VectorXcd result(2*l + 1);
    const int limit = 2*l+1;
    double mu;
    
    for (int m=0; m<limit; m++)
      {
	std::complex<double> res = 0.0;
	for (int i=0; i<limit; i++)
	  {
	    mu = static_cast<double>(i-l);
	    res += std::exp(mu* theta*std::complex<double>(0., -1.))*(this->jEigenvecs[l](i,m))*std::conj(this->jEigenvecs[l](i,n+l));
	  }
	result(m) = res;
      }
    
    return conj(result.array()*factor.array());
  }


  template <int l, int idx>
  VectorXcd Delta(double phi, double theta, double chi)
  {
    std::cerr << "Need to provide specialization!" << std::endl;
    exit(-1);
    return VectorXcd(1);
  }

private:
  MatrixXcd getEigenvecs(int j)
  {
    VectorXcd J = getJ(j).conjugate();
    VectorXcd J0 = VectorXcd::Zero(2*j+1);
  
    MatrixXcd Jb(2*j+1, 2);
    Jb << J0 , J;

    const int n = 2*j+1; // Size of the matrix
    const int kd = 1; // Number of super/sub-diagonals
    const int ldab = 2*j+1;

    // Array to store eigenvalues
    double w[n];

    // Array to store eigenvectors (size = n * n)
    MatrixXcd z(n, n);

    // Info variable to store status
    int info;

    info = LAPACKE_zhbev(LAPACK_ROW_MAJOR, 'V', 'L', n, kd, Jb.data(), ldab, w, z.data(), n);

    if (info == 0){
#ifdef DEBUG
      std::cout << "Eigenvalues:" << std::endl;
      for (int i = 0; i < n; i++) {
	std::cout << w[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Eigenvectors:" << std::endl;
      std::cout << z << std::endl;
#endif
    } else {
      std::cerr << "Error occurred in zhbev. Info: " << info << std::endl;
      exit(-1);
    }
    return z;
  }
  
  VectorXcd getJ(int j)
  {
    const VectorXcd idx = VectorXcd::LinSpaced(2*j + 1, 1-j, j + 1);
    const auto x0 = idx.array() + j;
    const auto x1 = (j+1) - idx.array();
    VectorXcd result = -(x1*x0).sqrt()/std::complex(0., 2.0);
    result(2*j) = 0.0;
    return result;
  } 
  int _lmax;
  std::vector<MatrixXcd> jEigenvecs;
};

