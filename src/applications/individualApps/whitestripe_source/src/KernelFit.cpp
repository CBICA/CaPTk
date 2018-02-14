// Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
// See LICENSE file (GPLv3)
// KernelFit/Source/KernelFit.cc

// This source file contains the definitions for the KernelFit objects.

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "KernelFit.h"

template<class T>
KernelFit1D<T>::KernelFit1D(const std::vector<T> &x, const std::vector<T> &y,
	const T &bandwidth){

	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), "
		"one or both input vectors are empty!");

	if ( x.size() != y.size() )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), input vectors "
		"must be equal in length!");

	if ( bandwidth <= 0.0 )
		throw KernelFitError("From KernelFit1D::KernelFit1D(), the bandwidth "
		"must be greater than zero!");

	_x = x;
	_y = y;
	_b = bandwidth * bandwidth; // squared ahead of time

}

template<class T>
std::vector<T> KernelFit1D<T>::Solve(const std::vector<T> &x){

	//
	// solve for the smooth profile through the data at all `x`
	//

	if ( x.empty() )
        throw KernelFitError("From KernelFit1D::Solve(), the input vector "
        "cannot be empty!");

	std::vector<T> f( x.size(), 0.0);

	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (int i = 0; i <  (int)x.size(); i++){

		T sum = 0.0;

		for (std::size_t j = 0; j < _x.size(); j++){

			T W   = Kernel(_x[j] - x[i]);
			f[i] += W * _y[j];
			sum  += W;
		}

		f[i] /= sum;
	}

	return f;
}

template<class T>
std::vector<T> KernelFit1D<T>::Solve(const std::vector<T> &x, T (*W)(T)){

    //
    // solve for the smooth profile through the data at all `x`
    // using an alternative kernel function `W`
    //

    if ( x.empty() )
        throw KernelFitError("From KernelFit1D::Solve(), the input vector "
        "cannot be empty!");

    std::vector<T> f( x.size(), 0.0);

    // omp_set_num_threads() should be called prior to here!
    #pragma omp parallel for shared(f)
    for (int i = 0; i <  (int)x.size(); i++){

        T sum = 0.0;

        for (std::size_t j = 0; j < _x.size(); j++){

            T WW  = W(_x[j] - x[i]);
            f[i] += WW * _y[j];
            sum  += WW;
        }

        f[i] /= sum;
    }

    return f;
}

template<class T>
std::vector<T> KernelFit1D<T>::StdDev(const std::vector<T> &x){

	//
    // Solve for the estimated standard deviation by evaluating
    // the profile *at* the raw data points.
    //

    if ( x.empty() )
        throw KernelFitError("From KernelFit1D::StdDiv(), the input vector "
        "cannot be empty!");

    // solve profile at data points
    std::vector<T> f = Solve( _x );

    // solve variance at data points
    std::vector<T> var(_x.size(), 0.0);
    for (std::size_t i = 0; i < _x.size(); i++)
      var[i] = (T)pow(_y[i] - f[i], 2.0);

    // solve for smooth curve through variance points
    KernelFit1D<T> profile(_x, var, _b);
    std::vector<T> stdev = profile.Solve(x);

    // take sqrt for standard deviation
    for (std::size_t i = 0; i < x.size(); i++)
        stdev[i] = sqrt(stdev[i]);

    return stdev;
}

template<class T>
KernelFit2D<T>::KernelFit2D(const std::vector<T> &x, const std::vector<T> &y,
	const std::vector<T> &z, const T &bandwidth){

	// Constructor for the KernelFit1D object. Save the `x`, `y` data
	// and set an initial value for the `bandwidth`.

	if ( x.empty() || y.empty() || z.empty() )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), one or more "
			"input vectors were empty!");

	if ( x.size() != y.size() || x.size() != z.size() )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), input vectors "
			"must be equal in length!");

	if ( bandwidth <= 0.0 )
		throw KernelFitError("From KernelFit2D::KernelFit2D(), the bandwidth "
			"must be greater than zero!");

	_x = x;
	_y = y;
	_z = z;
	_b = bandwidth * bandwidth; // square

}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Solve(const std::vector<T> &x,
	const std::vector<T> &y){

	//
	// solve for the smooth surface through the data at all (x, y)
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::Solve(), one or both of "
			"`x` and `y` were empty!");

	// initialize f[x][y] to zeros with proper dimensions
	std::vector< std::vector<T> > f(x.size(), std::vector<T>(y.size(), 0.0));

	// omp_set_num_threads() should be called prior to here!
	#pragma omp parallel for shared(f)
	for (int i = 0; i < (int)x.size(); i++)
	for (int j = 0; j < (int)y.size(); j++){

		T sum = 0.0;

		for (std::size_t k = 0; k < _x.size(); k++){

			T W      = Kernel(x[i] - _x[k], y[j] - _y[k]);
			f[i][j] += W * _z[k];
			sum     += W;
		}

		f[i][j] /= sum;
	}

	return f;
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::Solve(const std::vector<T> &x,
    const std::vector<T> &y, T (*W)(T, T)){

    //
    // solve for the smooth surface throught the xy data using alternative
    // kernel function `W`.
    //

    if ( x.empty() || y.empty() )
        throw KernelFitError("From KernelFit2D::Solve(), one or both of "
         "`x` and `y` were empty!");

    // initialize f[x][y] to zeros with proper dimensions
    std::vector< std::vector<T> > f( x.size(), std::vector<T>(y.size(), 0.0));

    // omp_set_num_threads() should be called prior to here!
    #pragma omp parallel for shared(f)
    for (int i = 0; i < (int)x.size(); i++)
    for (std::size_t j = 0; j < y.size(); j++){

        T sum = 0.0;

        for (std::size_t k = 0; k < _x.size(); k++){

            T WW     = W(x[i] - _x[k], y[j] - _y[k]);
            f[i][j] += WW * _z[k];
            sum     += WW;
        }

        f[i][j] /= sum;
    }

    return f;
}

template<class T>
std::vector< std::vector<T> > KernelFit2D<T>::StdDev(const std::vector<T> &x,
	const std::vector<T> &y){

	//
	// Solve for the estimated standard deviation by evaluating
	// the profile *at* the raw data points.
	//

	if ( x.empty() || y.empty() )
		throw KernelFitError("From KernelFit2D::StdDev(), one or both of the "
	    "input vectors were empty!");

	// initialize vector for profile at data points
	std::vector<T> f(_x.size(), 0.0);

	// solve profile at data points
	#pragma omp parallel for shared(f)
	for (int i = 0; i < (int)_x.size(); i++){

	    T sum = 0.0;

	    for (std::size_t j = 0; j < _x.size(); j++){

	        T W   = Kernel(_x[i] - _x[j], _y[i] - _y[j]);
	        f[i] += W * _z[j];
	        sum  += W;
	    }

	    f[i] /= sum;
	}

	// solve for variances at data points
	std::vector<T> var(_x.size(), 0.0);
	for (std::size_t i = 0; i < _x.size(); i++)
    var[i] = (T)pow(_z[i] - f[i], 2.0);

	// solve for smooth surface through variance points
	KernelFit2D<T> profile(_x, _y, var, _b);
	std::vector< std::vector<T> > stdev = profile.Solve(x, y);

	// take sqrt for standard deviation
	#pragma omp parallel for shared(stdev)
	for (int i = 0; i < (int)x.size(); i++)
	for (std::size_t j = 0; j < y.size(); j++)
		stdev[i][j] = sqrt(stdev[i][j]);

	return stdev;
}

// template classes
template class KernelFit1D<float>;
template class KernelFit2D<float>;
template class KernelFit1D<double>;
template class KernelFit2D<double>;
template class KernelFit1D<long double>;
template class KernelFit2D<long double>;
