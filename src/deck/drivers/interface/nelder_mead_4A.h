#ifndef __NELDER_MEAD_4A_H__
#define __NELDER_MEAD_4A_H__

// #include "MagneticDeck.h"
#include "FreeRTOS.h"
#include "task.h"
#include "MagneticDeck.h"

/*
   Multivariable optimization without derivatives/gradients.
   Implementation of Nelder-Mead simplex algorithm.
   http://www.scholarpedia.org/article/Nelder-Mead_algorithm

Do this:

#define NELDER_MEAD_IMPLEMENTATION

before you include this file in *one* C or C++ file to create the implementation.
You can disable the debug_log feature to remove some extra checks and includes.

#define NM_NO_DEBUG_LOG

This version is written by Justin Meiners (2021), derived from Matteo Maggioni's (2017) work.
The full MIT License is listed at the end of the file.
 */

// Types
//-----------------------------------------------------------------------------
// You can override these with #define if you really want to.
#define NM_RHO 1.0f

#define NM_CHI 2.0f

#define NM_GAMMA 0.5f

#define NM_SIGMA 0.5f
#define NM_REAL float
#define TEMP_POINT_COUNT_ 4 // default 4

#define DIMENSION 3

#define num_anchors_to_take_4A 4

typedef struct
{
    float Anchors_4A[num_anchors_to_take_4A][DIMENSION];
    float versore_orientamento_cf_4A[DIMENSION];
    float frequencies_4A[num_anchors_to_take_4A];
    float MeasuredVoltages_calibrated_4A[num_anchors_to_take_4A];
    float Gain_4A;
} myParams_t_4A;

// Cost function.

// Parameters:
// - number of variables
// - point (array of n values)
// - args is user  data
typedef NM_REAL (*nm_multivar_real_func_t_4A)(int, const NM_REAL *, void *);

// Parmeters:
// - tol_x: Terminate if any dimension of the simplex is smaller.
// - tol_fx: Terminate if we see improvement less than this amount.
// - max_iterations: Terminate if we exceed this number of iterations.
// - restarts: How many time to try improving after a termination.
// - debug_log: whether to show debug text.

typedef struct
{
    NM_REAL tol_x;
    NM_REAL tol_fx;
    int max_iterations;
    int restarts;
    int debug_log;
} nm_params_t_4A;

typedef struct
{
    int tol_satisfied;
    int iterations;
    NM_REAL min_fx;
} nm_result_t_4A;

// CONVENIENCE API
//-----------------------------------------------------------------------------

// These functions try to find the minimum using a variety of tricks and techniques on top of a simplex.
// They are supposed to be as much of a "black box" as possible.

// Parameters:
// - dimension: number of variables
// - initial: point to start search around (array of n values)
// - func: is a pointer to a cost function to optimize
// - args: optional user arguments passed to the function.
// - params: tolerance and iteration control.
// - out: the point which minimizes function (array of n values)
//
nm_result_t_4A nm_multivar_optimize_4A(
    int dimension,
    const NM_REAL *initial,
    const NM_REAL *initial_search_size,
    nm_multivar_real_func_t_4A func,
    void *args,
    const nm_params_t_4A *params,
    NM_REAL *out);

//-----------------------------------------------------------------------------
// DETAILED API
//-----------------------------------------------------------------------------

// If you want to write your own seach/restart strategy using simplex iteration
// this detailed API can help.

// An n-dimensional point x with it's function value fx.
typedef struct
{
    NM_REAL *x;
    NM_REAL fx;
} nm_simplex_pt_t_4A;

// An n-dimensional simplex with n+1 simplex points.
typedef struct
{
    int dimension;
    nm_simplex_pt_t_4A *p;
    NM_REAL *p_buffer;
    NM_REAL *temp_buffer;
} nm_simplex_t_4A;

void nm_simplex_init_4A(nm_simplex_t_4A *simplex, int dimension);
void nm_simplex_shutdown_4A(nm_simplex_t_4A *simplex);

// Guess the simplex size in each dimension using only the initial value.
void nm_guess_simplex_size_4A(int n, const NM_REAL *initial, NM_REAL *out_size);

// Place the simplex in a standard position around the initial point.
// This is roughly offseting the intiial point by each vector in the standard basis.
/*
   |\
   | \
   |  \
   0---
 */

void nm_simplex_position_around_4A(nm_simplex_t_4A *simplex, const NM_REAL *initial, const NM_REAL *size);

int nm_simplex_iterate_4A(
    nm_simplex_t_4A *simplex,
    nm_multivar_real_func_t_4A func,
    void *args,
    const nm_params_t_4A *params);

void nm_params_init_default_4A(nm_params_t_4A *params, int dimension);

/*
   Copyright 2017 Matteo Maggioni, 2022 Justin Meiners

   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#endif // __NELDER_MEAD_4A_H__
