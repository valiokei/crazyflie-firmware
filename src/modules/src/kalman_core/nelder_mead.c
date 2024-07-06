#include "math.h"
#include "nelder_mead.h"
#include <stdlib.h>
#include <string.h>

/*
 * Implementation of Nelder-Mead method for direct uncontrained optimization
 *   mdl:  model to optimize
 *   opt:  optimization settings
 *   smpl: simplex used for optimization (updated on exit)
 *   out:  minimizer found by the optimization (output argument)
 */
void nelder_mead(const model *mdl, const optimset *opt, simplex *smpl, point *out)
{
    // initialize Nelder-Mead parameters
    const real alpha = 1.0f;
    const real gamma = opt->adaptive ? 1.0f + 2.0f / (float)smpl->n : 2.0f;
    const real rho = opt->adaptive ? 0.75f - 0.5f / (float)smpl->n : 0.5f;
    const real sigma = opt->adaptive ? 1.0f - 1.0f / (float)smpl->n : 0.5f;

    // simplex contains n + 1 vertices, where n is the dimensionality of each point
    for (int i = 0; i < smpl->n + 1; i++)
    {
        cost(mdl, smpl->vertices + i);
        smpl->num_eval++;
    }
    sort(smpl);

    // internal labels -
    point *best = smpl->vertices;
    point *worst = smpl->vertices + smpl->n;

    while (!terminated(smpl, opt))
    {
        smpl->num_iter++;
        int shrink = 0;

        update_centroid(smpl);
        update_simplex(alpha, &smpl->centroid, &smpl->centroid, worst,
                       mdl, smpl, &smpl->reflected);

        if (smpl->reflected.y < best->y)
        {
            update_simplex(gamma, &smpl->centroid, &smpl->reflected, &smpl->centroid,
                           mdl, smpl, &smpl->expanded);
            if (smpl->expanded.y < smpl->reflected.y)
            {
                copy_point(smpl->n, &smpl->expanded, worst);
            }
            else
            {
                copy_point(smpl->n, &smpl->reflected, worst);
            }
        }
        else if (smpl->reflected.y < (worst - 1)->y)
        {
            copy_point(smpl->n, &smpl->reflected, worst);
        }
        else if (smpl->reflected.y < worst->y)
        {
            update_simplex(rho, &smpl->centroid, &smpl->reflected, &smpl->centroid,
                           mdl, smpl, &smpl->contracted);
            shrink = smpl->contracted.y >= smpl->reflected.y;
            if (!shrink)
            {

                copy_point(smpl->n, &smpl->contracted, worst);
            }
        }
        else
        {
            update_simplex(rho, &smpl->centroid, worst, &smpl->centroid,
                           mdl, smpl, &smpl->contracted);
            shrink = smpl->contracted.y > worst->y;
            if (!shrink)
            {

                copy_point(smpl->n, &smpl->contracted, worst);
            }
        }

        if (shrink)
        {

            for (int i = 1; i < smpl->n + 1; i++)
            {
                update_simplex(sigma, best, smpl->vertices + i, best,
                               mdl, smpl, smpl->vertices + i);
            }
        }
        sort(smpl);
    }

    // save best point in output argument
    copy_point(smpl->n, best, out);
}

/*
 * Euclidean distance between two points
 */
real distance(int n, const point *pnt0, const point *pnt1)
{
    real sum = 0.0f;
    for (int i = 0; i < n; i++)
    {
        sum += SQR(pnt0->x[i] - pnt1->x[i]);
    }
    return sqrtf(sum);
}

/*
 * Simplex sorting
 */
int compare(const void *arg1, const void *arg2)
{
    const real f1 = ((const point *)arg1)->y;
    const real f2 = ((const point *)arg2)->y;
    return (f1 > f2) - (f1 < f2);
}

void sort(simplex *smpl)
{
    qsort((void *)(smpl->vertices), (int)smpl->n + 1, sizeof(point), compare);
}

/*
 * Initial point at centroid, all vertices equally spaced, trial points allocated sss
 */
void init_simplex(int n, real scale, const point *inp, simplex *smpl)
{
    // simplex *smpl = malloc(sizeof(simplex));
    // smpl->vertices = malloc((n + 1) * sizeof(point));
    smpl->n = n;
    for (int i = 0; i < n + 1; i++)
    { // simplex vertices
        // smpl->vertices[i].x = malloc(n * sizeof(real));
        // ASSERT(smpl->vertices[i].x);
        for (int j = 0; j < n; j++)
        { // coordinates
            smpl->vertices[i].x[j] = 0.0f;
        }
    }
    real b = 0.0f;
    for (int j = 0; j < n; j++)
    {
        real c = sqrtf(1.0f - b);
        smpl->vertices[j].x[j] = c;
        real r = -(1.0f / n + b) / c;
        for (int i = j + 1; i < n + 1; i++)
        {
            smpl->vertices[i].x[j] = r;
        }
        b += SQR(r);
    }
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < n; j++)
        {
            smpl->vertices[i].x[j] = scale * smpl->vertices[i].x[j] + inp->x[j];
        }
    }
    // smpl->reflected = init_point(n);
    // smpl->expanded = init_point(n);
    // smpl->contracted = init_point(n);
    // smpl->centroid = init_point(n);
    smpl->num_iter = smpl->num_eval = 0;
    // return smpl;
}

/*
 * Get centroid (average position) of simplex
 */
void update_centroid(simplex *smpl)
{
    for (int j = 0; j < smpl->n; j++)
    {
        smpl->centroid.x[j] = 0.0f;
        for (int i = 0; i < smpl->n; i++)
        {
            smpl->centroid.x[j] += smpl->vertices[i].x[j];
        }
        smpl->centroid.x[j] /= smpl->n;
    }
}

/*
 * Extend a point from inp0 in the scaled direction of (inp1 - inp2)
 */
void update_simplex(real scale, const point *inp0, const point *inp1, const point *inp2,
                    const model *mdl, simplex *smpl, point *out)
{

    for (int j = 0; j < smpl->n; j++)
    {
        out->x[j] = inp0->x[j] + scale * (inp1->x[j] - inp2->x[j]);
    }
    cost(mdl, out);
    smpl->num_eval++;
}

/*
 * Compute tolerance on x (Euclidean distance between coordinates of best and worst point)
 */
real tolerance_x(const simplex *smpl)
{
    return distance(smpl->n, smpl->vertices, smpl->vertices + smpl->n);
}

/*
 * Compute tolerance on y (difference between values of best and worst point)
 */
real tolerance_y(const simplex *smpl)
{
    return smpl->vertices[smpl->n].y - smpl->vertices[0].y;
}

/*
 * Terminate or continue?
 */
int terminated(const simplex *smpl, const optimset *opt)
{
    return (
               smpl->num_eval >= opt->max_eval ||
               smpl->num_iter >= opt->max_iter) ||
           (tolerance_x(smpl) <= opt->tol_x &&
            tolerance_y(smpl) <= opt->tol_y);
}

/*
 * Point utilities
 */
// point *init_point(int n)
// {
//     // point *pnt = malloc(sizeof(point));
//     // ASSERT(pnt);
//     // pnt->x = malloc(n * sizeof(real));
//     // ASSERT(pnt->x);
//     point pnt;
//     return *pnt;
// }

void copy_point(int n, const point *inp, point *out)
{
    // memcpy(out->x, inp->x, n * sizeof(real));
    for (int i = 0; i < n; i++)
    {
        out->x[i] = inp->x[i];
    }
    out->y = inp->y;
}

// /*
//  * Memory utilities
// //  */
// void free_simplex(simplex *smpl)
// {
//     for (int i = 0; i < smpl->n; i++)
//     {
//         free(smpl->vertices[i].x);
//     }
//     free(smpl->vertices);
//     free_point(&smpl->reflected);
//     free_point(&smpl->expanded);
//     free_point(&smpl->contracted);
//     free_point(&smpl->centroid);
//     free(smpl);
// }

// void free_point(point *pnt)
// {
//     free(pnt->x);
//     free(pnt);
// }
