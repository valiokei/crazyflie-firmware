#include "nelder_mead_4A.h"
#include <stdlib.h>
#include <math.h>
#include "FreeRTOS.h"
#include "task.h"
#include "usec_time.h"

static nm_simplex_pt_t_4A my_p[DIMENSION + 1];
static NM_REAL my_p_buffer[(DIMENSION + 1) * DIMENSION];
static NM_REAL my_temp_buffer[TEMP_POINT_COUNT_ * DIMENSION];

void nm_params_init_default_4A(nm_params_t_4A *params, int dimension)
{
    params->tol_x = 0;
    params->tol_fx = 0.00001f;
    params->max_iterations = 1000 + 1000 * (int)sqrtf((float)dimension);

    params->restarts = 3;
    params->debug_log = 0;
}

void nm_simplex_init_4A(nm_simplex_t_4A *simplex, int dimension)
{
    int n = dimension;
    // simplex->p = malloc((n + 1) * sizeof(nm_simplex_pt_t_4A));
    simplex->p = my_p;
    simplex->dimension = n;

    // Simplex has n + 1 points where n is the dimension pf the data.
    // simplex->p_buffer = malloc((n + 1) * n * sizeof(NM_REAL));
    simplex->p_buffer = my_p_buffer;

    for (int i = 0; i < n + 1; i++)
    {
        simplex->p[i].x = simplex->p_buffer + i * n;
    }

    // simplex->temp_buffer = malloc(TEMP_POINT_COUNT_ * n * sizeof(NM_REAL));
    simplex->temp_buffer = my_temp_buffer;
}

void nm_simplex_shutdown_4A(nm_simplex_t_4A *simplex)
{
    // free(simplex->p_buffer);
    // free(simplex->p);
    // free(simplex->temp_buffer);
}

void nm_guess_simplex_size_4A(int n, const NM_REAL *initial, NM_REAL *out_size)
{
    for (int j = 0; j < n; ++j)
    {
        out_size[j] = (0.1f * fabsf(initial[j])) + 0.005f;
    }
}

void nm_simplex_position_around_4A(nm_simplex_t_4A *simplex, const NM_REAL *initial, const NM_REAL *size)
{
    int n = simplex->dimension;
    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            simplex->p[i].x[j] = initial[j] + ((i - 1 == j) ? size[j] : 0.0f);
        }
    }
}

// Simplex sorting
static int cmp_pt_value__4A(const void *arg1, const void *arg2)
{
    const nm_simplex_pt_t_4A *p1 = arg1;
    const nm_simplex_pt_t_4A *p2 = arg2;

    if (p1->fx < p2->fx)
    {
        return -1;
    }
    else if (p2->fx < p1->fx)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

static void simplex_sort_4A(nm_simplex_t_4A *simplex)
{
    qsort(simplex->p, simplex->dimension + 1, sizeof(nm_simplex_pt_t_4A), cmp_pt_value__4A);
}

// Get centroid (average position) of simplex
static void simplex_centroid_4A(const nm_simplex_t_4A *simplex, nm_simplex_pt_t_4A *centroid)
{
    int n = simplex->dimension;
    NM_REAL scale = 1.0f / (NM_REAL)n;
    // (last point is not included)
    for (int j = 0; j < n; j++)
    {
        centroid->x[j] = 0;
        for (int i = 0; i < n; i++)
        {
            centroid->x[j] += simplex->p[i].x[j];
        }
        centroid->x[j] *= scale;
    }
}

static void update_images_4A(nm_simplex_pt_t_4A *start, nm_simplex_pt_t_4A *end, int dimension, nm_multivar_real_func_t_4A func, void *args)
{
    while (start != end)
    {
        start->fx = func(dimension, start->x, args);
        ++start;
    }
}

static void simplex_shrink_4A(nm_simplex_t_4A *simplex)
{
    int n = simplex->dimension;
    for (int i = 1; i < n + 1; i++)
    {
        // for (int j = 0; j < n; j++)
        // {
        //     simplex->p[i].x[j] = simplex->p[0].x[j] + NM_SIGMA * (simplex->p[i].x[j] - simplex->p[0].x[j]);
        // }
        simplex->p[i].x[0] = simplex->p[0].x[0] + NM_SIGMA * (simplex->p[i].x[0] - simplex->p[0].x[0]);
        simplex->p[i].x[1] = simplex->p[0].x[1] + NM_SIGMA * (simplex->p[i].x[1] - simplex->p[0].x[1]);
        simplex->p[i].x[2] = simplex->p[0].x[2] + NM_SIGMA * (simplex->p[i].x[2] - simplex->p[0].x[2]);
    }
}

// Assess if simplex satisfies the minimization requirements

static int should_stop__4A(const nm_simplex_t_4A *simplex, int iterations, const nm_params_t_4A *params)
{

    if (iterations >= params->max_iterations && params->max_iterations > 0)
        return 1;

    int n = simplex->dimension;

    if (params->tol_fx > 0.0f)
    {
        NM_REAL best = fabsf(simplex->p[0].fx);
        NM_REAL worst = fabsf(simplex->p[n].fx);
        NM_REAL fractional_difference = 2.0f * (worst - best) / (worst + best);

        // printf("%f\n", fractional_difference);

        if (fractional_difference < params->tol_fx)
            return 1;
    }

    NM_REAL max_width = 0.0f;
    if (params->tol_x > 0.0f)
    {
        // for (int i = 1; i < n + 1; i++)
        // {
        //     for (int j = 0; j < n; j++)
        //     {
        //         const NM_REAL width = fabsf(simplex->p[0].x[j] - simplex->p[i].x[j]);
        //         if (width > max_width)
        //             max_width = width;
        //     }
        // }
        // return max_width < params->tol_x;
        for (int i = 1; i < 4; i++) // n + 1 dove n è 3
        {
            NM_REAL width;

            width = fabsf(simplex->p[0].x[0] - simplex->p[i].x[0]);
            if (width > max_width)
                max_width = width;

            width = fabsf(simplex->p[0].x[1] - simplex->p[i].x[1]);
            if (width > max_width)
                max_width = width;

            width = fabsf(simplex->p[0].x[2] - simplex->p[i].x[2]);
            if (width > max_width)
                max_width = width;
        }
        return max_width < params->tol_x;
    }
    else
    {
        return 0;
    }
}

// Update current point
static void update_point_4A(const nm_simplex_t_4A *simplex, const nm_simplex_pt_t_4A *centroid,
                            NM_REAL lambda, nm_simplex_pt_t_4A *point)
{
    int n = simplex->dimension;
    // for (int j = 0; j < n; j++)
    // {
    //     point->x[j] = (1.0f + lambda) * centroid->x[j] - lambda * simplex->p[n].x[j];
    // }
    point->x[0] = (1.0f + lambda) * centroid->x[0] - lambda * simplex->p[n].x[0];
    point->x[1] = (1.0f + lambda) * centroid->x[1] - lambda * simplex->p[n].x[1];
    point->x[2] = (1.0f + lambda) * centroid->x[2] - lambda * simplex->p[n].x[2];
}

// static void copy_point__4A(int n, const nm_simplex_pt_t_4A *src, nm_simplex_pt_t_4A *dst)
// {
//     for (int i = 0; i < n; ++i)
//     {
//         dst->x[i] = src->x[i];
//     }
//     dst->fx = src->fx;
// }

static void copy_point__4A(int n, const nm_simplex_pt_t_4A *src, nm_simplex_pt_t_4A *dst)
{
    dst->x[0] = src->x[0];
    dst->x[1] = src->x[1];
    dst->x[2] = src->x[2];
    dst->fx = src->fx;
}

static void swap_points__4A(nm_simplex_pt_t_4A *p1, nm_simplex_pt_t_4A *p2)
{
    nm_simplex_pt_t_4A temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

// #undef TEMP_POINT_COUNT

nm_result_t_4A nm_multivar_optimize_4A(
    int dimension,
    const NM_REAL *initial,
    const NM_REAL *initial_search_size,
    nm_multivar_real_func_t_4A func,
    void *args,
    const nm_params_t_4A *params,
    NM_REAL *out)
{

    nm_simplex_t_4A simplex;
    nm_simplex_init_4A(&simplex, dimension);

    nm_simplex_position_around_4A(&simplex, initial, initial_search_size);

    nm_result_t_4A r;
    r.iterations = nm_simplex_iterate_4A(&simplex, func, args, params);
    r.tol_satisfied = r.iterations < params->max_iterations;

    for (int restarts = 0; restarts <= params->restarts; ++restarts)
    {
        nm_simplex_position_around_4A(&simplex, simplex.p[0].x, initial_search_size);
        int iterations = nm_simplex_iterate_4A(&simplex, func, args, params);
        r.iterations += iterations;
        if (!r.tol_satisfied)
            r.tol_satisfied = r.iterations < params->max_iterations;
    }

    r.min_fx = simplex.p[0].fx;
    // save solution in output argument
    if (out)
    {
        // for (int i = 0; i < dimension; ++i)
        // {
        //     out[i] = simplex.p[0].x[i];
        // }
        out[0] = simplex.p[0].x[0];
        out[1] = simplex.p[0].x[1];
        out[2] = simplex.p[0].x[2];
    }
    nm_simplex_shutdown_4A(&simplex);
    return r;
}

int nm_simplex_iterate_4A(
    nm_simplex_t_4A *simplex,
    nm_multivar_real_func_t_4A func,
    void *args,
    const nm_params_t_4A *params)
{

    int n = simplex->dimension;
    // internal points
    nm_simplex_pt_t_4A point_r;
    nm_simplex_pt_t_4A point_e;
    nm_simplex_pt_t_4A point_c;
    nm_simplex_pt_t_4A centroid;

    point_r.x = simplex->temp_buffer + 0 * n;
    point_e.x = simplex->temp_buffer + 1 * n;
    point_c.x = simplex->temp_buffer + 2 * n;
    centroid.x = simplex->temp_buffer + 3 * n;

    // sort points in the simplex so that simplex.p[0] is the point having
    // minimum fx and simplex.p[n] is the one having the maximum fx
    update_images_4A(simplex->p, simplex->p + (n + 1), n, func, args);
    simplex_sort_4A(simplex);

    // compute the simplex centroid
    simplex_centroid_4A(simplex, &centroid);

    int iterations = 0;
    // uint32_t now_tick = xTaskGetTickCount();
    // uint64_t now = usecTimestamp();
    // continue minimization until stop conditions are met
    // TickType_t xTicksToWait = pdMS_TO_TICKS(9); // Converti 9 ms in tick

    // TickType_t xStartTick = xTaskGetTickCount(); // Ottieni il tick di partenza
    // vTaskDelay(M2T(1));
    while (!should_stop__4A(simplex, iterations, params))
    {

        if (iterations % 3 == 0)
        {
            vTaskDelay(M2T(1));
        }
        // uint64_t now = usecTimestamp();
        // if ((xTaskGetTickCount() - xStartTick) > xTicksToWait)
        // {
        //     // Cede il controllo ad altri task
        //     vTaskDelay(M2T(1));

        //     xStartTick = xTaskGetTickCount();
        // }

        int shrink = 0;

        // if (verbose)
        //     printf("Iteration %04d     ", iterations);

        update_point_4A(simplex, &centroid, NM_RHO, &point_r);
        point_r.fx = func(n, point_r.x, args);

        if (point_r.fx < simplex->p[0].fx)
        {
            update_point_4A(simplex, &centroid, NM_RHO * NM_CHI, &point_e);

            point_e.fx = func(n, point_e.x, args);
            if (point_e.fx < point_r.fx)
            {
                // expand
                // if (verbose)
                //     printf("expand          ");
                copy_point__4A(n, &point_e, simplex->p + n);
            }
            else
            {
                // reflect
                // if (verbose)
                //     printf("reflect         ");
                copy_point__4A(n, &point_r, simplex->p + n);
            }
        }
        else
        {
            if (point_r.fx < simplex->p[n - 1].fx)
            {
                // reflect
                // if (verbose)
                //     printf("reflect         ");
                copy_point__4A(n, &point_r, simplex->p + n);
            }
            else
            {
                if (point_r.fx < simplex->p[n].fx)
                {
                    update_point_4A(simplex, &centroid, NM_RHO * NM_GAMMA, &point_c);
                    point_c.fx = func(n, point_c.x, args);
                    if (point_c.fx <= point_r.fx)
                    {
                        // contract outside
                        // if (verbose)
                        //     printf("contract out    ");
                        copy_point__4A(n, &point_c, simplex->p + n);
                    }
                    else
                    {
                        // shrink
                        // if (verbose)
                        //     printf("shrink         ");
                        shrink = 1;
                    }
                }
                else
                {
                    update_point_4A(simplex, &centroid, -NM_GAMMA, &point_c);
                    point_c.fx = func(n, point_c.x, args);
                    if (point_c.fx <= simplex->p[n].fx)
                    {
                        // contract inside
                        // if (verbose)
                        //     printf("contract in     ");
                        copy_point__4A(n, &point_c, simplex->p + n);
                    }
                    else
                    {
                        // shrink
                        // if (verbose)
                        //     printf("shrink         ");
                        shrink = 1;
                    }
                }
            }
        }
        if (shrink)
        {
            simplex_shrink_4A(simplex);
            update_images_4A(simplex->p + 1, simplex->p + (n + 1), n, func, args);
            simplex_sort_4A(simplex);
        }
        else
        {
            for (int i = n - 1; i >= 0 && simplex->p[i + 1].fx < simplex->p[i].fx; i--)
            {
                swap_points__4A(simplex->p + (i + 1), simplex->p + i);
            }
        }
        simplex_centroid_4A(simplex, &centroid);
        ++iterations;

        // if (verbose)
        // {
        //     // print current minimum
        //     printf("[ ");
        //     for (int i = 0; i < n; i++)
        //     {
        //         printf("%.2f ", simplex->p[0].x[i]);
        //     }
        //     printf("]    %.2f \n", simplex->p[0].fx);
        // }

        // uint64_t end = usecTimestamp();
        // float diff = (float)(end - now) / 1000.0f;
        // float a = 0;
    }
    // uint32_t end_tick = xTaskGetTickCount();
    // uint64_t end = usecTimestamp();

    // float diff_in_S = (float)(end - now) / 1000000.0f;
    // float diff_in_ms = (float)(end - now) / 1000.0f;

    // float diff_in_S_tick = (float)(end_tick - now_tick) / configTICK_RATE_HZ;
    // float diff_in_ms_tick = (float)(end_tick - now_tick) * 1000.0f / configTICK_RATE_HZ;

    return iterations;
}
