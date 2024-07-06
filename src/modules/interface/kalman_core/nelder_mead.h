
#include "stabilizer_types.h"

#define DIMENSIONS 3

// #include "model.h"

typedef float real;

/*
 * "Low-noise" squaring for arguments with no side-effects
 */
#define SQR(x) ((x) * (x))

/*
 * Function definition
 */
typedef struct Model model;

/*
 * Point definition
 *   x: n-dimensional array with point coordinates
 *   y: value of a function f applied to the coordinates x, y = f(x)
 */
typedef struct Point
{
    real x[DIMENSIONS];
    real y;
} point;

/*
 * Initialize function
 */
void init_model(float *tag_pos_predicted, float *tag_or_versor, voltMeasurement_t *voltAnchor, model *mdl);

/*
 * Return expected number of dimensions
 */
int dimensions(void);

/*
 * Cost function implementation
 *   model: model to optimize
 *   point: point where to evaluate the function
 */
void cost(const model *, point *);

/*
 * Optimizer settings
 */
typedef struct Optimset
{
    int precision;         // significant figures in floats/exponentials
    int format;            // fixed or exponential floating point format
    int verbose;           // toggle verbose output during minimization
    real tol_x;            // tolerance on the simplex solutions coordinates
    real tol_y;            // tolerance on the function value
    unsigned int max_iter; // maximum number of allowed iterations
    unsigned int max_eval; // maximum number of allowed function evaluations
    int adaptive;          // simplex updates reduced for dimension > 2
    real scale;            // scaling factor of initial simplex
} optimset;

/*
 * The "simplex" containing an array of n + 1 points each of dimension n
 */
typedef struct Simplex
{
    int n;
    unsigned int num_iter, num_eval;
    point vertices[DIMENSIONS + 1];
    point reflected;
    point expanded;
    point contracted;
    point centroid;
} simplex;

/*
 * "Simplex" or "Amoeba" optimizer
 */
void nelder_mead(const model *, const optimset *, simplex *, point *);

/*
 * Utility functions
 */
real distance(int, const point *, const point *);

int compare(const void *, const void *);

void sort(simplex *);

void init_simplex(int, real, const point *, simplex *smpl);

void update_centroid(simplex *);

void update_simplex(real, const point *, const point *, const point *,
                    const model *, simplex *, point *);

real tolerance_x(const simplex *);

real tolerance_y(const simplex *);

int terminated(const simplex *, const optimset *);

// point *init_point(int);

void copy_point(int, const point *, point *);

void free_simplex(simplex *);

void free_point(point *);