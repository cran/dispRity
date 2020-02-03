/* 
 *  char.diff function
 * 
 *  Parts of this script are modifications of r-source/src/library/stats/src/distance.c
 *  from https://github.com/wch/r-source by Thomas Guillerme (guillert@tcd.ie)
 *
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2016   The R Core Team
 *  Copyright (C) 2002, 2004  The R Foundation
 *
 *  GNU General Public License available at https://www.R-project.org/Licenses/
 */


#include "dispRity.h"

// Calculating the Hamming character distance
static double R_Hamming(double *x, int nr, int nc, int i1, int i2, int translate)
{
    double diff, dist, vector1[nc], vector2[nc];
    int count, i, k;

    //Initialise the values
    count= 0;
    dist = 0;

    //Isolating the two comparable characters
    for(i = 0 ; i < nc ; i++) {
        if(both_non_NA(x[i1], x[i2])) {
            
            // Create the vectors
            vector1[count] = x[i1];
            vector2[count] = x[i2];

            //Increment the counter
            count++;
        }
        i1 += nr;
        i2 += nr;
    }

    // Normalising the characters
    if(translate) {
        Normalise_single_character(vector1, count);
        Normalise_single_character(vector2, count);
    }

    for(k = 0 ; k < count ; k++) { 
        diff = fabs(vector1[k] - vector2[k]);
        // Normalise the difference (Fitch like)
        if(diff > 1) {
            diff = 1;   
        }
        // Increment the differences
        if(!ISNAN(diff)) {
            dist += diff;
        }        
    }

    if(count == 0) {
        return NA_REAL;
    } else {
        if(translate) {
            dist = dist/(count - 1);
        } else {
            dist = dist/count;
        }
        return dist;
    }
}



// static double R_Hamming_bitwise(double *x, int nr, int nc, int i1, int i2, int translate)
// {
//     double diff, dist, character1[nc], character2[nc];
//     int count, i, k;

//     //Initialise the values
//     count= 0;
//     dist = 0;

//     //Isolating the two comparable characters
//     for(i = 0 ; i < nc ; i++) {
//         if(both_non_NA(x[i1], x[i2])) {
            
//             // Create the vectors
//             character1[count] = x[i1];
//             character2[count] = x[i2];

//             //Increment the counter
//             count++;
//         }
//         i1 += nr;
//         i2 += nr;
//     }

//     // // Normalising the characters
//     // if(translate) {
//     //     Normalise_single_character(vector1, count);
//     //     Normalise_single_character(vector2, count);
//     // }


//     // Using a pop count:
//     // https://en.wikipedia.org/wiki/Hamming_weight
//     int popcount64d(uint64_t diff)
//     {
//         int count;
//         for (count=0; x; count++)
//             x &= x - 1;
//         return count;
//     }












//     for(k = 0 ; k < count ; k++) { 
//         // Get the difference for the whole vector
//         diff = vector1[k] ^ vector2[k]
//         if(diff > 0) {
//             // Increment the difference
//             dist += 1 ;
//         }

//         diff = fabs(vector1[k] - vector2[k]);
//         // Normalise the difference (Fitch like)
//         if(diff > 1) {
//             diff = 1;   
//         }
//         // Increment the differences
//         if(!ISNAN(diff)) {
//             dist += diff;
//         }        
//     }

//     if(count == 0) {
//         return NA_REAL;
//     } else {
//         if(translate) {
//             dist = dist/(count - 1);
//         } else {
//             dist = dist/count;
//         }
//         return dist;
//     }
// }

// R_distance function (R::dist())
void R_distanceHamming(double *x, int *nr, int *nc, double *d, int *diag, int *method, int *translate)
{
    int dc, i, j;
    size_t  ij;  /* can exceed 2^31 - 1 */
    double (*distfun)(double*, int, int, int, int, int) = NULL;

    //Open MPI
#ifdef _OPENMP
    int nthreads;
#endif

    distfun = R_Hamming;

    dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */

#ifdef _OPENMP
    if (R_num_math_threads > 0)
    nthreads = R_num_math_threads;
    else
    nthreads = 1; /* for now */
    if (nthreads == 1) {
    /* do the nthreads == 1 case without any OMP overhead to see
       if it matters on some platforms */
    ij = 0;
    for(j = 0 ; j <= *nr ; j++)
        for(i = j+dc ; i < *nr ; i++)
        d[ij++] = distfun(x, *nr, *nc, i, j, *translate);
    }
    else
    /* This produces uneven thread workloads since the outer loop
       is over the subdiagonal portions of columns.  An
       alternative would be to use a loop on ij and to compute the
       i and j values from ij. */
#pragma omp parallel for num_threads(nthreads) default(none)    \
    private(i, j, ij)                       \
    firstprivate(nr, dc, d, method, distfun, nc, x)
    for(j = 0 ; j <= *nr ; j++) {
        ij = j * (*nr - dc) + j - ((1 + j) * j) / 2;
        for(i = j+dc ; i < *nr ; i++)
        d[ij++] = distfun(x, *nr, *nc, i, j, *translate);
    }
#else
    ij = 0;
    for(j = 0 ; j <= *nr ; j++)
    for(i = j+dc ; i < *nr ; i++)
        d[ij++] = distfun(x, *nr, *nc, i, j, *translate);
#endif
}


// R/C interface (former Diff)
SEXP C_diff_hamming(SEXP x, SEXP smethod, SEXP stranslate, SEXP attrs)
{
    // Define the variable
    SEXP result;
    int nr = nrows(x), nc = ncols(x), method = asInteger(smethod), translate = asInteger(stranslate);
    int diag = 0;
    R_xlen_t N;
    N = (R_xlen_t)nr * (nr-1)/2; /* avoid int overflow for N ~ 50,000 */
    PROTECT(result = allocVector(REALSXP, N));
    if(TYPEOF(x) != REALSXP) x = coerceVector(x, REALSXP);
    PROTECT(x);
    
    // Calculate the distance matrix
    R_distanceHamming(REAL(x), &nr, &nc, REAL(result), &diag, &method, &translate);
    
    // Wrap up the results
    SEXP names = PROTECT(getAttrib(attrs, R_NamesSymbol)); // Row/column names attributes
    for (int i = 0; i < LENGTH(attrs); i++) {
        setAttrib(result, install(translateChar(STRING_ELT(names, i))), VECTOR_ELT(attrs, i));
    }

    UNPROTECT(3);

    return result;
}
