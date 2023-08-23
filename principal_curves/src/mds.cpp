
#include "mds.hpp"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


using namespace std;

void CalculateMDS(double *SimilarityMatrix, int n, std::vector<std::pair<float, float> > *normalizedMDSpoints) {

    double* SquaredSimilarityMatrix = (double*)malloc(sizeof(double)*n*n);
    for(int i = 0; i < n; i++){
        for(int j = i; j < n; j++){

            SquaredSimilarityMatrix[i*n +j] = SimilarityMatrix[i*n +j]*SimilarityMatrix[i*n +j];
            SquaredSimilarityMatrix[j*n +i] = SimilarityMatrix[j*n +i]*SimilarityMatrix[j*n +i];
        }
    }

    CalculateMDS_Squared(SquaredSimilarityMatrix, n, normalizedMDSpoints);

    free(SquaredSimilarityMatrix);
}

void CalculateMDS_Squared(double* similarityMatrixSquared, int n, std::vector<std::pair<float,float>>* normalizedMDSpoints) {
    /*
    1. Set up the matrix of squared proximities P(2) = [p2].
    2. Apply the double centering:     ,
    where n is the number of objects.
    3. Extract the m largest positive eigenvalues λ1 . . . λm of B and the corresponding
    m eigenvectors e1 . . . em.
   */
    // squared proximities is in the weightedSimilarity measure
    gsl_matrix_view psquared = gsl_matrix_view_array (similarityMatrixSquared, n, n);

    gsl_matrix* J = gsl_matrix_alloc(n, n);

    for(int  i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            gsl_matrix_set(J, i, j, (-1.0/static_cast<float>(n))*1.0);
        }
        gsl_matrix_set(J, i, i, 1.0 - (1.0/static_cast<float>(n))*1.0);
    }


    gsl_matrix* t1 =  gsl_matrix_alloc(n,n);
    gsl_matrix* B = gsl_matrix_alloc(n,n);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, &psquared.matrix, J, 0.0, t1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, -0.5, J, t1, 0.0, B);

    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);

    gsl_eigen_symmv_workspace * w =  gsl_eigen_symmv_alloc (n);

    gsl_eigen_symmv (B, eval, evec, w);
    gsl_eigen_symmv_free (w);
    // positive eigenvalues
    gsl_eigen_symmv_sort (eval, evec,
            GSL_EIGEN_SORT_VAL_DESC);

    /*4. A m-dimensional spatial configuration of the n objects is derived from the
      coordinate matrix X = EmΛ1/2
      m , where Em is the matrix of m eigenvectors
      and Λm is the diagonal matrix of m eigenvalues of B, respectively.*/

    // eigen values has the diagonal matrix
    gsl_matrix* eigenvalues = gsl_matrix_alloc(2,2);
    gsl_matrix_set(eigenvalues, 0,1, 0);
    gsl_matrix_set(eigenvalues, 1,0, 0);
    for(int i =0; i < 2; i++)
    {
        // std::cout << "eigen val " << gsl_vector_get(eval,i) << std::endl;
        gsl_matrix_set(eigenvalues, i,i, sqrt(gsl_vector_get(eval,i)) );
    }
    for(int i = 0; i < n; i++) {
        printf("Eigenvalue %f\n", gsl_vector_get(eval, i));
    }
    //

    gsl_matrix* Em = gsl_matrix_alloc(n,2);
    for(int i = 0; i < 2; i++){
        // Get the first two eigen vectors
        gsl_vector_view evec_i = gsl_matrix_column (evec, i);

        // set each of them as a column
        for(int j = 0; j < n; j++){
            gsl_matrix_set(Em, j, i,  gsl_vector_get(&evec_i.vector,j) );
        }
    }

    gsl_matrix* finalPoints = gsl_matrix_alloc(n,2);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, Em, eigenvalues, 0.0, finalPoints);

    vector<pair<float,float>> MDSpoints;

    float MDSmin_dim1 = 1, MDSmax_dim1 = 1, MDSmin_dim2 = 1,MDSmax_dim2 = 1;

    for(int i = 0; i < n ; i++){

        float x = gsl_matrix_get(finalPoints, i, 0);
        float y = gsl_matrix_get(finalPoints, i, 1);

        if ( i == 0){
            MDSmin_dim1 = x;
            MDSmax_dim1 = x;
            MDSmin_dim2 = y;
            MDSmax_dim2 = y;
        }
        MDSmin_dim1 = min(MDSmin_dim1, x);
        MDSmax_dim1 = max(MDSmax_dim1, x);
        MDSmin_dim2 = min(MDSmin_dim2, y);
        MDSmax_dim2 = max(MDSmax_dim2, y);


        MDSpoints.push_back( make_pair(x,y));
    }
    normalizedMDSpoints->clear();

    for(int index = 0; index < n; index++){
        pair<float, float> nonNormalizedPoint = MDSpoints.at(index);
        // for range -1  to 1
        float x = ((nonNormalizedPoint.first - MDSmin_dim1 )/ ( MDSmax_dim1 - MDSmin_dim1));
        float y = ((nonNormalizedPoint.second - MDSmin_dim2 )/ ( MDSmax_dim2 - MDSmin_dim2));
        normalizedMDSpoints->push_back( make_pair(x,y));
    }

    gsl_matrix_free(Em);
    gsl_matrix_free(eigenvalues);
    gsl_matrix_free(t1);
    gsl_matrix_free(B);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
}

