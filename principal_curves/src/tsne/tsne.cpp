/*
 *
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY LAURENS VAN DER MAATEN ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */
//#include <stdafx.h>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime> // related to clock;
#include <iostream>
#include "vptree.h"
#include "sptree.h"
#include "tsne.hpp"
//#include "utils.h"

using namespace std;

// Perform t-SNE
// input, output: *X original data;
// input: N
// input, output: *Y output projection
// input: no_dims: the dimension of output;
// input: perplexity (from 5 to 50 is good)
// input theta 
// input rand_seed is the random seed if you want to skip random init rand_seed = -1; random, otherwise use it as a random seed
// input skip_random_init
// input max_iter
// input stop_lying_iter
// input mom_switch_iter
// when to call:
// where to call:
// Q: so that means if I want to apply this code, I have to modify them, because I already have distance matrix but don't have X, N, and D;
// Q: I also need to figure out how they initialize Y with no_dims // I got it, it doesn't require me anything to concern Y
// Q: my task is now to transform distance matrix to P, or row_P, column_P, val_P for sparse matrix; for the function ComputeGaussianPerplexity(X, N, D...) to Distance Matrix Compute
void TSNE::run(double* X, int N, int D, double* Y, int no_dims, double perplexity, double theta, int rand_seed,
        bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter) {

    // Set random seed
    if (skip_random_init != true) {
        if(rand_seed >= 0) {
            printf("Using random seed: %d\n", rand_seed);
            srand((unsigned int) rand_seed);
        } else {
            printf("Using current time as random seed...\n");
            srand(time(NULL));
        }
    }// Q: alles gut but what is the randomness used for? for initialize Y in lower dimensional space

    // Determine whether we are using an exact algorithm
    if(N - 1 < 3 * perplexity) { printf("Perplexity too large for the number of data points!\n"); exit(1); } // Q: what is N? and What is perplexity. why do you this condition
    printf("Using no_dims = %d, perplexity = %f, and theta = %f\n", no_dims, perplexity, theta);// Q: what is no_dims, and theta?
    bool exact = (theta == .0) ? true : false;

    // Set learning parameters
    float total_time = .0; // What is this?
    clock_t start, end; // Q: What are they?
    double momentum = .5, final_momentum = .8; // What are they? Why do you use them?
    double eta = 200.0; // What is this? Why do you use it?

    // Allocate some memory
    double* dY    = (double*) malloc(N * no_dims * sizeof(double)); // What is this? It's an array include N*no_dims of double, no_dims could be the number of lower dimensional space
    double* uY    = (double*) malloc(N * no_dims * sizeof(double)); // Then I guess they could be to store the original data points into lower dimension
    double* gains = (double*) malloc(N * no_dims * sizeof(double));
    if(dY == NULL || uY == NULL || gains == NULL) { printf("Memory allocation failed!\n"); exit(1); } // Q: just to check if allocating memory is successful!
    for(int i = 0; i < N * no_dims; i++)    uY[i] =  .0; // Q: initializing uY = 0.0;
    for(int i = 0; i < N * no_dims; i++) gains[i] = 1.0; // Q: intializing gains = 1.0;

    // Normalize input data (to prevent numerical problems)
    printf("Computing input similarities...\n");
    start = clock(); // Q: start the clock here to count time
    zeroMean(X, N, D); // Q: What is the function: to modify the original point to have zero mean;
    double max_X = .0;
    for(int i = 0; i < N * D; i++) {
        if(fabs(X[i]) > max_X) max_X = fabs(X[i]);
    } // it's like max_X = max|X(i)|; for all element in every coordinates;
    for(int i = 0; i < N * D; i++) X[i] /= max_X; // Q: Normalized every element in every coordinates into [-1.0, 1.0];

    // Compute input similarities for exact t-SNE
    double* P; unsigned int* row_P; unsigned int* col_P; double* val_P;
    if(exact) { // Q: if theta == 0.0;

        // Compute similarities
        printf("Exact?");
        P = (double*) malloc(N * N * sizeof(double)); // Q: it seems like P is a matrix having size N x N 
        if(P == NULL) { printf("Memory allocation failed!\n"); exit(1); } // Q: check if allocating memory successfully
        computeGaussianPerplexity(X, N, D, P, perplexity); // Q: what is this function? After computing this function I will have the output P;

        // Symmetrize input similarities
        printf("Symmetrizing...\n");
        int nN = 0;
        for(int n = 0; n < N; n++) {
            int mN = (n + 1) * N;
            for(int m = n + 1; m < N; m++) {
                P[nN + m] += P[mN + n];
                P[mN + n]  = P[nN + m];
                mN += N;
            }
            nN += N;
        }
        double sum_P = .0;
        for(int i = 0; i < N * N; i++) sum_P += P[i];
        for(int i = 0; i < N * N; i++) P[i] /= sum_P;
    }

    // Compute input similarities for approximate t-SNE
    else {

        // Compute asymmetric pairwise input similarities
        computeGaussianPerplexity(X, N, D, &row_P, &col_P, &val_P, perplexity, (int) (3 * perplexity));

        // Symmetrize input similarities
        symmetrizeMatrix(&row_P, &col_P, &val_P, N);
        double sum_P = .0;
        for(int i = 0; i < row_P[N]; i++) sum_P += val_P[i];
        for(int i = 0; i < row_P[N]; i++) val_P[i] /= sum_P;
    }
    end = clock(); // Q: End clock here

    // Lie about the P-values
    if(exact) { for(int i = 0; i < N * N; i++)        P[i] *= 12.0; }
    else {      for(int i = 0; i < row_P[N]; i++) val_P[i] *= 12.0; } // Why is 12.0 needed to multiply with P? and why 12.0? Is this a magic number?

    // Initialize solution (randomly)
    if (skip_random_init != true) {
        for(int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;
    } // Q: There we go the initialize solution

    // Perform main training loop
    if(exact) printf("Input similarities computed in %4.2f seconds!\nLearning embedding...\n", (float) (end - start) / CLOCKS_PER_SEC);
    else printf("Input similarities computed in %4.2f seconds (sparsity = %f)!\nLearning embedding...\n", (float) (end - start) / CLOCKS_PER_SEC, (double) row_P[N] / ((double) N * (double) N));
    start = clock();

    for(int iter = 0; iter < max_iter; iter++) {

        // Compute (approximate) gradient
        if(exact) computeExactGradient(P, Y, N, no_dims, dY); // it has Y and no_dims;
        else computeGradient(P, row_P, col_P, val_P, Y, N, no_dims, dY, theta);

        // Update gains
        for(int i = 0; i < N * no_dims; i++) gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
        for(int i = 0; i < N * no_dims; i++) if(gains[i] < .01) gains[i] = .01;

        // Perform gradient update (with momentum and gains)
        for(int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
        for(int i = 0; i < N * no_dims; i++)  Y[i] = Y[i] + uY[i];

        // Make solution zero-mean
        zeroMean(Y, N, no_dims);

        // Stop lying about the P-values after a while, and switch momentum
        if(iter == stop_lying_iter) {
            if(exact) { for(int i = 0; i < N * N; i++)        P[i] /= 12.0; }
            else      { for(int i = 0; i < row_P[N]; i++) val_P[i] /= 12.0; }
        }
        if(iter == mom_switch_iter) momentum = final_momentum;

        // Print out progress
        if (iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
            end = clock();
            double C = .0;
            if(exact) C = evaluateError(P, Y, N, no_dims);
            else      C = evaluateError(row_P, col_P, val_P, Y, N, no_dims, theta);  // doing approximate computation here!
            if(iter == 0)
                printf("Iteration %d: error is %f\n", iter + 1, C);
            else {
                total_time += (float) (end - start) / CLOCKS_PER_SEC;
                printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n", iter, C, (float) (end - start) / CLOCKS_PER_SEC);
            }
            start = clock();
        }
    }
    end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;

    // Clean up memory
    free(dY);
    free(uY);
    free(gains);
    if(exact) free(P);
    else {
        free(row_P); row_P = NULL;
        free(col_P); col_P = NULL;
        free(val_P); val_P = NULL;
    }
    printf("Fitting performed in %4.2f seconds.\n", total_time);
}

// Quynh's version of run with the input of square distance matrix DD
// this version only covers the exact case. The approximate case will be studied and updated later.
// Still have no guildline how to use this function? What are the inputs? What are the output?
// When to call the function???
void TSNE::Q_run(double* DD, int N, double* &Y, int no_dims, double perplexity, int rand_seed,
        bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter)
{
    // Set random seed
    if (skip_random_init != true) {
        if(rand_seed >= 0) {
            printf("Using random seed: %d\n", rand_seed);
            srand((unsigned int) rand_seed);
        } else {
            printf("Using current time as random seed...\n");
            srand(time(NULL));
        }
    }// Q: alles gut but what is the randomness used for? for initialize Y in lower dimensional space

    if(N - 1 < 3 * perplexity) { printf("Perplexity too large for the number of data points!\n"); exit(1); } // Q: what is N? and What is perplexity. why do you this condition

    // Set learning parameters
    float total_time = .0; // What is this?
    clock_t start, end; // Q: What are they?
    double momentum = .5, final_momentum = .8; // What are they? Why do you use them?
    double eta = 200.0; // What is this? Why do you use it?

    // Allocate some memory
    double* dY    = (double*) malloc(N * no_dims * sizeof(double)); // What is this? It's an array include N*no_dims of double, no_dims could be the number of lower dimensional space
    double* uY    = (double*) malloc(N * no_dims * sizeof(double)); // Then I guess they could be to store the original data points into lower dimension
    double* gains = (double*) malloc(N * no_dims * sizeof(double));
    if(dY == NULL || uY == NULL || gains == NULL) { printf("Memory allocation failed!\n"); exit(1); } // Q: just to check if allocating memory is successful!
    for(int i = 0; i < N * no_dims; i++)    uY[i] =  .0; // Q: initializing uY = 0.0;
    for(int i = 0; i < N * no_dims; i++) gains[i] = 1.0; // Q: intializing gains = 1.0;

    // Normalize input data (to prevent numerical problems)
    printf("Computing input similarities...\n");
    start = clock(); // Q: start the clock here to count time

    // Compute input similarities for exact t-SNE
    double* P;
    // Compute similarities
    P = (double*) malloc(N * N * sizeof(double)); // Q: it seems like P is a matrix having size N x N 
    if(P == NULL) { printf("Memory allocation failed!\n"); exit(1); } // Q: check if allocating memory successfully
    computeGaussianPerplexity(DD, N,  P, perplexity); // Q: it's my version of computing gaussian perplexity

    // Symmetrize input similarities
    printf("Symmetrizing...\n");
    int nN = 0;
    for(int n = 0; n < N; n++) {
        int mN = (n + 1) * N;
        for(int m = n + 1; m < N; m++) {
            P[nN + m] += P[mN + n];
            P[mN + n]  = P[nN + m];
            mN += N;
        }
        nN += N;
    }
    double sum_P = .0;
    for(int i = 0; i < N * N; i++) sum_P += P[i];
    for(int i = 0; i < N * N; i++) P[i] /= sum_P;

    end = clock(); // Q: End clock here

    // Lie about the P-values
    for(int i = 0; i < N * N; i++) P[i] *= 12.0; 

    Y = (double*) malloc(N * no_dims * sizeof(double)); // Q's code
    if(Y == NULL) { printf("Memory allocation failed! \n"); exit(1); }; // Q's code


    // Initialize solution (randomly)
    if (skip_random_init != true) {
        for(int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;
    } // Q: There we go the initialize solution which I should choose false, because I want it to be random;

    // Perform main training loop
    printf("Input similarities computed in %4.2f seconds!\nLearning embedding...\n", (float) (end - start) / CLOCKS_PER_SEC);

    start = clock();
    for(int iter = 0; iter < max_iter; iter++) {

        // Compute (approximate) gradient
        computeExactGradient(P, Y, N, no_dims, dY); // it has Y and no_dims

        // Update gains
        for(int i = 0; i < N * no_dims; i++) gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
        for(int i = 0; i < N * no_dims; i++) if(gains[i] < .01) gains[i] = .01;

        // Perform gradient update (with momentum and gains)
        for(int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
        for(int i = 0; i < N * no_dims; i++)  Y[i] = Y[i] + uY[i];

        // Make solution zero-mean
        zeroMean(Y, N, no_dims);

        // Stop lying about the P-values after a while, and switch momentum
        if(iter == stop_lying_iter) {
            for(int i = 0; i < N * N; i++) P[i] /= 12.0; 
        }
        if(iter == mom_switch_iter) momentum = final_momentum;

        // Print out progress
        if (iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
            end = clock();
            double C = .0;
            C = evaluateError(P, Y, N, no_dims);
            if(iter == 0)
                printf("Iteration %d: error is %f\n", iter + 1, C);
            else {
                total_time += (float) (end - start) / CLOCKS_PER_SEC;
                printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n", iter, C, (float) (end - start) / CLOCKS_PER_SEC);
            }
            start = clock();
        }
    }
    end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;

    // Clean up memory
    free(dY);
    free(uY);
    free(gains);
    free(P);
    printf("Fitting performed in %4.2f seconds.\n", total_time);
}

// Quynh's function which transfer the distance matrix to square distance matrix 
// input: Distance matrix in form of constant matrix
// output: double* DD, and int &N which N is the number of element which is also the size of DD; and DD is the output square distance matrix
void TSNE::Q_transfer_DistMat_To_SquareDistMat_Pointer(const vector<vector<float>> &DistanceMatrix, double* DD, int &N)
{
    N = (int) DistanceMatrix.size(); // get the size of the distance matrix and output to N;
    vector<vector<float>> _distanceMat_tSNE;
    for(int i = 0; i < N; i++)
    {
        vector<float> _row_I;
        for(int j = 0; j < N; j ++)
        {
            _row_I.push_back(DistanceMatrix[i][j]);
        }
        _distanceMat_tSNE.push_back(_row_I);
    }// it's not necessary if I use Matrixd in class TSNE

    // normalized _distanceMat_tSNE first
    /*double _maxDis = -DBL_MAX;
      double _minDis =  DBL_MAX;
      for(int i = 0; i < N; i++)
      {
      for(int j = 0; j < N; j ++)
      {
      if(i != j)
      {
      if(_maxDis < _distanceMat_tSNE[i][j]) _maxDis = _distanceMat_tSNE[i][j];
      if(_minDis > _distanceMat_tSNE[i][j]) _minDis = _distanceMat_tSNE[i][j];
      }
      }
      }
    // modified the distance matrix
    if(_maxDis != _minDis)
    {
    double Diff = 1.0/(_maxDis - _minDis);
    for(int i = 0; i <N; i ++)
    {
    for(int j = 0; j < N; j ++)
    {
    if(i != j)
    {
    _distanceMat_tSNE[i][j] = 1.0 - (_distanceMat_tSNE[i][j] - _minDis)*Diff;
    }
    }
    }
    }// Done modified the distance matrix;*/

    //Now assign element to DD
    for(int i = 0; i < N; i ++) // i is the row index
    {
        for(int j = 0; j < N; j ++) // j is the column index
        {
            DD[i*N + j] = (_distanceMat_tSNE[i][j])*(_distanceMat_tSNE[i][j]);// the result stay the same if one transpose i to j because Distnacematrix is a symetric one.
        }
    }


}

// Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
void TSNE::computeGradient(double* P, unsigned int* inp_row_P, unsigned int* inp_col_P, double* inp_val_P, double* Y, int N, int D, double* dC, double theta)
{

    // Construct space-partitioning tree on current map
    SPTree* tree = new SPTree(D, Y, N);

    // Compute all terms required for t-SNE gradient
    double sum_Q = .0;
    double* pos_f = (double*) calloc(N * D, sizeof(double));
    double* neg_f = (double*) calloc(N * D, sizeof(double));
    if(pos_f == NULL || neg_f == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    tree->computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);
    for(int n = 0; n < N; n++) tree->computeNonEdgeForces(n, theta, neg_f + n * D, &sum_Q);

    // Compute final t-SNE gradient
    for(int i = 0; i < N * D; i++) {
        dC[i] = pos_f[i] - (neg_f[i] / sum_Q);
    }
    free(pos_f);
    free(neg_f);
    delete tree;
}

// Compute gradient of the t-SNE cost function (exact)
void TSNE::computeExactGradient(double* P, double* &Y, int N, int D, double* &dC) {

    // Make sure the current gradient contains zeros
    for(int i = 0; i < N * D; i++) dC[i] = 0.0;

    // Compute the squared Euclidean distance matrix
    double* DD = (double*) malloc(N * N * sizeof(double));
    if(DD == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    computeSquaredEuclideanDistance(Y, N, D, DD);

    // Compute Q-matrix and normalization sum
    double* Q    = (double*) malloc(N * N * sizeof(double));
    if(Q == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    double sum_Q = .0;
    int nN = 0;
    for(int n = 0; n < N; n++) {
        for(int m = 0; m < N; m++) {
            if(n != m) {
                Q[nN + m] = 1 / (1 + DD[nN + m]);
                sum_Q += Q[nN + m];
            }
        }
        nN += N;
    }

    // Perform the computation of the gradient
    nN = 0;
    int nD = 0;
    for(int n = 0; n < N; n++) {
        int mD = 0;
        for(int m = 0; m < N; m++) {
            if(n != m) {
                double mult = (P[nN + m] - (Q[nN + m] / sum_Q)) * Q[nN + m];
                for(int d = 0; d < D; d++) {
                    dC[nD + d] += (Y[nD + d] - Y[mD + d]) * mult;
                }
            }
            mD += D;
        }
        nN += N;
        nD += D;
    }

    // Free memory
    free(DD); DD = NULL;
    free(Q);  Q  = NULL;
}


// Evaluate t-SNE cost function (exactly)
double TSNE::evaluateError(double* P, double* &Y, int N, int D) {

    // Compute the squared Euclidean distance matrix
    double* DD = (double*) malloc(N * N * sizeof(double));
    double* Q = (double*) malloc(N * N * sizeof(double));
    if(DD == NULL || Q == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    computeSquaredEuclideanDistance(Y, N, D, DD);

    // Compute Q-matrix and normalization sum
    int nN = 0;
    double sum_Q = DBL_MIN;
    for(int n = 0; n < N; n++) {
        for(int m = 0; m < N; m++) {
            if(n != m) {
                Q[nN + m] = 1 / (1 + DD[nN + m]);
                sum_Q += Q[nN + m];
            }
            else Q[nN + m] = DBL_MIN;
        }
        nN += N;
    }
    for(int i = 0; i < N * N; i++) Q[i] /= sum_Q;

    // Sum t-SNE error
    double C = .0;
    for(int n = 0; n < N * N; n++) {
        C += P[n] * log((P[n] + FLT_MIN) / (Q[n] + FLT_MIN));
    }

    // Clean up memory
    free(DD);
    free(Q);
    return C;
}

// Evaluate t-SNE cost function (approximately)
double TSNE::evaluateError(unsigned int* row_P, unsigned int* col_P, double* val_P, double* Y, int N, int D, double theta)
{

    // Get estimate of normalization term
    SPTree* tree = new SPTree(D, Y, N);
    double* buff = (double*) calloc(D, sizeof(double));
    double sum_Q = .0;
    for(int n = 0; n < N; n++) tree->computeNonEdgeForces(n, theta, buff, &sum_Q);

    // Loop over all edges to compute t-SNE error
    int ind1, ind2;
    double C = .0, Q;
    for(int n = 0; n < N; n++) {
        ind1 = n * D;
        for(int i = row_P[n]; i < row_P[n + 1]; i++) {
            Q = .0;
            ind2 = col_P[i] * D;
            for(int d = 0; d < D; d++) buff[d]  = Y[ind1 + d];
            for(int d = 0; d < D; d++) buff[d] -= Y[ind2 + d];
            for(int d = 0; d < D; d++) Q += buff[d] * buff[d];
            Q = (1.0 / (1.0 + Q)) / sum_Q;
            C += val_P[i] * log((val_P[i] + FLT_MIN) / (Q + FLT_MIN));
        }
    }

    // Clean up memory
    free(buff);
    delete tree;
    return C;
}


// Compute input similarities with a fixed perplexity
// input: so I guess X, N, D is the original data in which N is the number of vectors which are D dimensional, and X i the array storing the data.
// input: perplexity is from the original paper, it's the cost function parameter
// output: P is the matrix of perplexity, this optimize the perplexity 
void TSNE::computeGaussianPerplexity(double* X, int N, int D, double* P, double perplexity) {

    // Compute the squared Euclidean distance matrix
    double* DD = (double*) malloc(N * N * sizeof(double));
    if(DD == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    computeSquaredEuclideanDistance(X, N, D, DD); // Q: that's I have already, I need to transform the distance matrix to square distance matrix only inform of D, DD;

    // Compute the Gaussian kernel row by row
    int nN = 0;
    for(int n = 0; n < N; n++) {

        // Initialize some variables
        bool found = false;
        double beta = 1.0; // 1/2sigma^2 = 1.0 for all i, j?? No, that's why they have a loop to find a good perplexity
        double min_beta = -DBL_MAX;
        double max_beta =  DBL_MAX;
        double tol = 1e-5; // tolerate? to help to figure it out what is a good perplexity
        double sum_P;

        // Iterate until we found a good perplexity
        int iter = 0;
        while(!found && iter < 200) {

            // Compute Gaussian kernel row
            for(int m = 0; m < N; m++) P[nN + m] = exp(-beta * DD[nN + m]); // Q: P(n, m) = exp( -beta * DD(n,m));
            P[nN + n] = DBL_MIN;

            // Compute entropy of current row
            sum_P = DBL_MIN; // Can DBL_MIN be added to another number? Yes, it is considered as 0 when being added to another number
            for(int m = 0; m < N; m++) sum_P += P[nN + m];
            double H = 0.0; // What is H? H = beta * sum( D(n, m) * P(n, m) );
            for(int m = 0; m < N; m++) H += beta * (DD[nN + m] * P[nN + m]);
            H = (H / sum_P) + log(sum_P); 

            // Evaluate whether the entropy is within the tolerance level
            double Hdiff = H - log(perplexity);
            if(Hdiff < tol && -Hdiff < tol) { // A good perplexity needs to satisfy this condition
                found = true;
            }
            else {
                if(Hdiff > 0) {
                    min_beta = beta;
                    if(max_beta == DBL_MAX || max_beta == -DBL_MAX)
                        beta *= 2.0;
                    else
                        beta = (beta + max_beta) / 2.0;
                }
                else {
                    max_beta = beta;
                    if(min_beta == -DBL_MAX || min_beta == DBL_MAX)
                        beta /= 2.0;
                    else
                        beta = (beta + min_beta) / 2.0;
                }
            }

            // Update iteration counter
            iter++;
        }

        // Row normalize P
        for(int m = 0; m < N; m++) P[nN + m] /= sum_P;
        nN += N;
    }

    // Clean up memory
    free(DD); DD = NULL;
}

// Q: compute gaussian perplexity given the square distance matrix;
// input: N is the number of elements which is also the size of matrix DD;
// input: double* DD is the matrix of square distance among elements.
// input: perplexity is the user choice of perplexity
// output: the similarity matrix P for t-SNE in exact case which is already allocated memory before which has N*N*sizeof(double) which is same size with DD;
void TSNE::computeGaussianPerplexity(double* DD, int N, double* &P, double perplexity)
{
    // Compute the Gaussian kernel row by row
    int nN = 0;
    for(int n = 0; n < N; n++) {

        // Initialize some variables
        bool found = false;
        double beta = 1.0; // 1/2sigma^2 = 1.0 for all i, j?? No, that's why they have a loop to find a good sigma
        double min_beta = -DBL_MAX;
        double max_beta =  DBL_MAX;
        double tol = 1e-5; // tolerate? to help to figure it out what is a good perplexity
        double sum_P;

        // Iterate until we found a good perplexity
        int iter = 0;
        while(!found && iter < 200) {

            // Compute Gaussian kernel row
            for(int m = 0; m < N; m++) P[nN + m] = exp(-beta * DD[nN + m]); // Q: P(n, m) = exp( -beta * DD(n,m));
            P[nN + n] = DBL_MIN;

            // Compute entropy of current row
            sum_P = DBL_MIN; // Can DBL_MIN be added to another number? Yes, it is considered as 0 when being added to another number
            for(int m = 0; m < N; m++) sum_P += P[nN + m];
            double H = 0.0; // What is H? H = beta * sum( D(n, m) * P(n, m) );
            for(int m = 0; m < N; m++) H += beta * (DD[nN + m] * P[nN + m]);
            H = (H / sum_P) + log(sum_P); 

            // Evaluate whether the entropy is within the tolerance level
            double Hdiff = H - log(perplexity);
            if(Hdiff < tol && -Hdiff < tol) { // A good perplexity needs to satisfy this condition
                found = true;
            }
            else {
                if(Hdiff > 0) {
                    min_beta = beta;
                    if(max_beta == DBL_MAX || max_beta == -DBL_MAX)
                        beta *= 2.0;
                    else
                        beta = (beta + max_beta) / 2.0;
                }
                else {
                    max_beta = beta;
                    if(min_beta == -DBL_MAX || min_beta == DBL_MAX)
                        beta /= 2.0;
                    else
                        beta = (beta + min_beta) / 2.0;
                }
            }

            // Update iteration counter
            iter++;
        }

        // Row normalize P
        for(int m = 0; m < N; m++) P[nN + m] /= sum_P;
        nN += N;
    }
}


// Compute input similarities with a fixed perplexity using ball trees (this function allocates memory another function should free)
// Q: it seems to be a little bit more complicated over here for the sparse matrix;
void TSNE::computeGaussianPerplexity(double* X, int N, int D, unsigned int** _row_P, unsigned int** _col_P, double** _val_P, double perplexity, int K) {

    if(perplexity > K) printf("Perplexity should be lower than K!\n");

    // Allocate the memory we need
    *_row_P = (unsigned int*)    malloc((N + 1) * sizeof(unsigned int));
    *_col_P = (unsigned int*)    calloc(N * K, sizeof(unsigned int));
    *_val_P = (double*) calloc(N * K, sizeof(double));
    if(*_row_P == NULL || *_col_P == NULL || *_val_P == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    unsigned int* row_P = *_row_P;
    unsigned int* col_P = *_col_P; 
    double* val_P = *_val_P;
    double* cur_P = (double*) malloc((N - 1) * sizeof(double));
    if(cur_P == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    row_P[0] = 0;
    for(int n = 0; n < N; n++) row_P[n + 1] = row_P[n] + (unsigned int) K; // Q: initialize 0, K, 2K, 3K, ...,(N-1)K

    // Build ball tree on data set
    VpTree<DataPoint, euclidean_distance>* tree = new VpTree<DataPoint, euclidean_distance>();
    vector<DataPoint> obj_X(N, DataPoint(D, -1, X));
    for(int n = 0; n < N; n++) obj_X[n] = DataPoint(D, n, X + n * D);
    tree->create(obj_X);

    // Loop over all points to find nearest neighbors
    printf("Building tree...\n");
    vector<DataPoint> indices;
    vector<double> distances;
    for(int n = 0; n < N; n++) {

        if(n % 10000 == 0) printf(" - point %d of %d\n", n, N);

        // Find nearest neighbors
        indices.clear();
        distances.clear();
        tree->search(obj_X[n], K + 1, &indices, &distances);

        // Initialize some variables for binary search
        bool found = false;
        double beta = 1.0;
        double min_beta = -DBL_MAX;
        double max_beta =  DBL_MAX;
        double tol = 1e-5;

        // Iterate until we found a good perplexity
        int iter = 0; double sum_P;
        while(!found && iter < 200) {

            // Compute Gaussian kernel row
            for(int m = 0; m < K; m++) cur_P[m] = exp(-beta * distances[m + 1] * distances[m + 1]);

            // Compute entropy of current row
            sum_P = DBL_MIN;
            for(int m = 0; m < K; m++) sum_P += cur_P[m];
            double H = .0;
            for(int m = 0; m < K; m++) H += beta * (distances[m + 1] * distances[m + 1] * cur_P[m]);
            H = (H / sum_P) + log(sum_P);

            // Evaluate whether the entropy is within the tolerance level
            double Hdiff = H - log(perplexity);
            if(Hdiff < tol && -Hdiff < tol) {
                found = true;
            }
            else {
                if(Hdiff > 0) {
                    min_beta = beta;
                    if(max_beta == DBL_MAX || max_beta == -DBL_MAX)
                        beta *= 2.0;
                    else
                        beta = (beta + max_beta) / 2.0;
                }
                else {
                    max_beta = beta;
                    if(min_beta == -DBL_MAX || min_beta == DBL_MAX)
                        beta /= 2.0;
                    else
                        beta = (beta + min_beta) / 2.0;
                }
            }

            // Update iteration counter
            iter++;
        }

        // Row-normalize current row of P and store in matrix
        for(unsigned int m = 0; m < K; m++) cur_P[m] /= sum_P;
        for(unsigned int m = 0; m < K; m++) {
            col_P[row_P[n] + m] = (unsigned int) indices[m + 1].index();
            val_P[row_P[n] + m] = cur_P[m];
        }
    }

    // Clean up memory
    obj_X.clear();
    free(cur_P);
    delete tree;
}


// Symmetrizes a sparse matrix
// input, output: **_row_P could be a matrix, no, it's a pointer to *row_P is an array storing the row index of a sparse matrix
// input, output: **_col_P could be a matrix, no,.........................
// input, output: **_val_P could be a matrix, no, ...................................
// input: N is a number of input vector, no, N could be the number of non-zero elements
// input: there is no information about the dimension
void TSNE::symmetrizeMatrix(unsigned int** _row_P, unsigned int** _col_P, double** _val_P, int N) {

    // Get sparse matrix
    unsigned int* row_P = *_row_P; // row_P is a pointer which points to the content of _row_P which is a pointer; _row_P is a pointer to a pointer;
    unsigned int* col_P = *_col_P; // col_P points to a pointer which is the content of _col_P;
    double* val_P = *_val_P; // val_P points to a pointer which is the content pointed by _val_P; // Okay, now  they get the sparse matrix;

    //Q: a sparse matrix can be represented using three arrays, one for row, one for column, and one for value;
    //Q: the first row store number of rows, number of column, and number of non-zero elements.

    // Count number of elements and row counts of symmetric matrix
    int* row_counts = (int*) calloc(N, sizeof(int)); // What is it? Why is it a vector of N dimensions?
    if(row_counts == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    for(int n = 0; n < N; n++) {
        for(int i = row_P[n]; i < row_P[n + 1]; i++) { // that means row_P[n], row_P[n + 1] are integer ; in case row_P[n + 1] = row_P[n] this will not excecuted;
            //numbers or row_P is just a normal pointer; but how are you sure that row_P[n] < row_P[n + 1], they are!!! Row can be sorted ascendingly 

            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) { // Explain it to me? What are _row_P, _col_P,  _val_P????=> posibly it's the way they store a sparse matrix I found the answer
                if(col_P[m] == n) present = true;
            }
            if(present) row_counts[n]++;
            else {
                row_counts[n]++;
                row_counts[col_P[i]]++;
            }
        }
    }
    int no_elem = 0;
    for(int n = 0; n < N; n++) no_elem += row_counts[n];

    // Allocate memory for symmetrized matrix
    unsigned int* sym_row_P = (unsigned int*) malloc((N + 1) * sizeof(unsigned int));
    unsigned int* sym_col_P = (unsigned int*) malloc(no_elem * sizeof(unsigned int));
    double* sym_val_P = (double*) malloc(no_elem * sizeof(double));
    if(sym_row_P == NULL || sym_col_P == NULL || sym_val_P == NULL) { printf("Memory allocation failed!\n"); exit(1); }

    // Construct new row indices for symmetric matrix
    sym_row_P[0] = 0;
    for(int n = 0; n < N; n++) sym_row_P[n + 1] = sym_row_P[n] + (unsigned int) row_counts[n];

    // Fill the result matrix
    int* offset = (int*) calloc(N, sizeof(int));
    if(offset == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    for(int n = 0; n < N; n++) {
        for(unsigned int i = row_P[n]; i < row_P[n + 1]; i++) {                                  // considering element(n, col_P[i])

            // Check whether element (col_P[i], n) is present
            bool present = false;
            for(unsigned int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++) {
                if(col_P[m] == n) {
                    present = true;
                    if(n <= col_P[i]) {                                                 // make sure we do not add elements twice
                        sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                        sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                        sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i] + val_P[m];
                        sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i] + val_P[m];
                    }
                }
            }

            // If (col_P[i], n) is not present, there is no addition involved
            if(!present) {
                sym_col_P[sym_row_P[n]        + offset[n]]        = col_P[i];
                sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                sym_val_P[sym_row_P[n]        + offset[n]]        = val_P[i];
                sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = val_P[i];
            }

            // Update offsets
            if(!present || (present && n <= col_P[i])) {
                offset[n]++;
                if(col_P[i] != n) offset[col_P[i]]++;
            }
        }
    }

    // Divide the result by two
    for(int i = 0; i < no_elem; i++) sym_val_P[i] /= 2.0;

    // Return symmetrized matrices
    free(*_row_P); *_row_P = sym_row_P;
    free(*_col_P); *_col_P = sym_col_P;
    free(*_val_P); *_val_P = sym_val_P;

    // Free up some memery
    free(offset); offset = NULL;
    free(row_counts); row_counts  = NULL;
}

// Compute squared Euclidean distance matrix
// The input X, N, D is The original data array X which forms N vectors in D dimensional space;
// The output is the square distance matrix DD which compute the square distance matrix among elements in the original data;
// The n x n = 0, the diagonal line
void TSNE::computeSquaredEuclideanDistance(double* X, int N, int D, double* &DD) {
    const double* XnD = X; // XnD is X now
    for(int n = 0; n < N; ++n, XnD += D) { // for each n, XnD will move to new element, that I understood and it's clever, I learn a lot from them about the way to store the vector in one dimension array
        const double* XmD = XnD + D; // Q: XmD is the predesessor vector (not successor) of XnD 
        double* curr_elem = &DD[n*N + n]; // Q: current element = DD at n x n entry; DD is a N x N matrix
        *curr_elem = 0.0;// Q: and now they modified  DD at n x n to be 0; But that's good to know that at least they're doing something in this function;w
        double* curr_elem_sym = curr_elem + N; // current element sym = current element + N, the next row which is not the diagonal line any more
        for(int m = n + 1; m < N; ++m, XmD+=D, curr_elem_sym+=N) { // a loop from n + 1 to the end, the current element sym increases accordingly, and it's the symetric element of;
            *(++curr_elem) = 0.0; // make the current element = 0 after increasing it to 1
            for(int d = 0; d < D; ++d) {
                *curr_elem += (XnD[d] - XmD[d]) * (XnD[d] - XmD[d]); // the current element  get square Euclidean distance; DD = distance * distance matrix;
            } // Done computing D[n, m]
            *curr_elem_sym = *curr_elem; // The current element symetric will be equal to the current element;
        }
    }
}


// Makes data zero-mean
// Q: the input array X include N*D elements which stand for N element in D dimensional spaces;
// Q: the function calculates the mean vector of the data and all the vector then get substracted by the mean so that the cloud point has zero-mean;
// Q: the function modifies the data X;
void TSNE::zeroMean(double* &X, int N, int D) {

    // Compute data mean
    double* mean = (double*) calloc(D, sizeof(double));
    if(mean == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    int nD = 0;
    for(int n = 0; n < N; n++) {
        for(int d = 0; d < D; d++) {
            mean[d] += X[nD + d];
        }
        nD += D;
    }
    for(int d = 0; d < D; d++) {
        mean[d] /= (double) N;
    }

    // Subtract data mean
    nD = 0;
    for(int n = 0; n < N; n++) {
        for(int d = 0; d < D; d++) {
            X[nD + d] -= mean[d];
        }
        nD += D;
    }
    free(mean); mean = NULL;
}


// Generates a Gaussian random number
double TSNE::randn() {
    double x, y, radius;
    do {
        x = 2 * (rand() / ((double) RAND_MAX + 1)) - 1; // rand()/ (RAND_MAX + 1)  in [0, 1) then x will be in [-1, 1)
        y = 2 * (rand() / ((double) RAND_MAX + 1)) - 1; // y will be in [-1, 1)
        radius = (x * x) + (y * y); // radius is square of distance of the point in the square [-1, 1) x [-1, 1);
    } while((radius >= 1.0) || (radius == 0.0)); // Do it until get a point which is different from 0 and inside the unit circle
    radius = sqrt(-2 * log(radius) / radius); // Q: 
    x *= radius; //Q: x = x*sqrt(-2 * log(radius)/radius);
    y *= radius;
    return x;
} // Q: okay, assume it works at this point 7_9_2017

// Function that loads data from a t-SNE file
// Note: this function does a malloc that should be freed elsewhere
bool TSNE::load_data(double** data, int* n, int* d, int* no_dims, double* theta, double* perplexity, int* rand_seed, int* max_iter) {

    // Open file, read first 2 integers, allocate memory, and read the data
    FILE *h;
    if((h = fopen("data.dat", "r+b")) == NULL) {
        printf("Error: could not open data file.\n");
        return false;
    }
    fread(n, sizeof(int), 1, h);											// number of datapoints
    fread(d, sizeof(int), 1, h);											// original dimensionality
    fread(theta, sizeof(double), 1, h);										// gradient accuracy
    fread(perplexity, sizeof(double), 1, h);								// perplexity
    fread(no_dims, sizeof(int), 1, h);                                      // output dimensionality
    fread(max_iter, sizeof(int),1,h);                                       // maximum number of iterations
    *data = (double*) malloc(*d * *n * sizeof(double));
    if(*data == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    fread(*data, sizeof(double), *n * *d, h);                               // the data
    if(!feof(h)) fread(rand_seed, sizeof(int), 1, h);                       // random seed
    fclose(h);
    printf("Read the %i x %i data matrix successfully!\n", *n, *d);
    return true;
}

// Function that saves map to a t-SNE file
void TSNE::save_data(double* data, int* landmarks, double* costs, int n, int d) {

    // Open file, write first 2 integers and then the data
    FILE *h;
    if((h = fopen("result.dat", "w+b")) == NULL) {
        printf("Error: could not open data file.\n");
        return;
    }
    fwrite(&n, sizeof(int), 1, h);
    fwrite(&d, sizeof(int), 1, h);
    fwrite(data, sizeof(double), n * d, h);
    fwrite(landmarks, sizeof(int), n, h);
    fwrite(costs, sizeof(double), n, h);
    fclose(h);
    printf("Wrote the %i x %i data matrix successfully!\n", n, d);
}


void TSNE::Q_normalize_tSNE(double* &_my_tSNE, int _my_no_entries) // modified this to be able to fit into [0, 1];
{
    double _minX = DBL_MAX; 
    double _maxX = -DBL_MAX;
    double _minY = DBL_MAX;
    double _maxY = -DBL_MAX;
    for(int i = 0; i < _my_no_entries; i ++)
    {
        if(_minX > _my_tSNE[2*i]) _minX = _my_tSNE[2*i];
        if(_maxX < _my_tSNE[2*i]) _maxX = _my_tSNE[2*i];
        if(_minY > _my_tSNE[2*i+1]) _minY = _my_tSNE[2*i+1];
        if(_maxY < _my_tSNE[2*i+1]) _maxY = _my_tSNE[2*i+1];
    }

    if((_maxX != _minX) && (_maxY != _minY))
    {
        double _diffX = 1.0/(_maxX - _minX);
        double _diffY = 1.0/(_maxY - _minY);

        for(int i = 0; i < _my_no_entries; i ++)
        {
            _my_tSNE[2*i] = (_my_tSNE[2*i] - _minX)*_diffX;
            _my_tSNE[2*i + 1] = (_my_tSNE[2*i + 1] - _minY)*_diffY;
        }
    }// Done normalizing;
}

// get 51 co-activation matrix and compute t-SNE
// as p[0] is random, but p[t] is the initializing of p[t+1] with t >= 1;
/*
void TSNE::Q_t_SNE_MultipleFile(int numberOfFile, int numberOfRow)
{
    // Init_p_0 = random();
    double* Y;
    Y = (double*) malloc(numberOfRow * 2 * sizeof(double)); // Q's code
    if(Y == NULL) { printf("Memory allocation failed! \n"); exit(1); }; // Q's code
    // Initialize solution (randomly)
    vector<point> Input;
    for(int i = 0; i < numberOfRow * 2; i++) 
    {
        Y[i] = randn() * .0001;
        //point PointI;
        //Input.push_back(PointI);
    }
    Q_transform_tSNE2VectorOfPoint(Y, numberOfRow, Input);
    //cout << "The size of the input is " << (int) Input.size() << endl;
    vector<vector<float>> DisMat;
    // Utils::Read first_file;
    Utils::read_DistanceMat("DisMat_0", DisMat, numberOfRow);
    vector<point> Output;
    // Q_run(first_file, Init_p_0, &output);
    Q_run_withCo_activation(DisMat, Input, Output);
    // Q_write2File(output);
    Q_write2File_tSNE(Output, "tSNE_0");
    // Init_p = output;
    Input.clear();
    Input = Output;
    DisMat.clear();
    Output.clear();
    for(int i = 1; i < numberOfFile; i ++)
    {
        ostringstream ReadFilename;
        ReadFilename << "DisMat_" << i;
        // Read file_I;
        Utils::read_DistanceMat(ReadFilename.str(), DisMat, numberOfRow);
        // Q_run(file_I, Init_p, &outputI);
        Q_run_withCo_activation(DisMat, Input, Output);
        // Q_write2File(outputI);
        ostringstream WriteFilename;
        WriteFilename << "tSNE_" << i;
        Q_write2File_tSNE(Output, WriteFilename.str());
        // Init_p = outputI;
        Input.clear(); // Q: Don't know if I have to do that
        Input = Output;
        DisMat.clear();
        Output.clear();
    }
}

*/
// same functionality as the above function but compute t-SNE 1D instead
// Y is not yet allocated
/*
void TSNE::Q_t_SNE_MultipleFile_1D(int numberOfFile, int NumberOfRow, double perplexity) // lala
{
    // Init_p_0 = random();
    vector<double> Input;
    //Y = (double*) malloc(numberOfRow * sizeof(double)); // Q's code
    //if(Y == NULL) { printf("Memory allocation failed! \n"); exit(1); }; // Q's code
    // Initialize solution (randomly)
    for(int i = 0; i < NumberOfRow; i++) 
    {
        Input.push_back(randn() * .0001);
        //point PointI;
        //Input.push_back(PointI);
    }

    vector<vector<float>> DisMat;
    // Utils::Read first_file;
    Utils::read_DistanceMat("DisMat_0", DisMat, NumberOfRow);
    vector<double> Output;
    // Q_run(first_file, Init_p_0, &output);
    Q_run_withCo_activation_1D(DisMat, Input, Output, perplexity ); // I have to re-write this function => done hehe I'm smart:P
    // Q_write2File(output);
    Q_write2File_tSNE_1D(Output, "1D_tSNE_0"); // I have to re-write this function and the function to read 1D tsne from file, and the function to normalize after reading from file;
    // Init_p = output;
    Input.clear();
    Input = Output;
    DisMat.clear();
    Output.clear();
    for(int i = 1; i < numberOfFile; i ++)
    {
        ostringstream ReadFilename;
        ReadFilename << "DisMat_" << i;
        // Read file_I;
        Utils::read_DistanceMat(ReadFilename.str(), DisMat, NumberOfRow);
        Q_run_withCo_activation_1D(DisMat, Input, Output, perplexity);
        // Q_write2File(outputI);
        ostringstream WriteFilename;
        WriteFilename << "1D_tSNE_" << i;
        Q_write2File_tSNE_1D(Output, WriteFilename.str());
        // Init_p = outputI;
        Input.clear(); // Q: Don't know if I have to do that
        Input = Output;
        DisMat.clear();
        Output.clear();
    }
}
*/
// The input as well as output Y is initialized already;

void TSNE::Q_run_withY_allocated(double* DD, int N, double* &Y, int no_dims, double perplexity, int rand_seed,
        bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter)
{
    // Set random seed
    if (skip_random_init != true) {
        if(rand_seed >= 0) {
            printf("Using random seed: %d\n", rand_seed);
            srand((unsigned int) rand_seed);
        } else {
            printf("Using current time as random seed...\n");
            srand(time(NULL));
        }
    }// Q: alles gut but what is the randomness used for? for initialize Y in lower dimensional space

    if(N - 1 < 3 * perplexity) { printf("Perplexity too large for the number of data points!\n"); exit(1); } // Q: what is N? and What is perplexity. why do you this condition

    // Set learning parameters
    float total_time = .0; // What is this?
    clock_t start, end; // Q: What are they?
    double momentum = .5, final_momentum = .8; // What are they? Why do you use them?
    double eta = 200.0; // What is this? Why do you use it?

    // Allocate some memory
    double* dY    = (double*) malloc(N * no_dims * sizeof(double)); // What is this? It's an array include N*no_dims of double, no_dims could be the number of lower dimensional space
    double* uY    = (double*) malloc(N * no_dims * sizeof(double)); // Then I guess they could be to store the original data points into lower dimension
    double* gains = (double*) malloc(N * no_dims * sizeof(double));
    if(dY == NULL || uY == NULL || gains == NULL) { printf("Memory allocation failed!\n"); exit(1); } // Q: just to check if allocating memory is successful!
    for(int i = 0; i < N * no_dims; i++)    uY[i] =  .0; // Q: initializing uY = 0.0;
    for(int i = 0; i < N * no_dims; i++) gains[i] = 1.0; // Q: intializing gains = 1.0;

    // Normalize input data (to prevent numerical problems)
    printf("Computing input similarities...\n");
    start = clock(); // Q: start the clock here to count time

    // Compute input similarities for exact t-SNE
    double* P;
    // Compute similarities
    P = (double*) malloc(N * N * sizeof(double)); // Q: it seems like P is a matrix having size N x N 
    if(P == NULL) { printf("Memory allocation failed!\n"); exit(1); } // Q: check if allocating memory successfully
    computeGaussianPerplexity(DD, N,  P, perplexity); // Q: it's my version of computing gaussian perplexity

    // Symmetrize input similarities
    printf("Symmetrizing...\n");
    int nN = 0;
    for(int n = 0; n < N; n++) {
        int mN = (n + 1) * N;
        for(int m = n + 1; m < N; m++) {
            P[nN + m] += P[mN + n];
            P[mN + n]  = P[nN + m];
            mN += N;
        }
        nN += N;
    }
    double sum_P = .0;
    for(int i = 0; i < N * N; i++) sum_P += P[i];
    for(int i = 0; i < N * N; i++) P[i] /= sum_P;

    end = clock(); // Q: End clock here

    // Lie about the P-values
    for(int i = 0; i < N * N; i++) P[i] *= 12.0; 

    //Y = (double*) malloc(N * no_dims * sizeof(double)); // Q's code
    //if(Y == NULL) { printf("Memory allocation failed! \n"); exit(1); }; // Q's code


    // Initialize solution (randomly)
    //if (skip_random_init != true) {
    //for(int i = 0; i < N * no_dims; i++) Y[i] = randn() * .0001;
    //} // Q: There we go the initialize solution which I should choose false, because I want it to be random;

    // Perform main training loop
    printf("Input similarities computed in %4.2f seconds!\nLearning embedding...\n", (float) (end - start) / CLOCKS_PER_SEC);

    start = clock();
    for(int iter = 0; iter < max_iter; iter++) {

        // Compute (approximate) gradient
        computeExactGradient(P, Y, N, no_dims, dY); // it has Y and no_dims

        // Update gains
        for(int i = 0; i < N * no_dims; i++) gains[i] = (sign(dY[i]) != sign(uY[i])) ? (gains[i] + .2) : (gains[i] * .8);
        for(int i = 0; i < N * no_dims; i++) if(gains[i] < .01) gains[i] = .01;

        // Perform gradient update (with momentum and gains)
        for(int i = 0; i < N * no_dims; i++) uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
        for(int i = 0; i < N * no_dims; i++)  Y[i] = Y[i] + uY[i];

        // Make solution zero-mean
        zeroMean(Y, N, no_dims);

        // Stop lying about the P-values after a while, and switch momentum
        if(iter == stop_lying_iter) {
            for(int i = 0; i < N * N; i++) P[i] /= 12.0; 
        }
        if(iter == mom_switch_iter) momentum = final_momentum;

        // Print out progress
        if (iter > 0 && (iter % 50 == 0 || iter == max_iter - 1)) {
            end = clock();
            double C = .0;
            C = evaluateError(P, Y, N, no_dims);
            if(iter == 0)
                printf("Iteration %d: error is %f\n", iter + 1, C);
            else {
                total_time += (float) (end - start) / CLOCKS_PER_SEC;
                //  printf("Iteration %d: error is %f (50 iterations in %4.2f seconds)\n", iter, C, (float) (end - start) / CLOCKS_PER_SEC);
            }
            start = clock();
        }
    }
    end = clock(); total_time += (float) (end - start) / CLOCKS_PER_SEC;

    // Clean up memory
    free(dY);
    free(uY);
    free(gains);
    free(P);
    printf("Fitting performed in %4.2f seconds.\n", total_time);
}



// What is Input_Allocated???
// input: Co_mat is the co-activation matrix
// input: Input_Allocated?? What is this??? Ah, it's the output of the previous tsne, and consider to be the input
// ot this tsne run to get the output;
/*
void TSNE::Q_run_withCo_activation(const vector<vector<float>> &Co_mat, const vector<point> &Input_Allocated, vector<point> &output)
{
    double* Y;
    Q_transfer_vectorPoint2Pointer(Input_Allocated,  Y);
    int _size = (int)Input_Allocated.size();
    double* DD;
    // Compute the squared Euclidean distance matrix
    DD = (double*) malloc(_size * _size * sizeof(double));
    if(DD == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    Q_transfer_DistMat_To_SquareDistMat_Pointer(Co_mat, DD, _size);
    Q_run_withY_allocated(DD, _size, Y, 2, 2.0, -1, false, 100000, 2500, 2500);//, int rand_seed,
    //bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter)
    //Q_normalize_tSNE(Y, _size);
    Q_transform_tSNE2VectorOfPoint(Y, _size, output);
    free(DD);
    free(Y);
}

/// same as function above but for 1D;
void TSNE::Q_run_withCo_activation_1D(const std::vector<std::vector<float>> &DisMat, const std::vector<double> &Input, std::vector<double> &Output, double perplexity)
{
    double* Y;
    int _size = (int) Input.size();
    Y = (double*) malloc(_size*sizeof(double));
    if(Y == NULL) { printf("Memory allocation failed!\n"); exit(1);}
    for(int i = 0; i < _size; i ++ )
    {
        Y[i] = Input[i];
    }
    double* DD;
    // Compute the squared Euclidean distance matrix
    DD = (double*) malloc(_size * _size * sizeof(double));
    if(DD == NULL) { printf("Memory allocation failed!\n"); exit(1); }
    Q_transfer_DistMat_To_SquareDistMat_Pointer(DisMat, DD, _size);
    Q_run_withY_allocated(DD, _size, Y, 1, perplexity, -1, false, 100000, 2500, 2500);//, int rand_seed,
    //bool skip_random_init, int max_iter, int stop_lying_iter, int mom_switch_iter)
    if(!Output.empty()) Output.clear();
    for(int i = 0; i < _size; i ++)
    {
        Output.push_back(Y[i]);
    }
    free(Y);
}
void TSNE::Q_write2File_tSNE_1D(const std::vector<double> &Output, std::string filename) const // to write to filename the tsne_1D
{
    ofstream dev;
    dev.open(filename, ios::out);

    if(!dev.good())
    {
        cerr << "Cannot open the file." << endl;
        exit( 1 );
    }

    int Length_Vector = (int) Output.size();// (int) MDS.size();

    for(int i = 0; i < Length_Vector; i ++ )
    {
        if( i != Length_Vector - 1 )dev << Output[i] << " ";
        else if ( i == Length_Vector - 1 ) dev << Output[i] << endl;
    }
    dev.close();
    cout << "Write t-SNE to file successfully!" << endl;
}
*/

/*
void TSNE::Q_loadFile_tSNE_1D(std::vector<double> &output, std::string filename) // to read tSNE_1D from file
{
    if(!output.empty()) output.clear(); // erarse existing data.

    ifstream dev;
    dev.open( filename );

    if(!dev.good())
    {
        cerr << "Cannot open the file." << endl;
        exit( 1 );
    }

    while(dev.good())
    {
        string RawString;
        getline(dev, RawString, '\n');
        vector<string> item;
        item = Utils::tokenize(RawString, ' ');
        if(item.empty()) continue;
        for(int i = 0; i < (int)item.size(); i ++) output.push_back(Utils::stringtoDouble(item[i]));
    }

    dev.close();
    cout << "The size of output is " << (int) output.size() << endl;
    cout << "Loading file t_SNE 1D successfully." << endl;
}
*/

void TSNE::Q_normalize_tSNE_1D(vector<double> &input_tSNE) // modified this to fit into [0, 1];
{
    int _size = (int) input_tSNE.size();
    double _min = input_tSNE[0];
    double _max = input_tSNE[0];
    for(int i = 0; i < _size; i ++)
    {
        if( input_tSNE[i] < _min ) _min = input_tSNE[i];
        if( input_tSNE[i] > _max ) _max = input_tSNE[i];
    }
    if(_max != _min)
    {
        double _diff = 1.0/(_max - _min);
        for(int i = 0; i < _size; i ++)
        {
            input_tSNE[i] = (input_tSNE[i] - _min)*_diff;
        }
    }
} // done normalizing 1D;
