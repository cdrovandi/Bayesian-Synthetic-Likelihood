#include <math.h>
#include <string.h>
#include <mex.h>
#include "mt19937ar.h"

// simulate from U(a,b)
double unifrnd(double a,double b){
	return ((b-a)*genrand_real3() + a);
}

// sample from categorical distribution with equal probabilities
int sample_equal(int n){ 
    double cumsum = 0.0;
    double inc;
    int i;
    double p = unifrnd(0.0, 1.0);
    inc = 1.0/(1.0*n);
    for (i = 0; i < n; i++){
        cumsum += inc;
        if (p <= cumsum){
            return(i);
        }
    }
}

// the model simulation
int simulate(int *x, int *rows, int *cols, double Pm, double Pp, int N, int Ninit, int nrow, int ncol, int sim_iters){
        int i, r, row, col, row_prop, col_prop, K, Kfix, s; 
		double u;
        int num_motility_prop = 0, num_prolif_prop = 0;
        
        K=N;
        
        // simulate over discrete time steps
        for (s=0; s<sim_iters; s++){
            // potential motility event for each cell
            Kfix = K;
            for (i=0; i<Kfix; i++){
                if (unifrnd(0.0,1.0)<Pm){ //attempt a move with some probability (inputted probability of motility)
                    num_motility_prop++;
                    r = sample_equal(K);
                    row = rows[r];
                    col = cols[r];
                    u = unifrnd(0.0,1.0); //randomly choose location
                    if (u < 0.25){
                        row_prop = row-1;
                        col_prop = col;
                    }
                    else if (u > 0.25 && u < 0.5){
                        row_prop = row+1;
                        col_prop = col;
                    }
                    else if (u > 0.5 && u < 0.75){
                        row_prop = row;
                        col_prop = col-1;
                    }else{
                        row_prop = row;
                        col_prop = col+1;
                    }

                    if (row_prop >= 0 && row_prop <= (nrow-1) && col_prop >= 0 && col_prop <= (ncol-1)){ //if chosen location is empty and within bounds then move is successful
                        if (x[row_prop + nrow*col_prop]==0){
                            // move
                            x[row_prop + nrow*col_prop]=1;
                            x[row + nrow*col] = 0;
                            rows[r] = row_prop;
                            cols[r] = col_prop;
                        }
                    }

                }


            }

            // potential proliferation event for each cell
            for (i=0; i<Kfix; i++){
                if (unifrnd(0.0,1.0)<Pp){ //attempt a proliferation event with some probability (inputted probability of proliferation)
                    num_prolif_prop++;
                    r = sample_equal(Kfix);
                    row = rows[r];
                    col = cols[r];
                    u = unifrnd(0.0,1.0);  //randomly choose location and if empty then move is successful
                    if (u < 0.25){
                        row_prop = row-1;
                        col_prop = col;
                    }
                    else if (u > 0.25 && u < 0.5){
                        row_prop = row+1;
                        col_prop = col;
                    }
                    else if (u > 0.5 && u < 0.75){
                        row_prop = row;
                        col_prop = col-1;
                    }else{
                        row_prop = row;
                        col_prop = col+1;
                    }

                    if (row_prop >= 0 && row_prop <= (nrow-1) && col_prop >= 0 && col_prop <= (ncol-1)){ //if chosen location is empty and within bounds then proliferation is successful
                        if (x[row_prop + nrow*col_prop]==0){
                            // proliferate
                            x[row_prop + nrow*col_prop]=1;
                            rows[K] = row_prop;
                            cols[K] = col_prop;
                            K++;
                            if (K >= 2*Ninit){
                                // then too many birthing of cells to store in array (this parameter will end up being rejected)
                                return(K);
                            }
                        }
                    }

                }
            }
        
        }
        return(K);
}


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) {

	int *x, *xc, *rows, *cols, *xr, *rowsr, *colsr, *colsc, *rowsc; 
    double Pm, Pp;
	int nrow, ncol, i, j, N, Ninit, K, seed, sim_iters, num_obs, k, pos = 0;
    int  ndim = 3, *dims;
    N=0;

    // read in the data
	x = (int*)mxGetData(prhs[0]);
    rows = (int*)mxGetData(prhs[1]);
    cols = (int*)mxGetData(prhs[2]);
    sim_iters = (int)mxGetScalar(prhs[5]);
    num_obs = (int)mxGetScalar(prhs[6]);
    seed = (int)mxGetScalar(prhs[7]); // random seed for simulation
	nrow = (int)mxGetM(prhs[0]);
    ncol = (int)mxGetN(prhs[0]);
    Ninit = mxGetNumberOfElements(prhs[1]);
    N = Ninit;
    
    // set the seed
    init_genrand(seed);
    
    // model parameters
    Pm = mxGetScalar(prhs[3]); Pp = mxGetScalar(prhs[4]);
    
    dims = (int*)malloc(3*sizeof(int));
    dims[0] = nrow; dims[1] = ncol; dims[2] = num_obs;
   	
	// allocate double the space required 
    rowsc = (int*)mxMalloc((Ninit*2)*sizeof(int));
    colsc = (int*)mxMalloc((Ninit*2)*sizeof(int));
    
    // output (collection of binary matrices)
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    xr  = (int*)mxGetData(plhs[0]);
	
	memcpy(rowsc, rows, Ninit*sizeof(int));
	memcpy(colsc, cols, Ninit*sizeof(int));
    
    xc = (int*)mxMalloc((nrow*ncol)*sizeof(int));
    memcpy(xc, x, (nrow*ncol)*sizeof(int));
    
    // perform simulation over num_obs time points
    for (k=0; k<num_obs; k++){
         if (N < 2*Ninit){ 
            N = simulate(xc, rowsc, colsc, Pm, Pp, N, Ninit, nrow, ncol, sim_iters);
         }
         // store result in output object
         for (i=0; i<nrow; i++){
             for (j=0; j<ncol; j++){
                 xr[pos + i + nrow*j] = xc[i + nrow*j];
             }
         }
         pos += nrow*ncol;
     }
    
    // clean memory
    free(dims);
}
