/*
    BVP_LCX_jac.C 
        MEX file corresponding to BVPjac.m
        Does the evaluation of the jacobian of the BVP
        
    calling syntax:
        result = BVP_LCX_jac(lds.func,x,p,T,pars,nc,lds,p2,lds.Jacobian,lds.ActiveParams,lds.JacobianP)
*/

#include<math.h>
#include<mex.h>
#include<matrix.h>
#include<stdio.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


    /* Declarations */
	/* ------------ */

	double *x,*p, *nc;
	mxArray *thisfield;
	int ntst, ncol, nphase, *ActiveParams, ncoords, nfreep, tps;
	double *upoldp, *mesh, *wt, *wp, *pwi;

	double *dt, *wploc, T;

	double *pr;
	mwIndex *ir, *jc;
    
    mxArray *evalrhs[1000], *jacrhs[5];
	mxArray *evallhs[1], *jaclhs[1];
	double *xtmp, *xtmpn, *xtmpm;

	int filled, elementcounter, remm;
	int i,j,k,l,l2;	/* Indexation variables */

	int *range1, *range2, *range3, *range4;

	double *jac, *jacp, *icjac, *sysjac, *sysjacp, *zero, *ptmp;
	double *frhs, *frhstmp;

	double *Tcol, *freepcols, *tempmatrix;
	double tmpperiod;

	/* Initializations */
	/* --------------- */
	    
	/* Retrieve parameters. */
	x = mxGetPr(prhs[1]);
	p = mxGetPr(prhs[2]);
    nc = mxGetPr(prhs[5]);
    
    /* LDS FIELDS */
	thisfield = mxGetFieldByNumber(prhs[6],0,10);
	nphase = *(mxGetPr(thisfield));		/* Size of one point */
    thisfield = mxGetFieldByNumber(prhs[6],0,11);
    nfreep = mxGetNumberOfElements(thisfield);/* number of free parameters */
    ActiveParams = mxCalloc(nfreep,sizeof(int));
    frhstmp = mxGetPr(thisfield);
    for (i=0; i<nfreep; i++) 
        *(ActiveParams+i) = (int)(*(frhstmp+i));         
	thisfield = mxGetFieldByNumber(prhs[6],0,13);
	ntst = *(mxGetPr(thisfield));	/* Number of mesh intervals */
	thisfield = mxGetFieldByNumber(prhs[6],0,14);
	ncol = *(mxGetPr(thisfield));		/* Number of collocation points */
    thisfield = mxGetFieldByNumber(prhs[6],0,17);
    tps = *(mxGetPr(thisfield));      
    thisfield = mxGetFieldByNumber(prhs[6],0,18);
    ncoords = *(mxGetPr(thisfield));    
	thisfield = mxGetFieldByNumber(prhs[6],0,22);
	mesh = mxGetPr(thisfield);			/* Current mesh coordinates */
    thisfield = mxGetFieldByNumber(prhs[6],0,24);
    dt = mxGetPr(thisfield); /* Interval widths */   
	thisfield = mxGetFieldByNumber(prhs[6],0,25);
	upoldp = mxGetPr(thisfield);	/* Derivative of cycle at old mesh coordinates */
	thisfield = mxGetFieldByNumber(prhs[6],0,29);
	wt = mxGetPr(thisfield);		/* Weights of collocation points */
	thisfield = mxGetFieldByNumber(prhs[6],0,31);
	T = *(mxGetPr(thisfield));		/* Weights of collocation points */
	thisfield = mxGetFieldByNumber(prhs[6],0,39);
	wp = mxGetPr(thisfield);	/* Derivative weights of collocation points */
    /* Kronecker product of the derivative weights and the identity matrix */
	wploc = mxCalloc(mxGetN(thisfield)*mxGetM(thisfield),sizeof(double));
	thisfield = mxGetFieldByNumber(prhs[6],0,41);
	pwi = mxGetPr(thisfield);	/* Extension of weights */
	    
    /* Sparse matrix as returnvalue */
	plhs[0] = mxCreateSparse(ncoords+1,ncoords+1+nfreep,ncoords*ncoords,mxREAL);
	pr = mxGetPr(plhs[0]);
	ir = mxGetIr(plhs[0]);
	jc = mxGetJc(plhs[0]);
	*jc = 0;
	   
    
	/* Parameters for rhs-evaluation-call to Matlab */
	evalrhs[0] = (struct mxArray_tag*) prhs[0];
	evalrhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    zero = mxGetPr(evalrhs[1]);
    *zero = 0;
	/*evalrhs[2] = mxCreateDoubleMatrix(nphase+ActiveParams,1,mxREAL);*/
    evalrhs[2] = mxCreateDoubleMatrix(nphase,1,mxREAL);
	xtmp = mxGetPr(evalrhs[2]); 
    for (i=0; i<mxGetNumberOfElements(prhs[2]); i++) {
        evalrhs[i+3] = mxCreateDoubleMatrix(1,1,mxREAL);
        ptmp = mxGetPr(evalrhs[i+3]);   
        *ptmp = *(p+i);
    }
	/* Parameters for jacobian-call to Matlab */
    jacrhs[0] = (struct mxArray_tag*) prhs[0];
    jacrhs[1] = (struct mxArray_tag*) prhs[8];
    jacrhs[2] = evalrhs[2];
    jacrhs[3] = (struct mxArray_tag*) prhs[7];
    jacrhs[4] = mxCreateDoubleMatrix(nfreep,1,mxREAL);
    xtmpn = mxGetPr(jacrhs[4]);
    for (i=0; i<nfreep; i++)
        *(xtmpn+i) = *(ActiveParams+i);

	filled = 0;			/* Help-variable that will be used in storage-procedure */
	elementcounter = 0;	/* Counts number of elements already stored in sparse matrix */
    tmpperiod = *(mxGetPr(prhs[3]));
   
	/* Other memory allocations */
	/* ------------------------ */    
	range1 = mxCalloc(ncol+1,sizeof(int));
	range2 = mxCalloc(ncol*nphase,sizeof(int));
	range3 = mxCalloc((ncol+1)*nphase,sizeof(int));
	range4 = mxCalloc(nphase,sizeof(int));
	tempmatrix = mxCalloc(nphase*nphase*ncol,sizeof(double));
	
	jac = mxCalloc(nphase*nphase,sizeof(double));
	jacp = mxCalloc(nphase*nfreep,sizeof(double));
	frhs = mxCalloc(nphase*ncol,sizeof(double));
	icjac = mxCalloc(ncoords,sizeof(double));
	sysjac = mxCalloc(nphase*ncol*(ncol+1)*nphase,sizeof(double));
	sysjacp = mxCalloc(nphase*ncol*nfreep,sizeof(double));

	Tcol = mxCalloc((tps-1)*nphase,sizeof(double));
	freepcols = mxCalloc((tps-1)*nfreep*nphase,sizeof(double));
    

	/* Compute third component: the integral constraint */
	/* ------------------------------------------------ */
	/* Storage in sparse matrix is done later on */

	/* Define some ranges */
	for (i=0; i<(ncol+1); i++) {
		*(range1+i) = i;		
		for (j=0; j<nphase; j++) 
			*(range3+i*nphase+j) = i*nphase+j;      
	}

    
	for (i=0; i<ntst; i++) {
		/* Compute elements of third component */
		for (j=0; j<(ncol+1)*nphase; j++) {
    		*(icjac + *(range3+j)) = *(icjac + *(range3+j)) + (*(dt+i)) * (*(upoldp+(*range1)*nphase+j)) * (*(pwi+j));
        }
		/* Shift the ranges to next intervals */
		for (j=0; j<ncol+1; j++) {
			*(range1+j) = *(range1+j) + ncol;
			for (k=0; k<nphase; k++)
				*(range3+j*nphase+k) = *(range3+j*nphase+k) + ncol*nphase;
		}	
	}
    
	   
	/* Compute first component */
	/* ----------------------- */
	
	/* Define some ranges*/
	for (i=0; i<(ncol+1); i++) {
		*(range1+i) = i;	
		if (i < ncol) 
			for (j=0; j<nphase; j++) {
				*(range2+i*nphase+j) = i*nphase+j;
				*(range3+i*nphase+j) = i*nphase+j;
			}
		else 
			for (j=0; j<nphase; j++)
				*(range3+i*nphase+j) = i*nphase+j;
	}    

	/* Actual computation of component elements */
	for (i=0; i<ntst; i++) {
		
		/* Define a new range */
		for (j=0; j<nphase; j++)
			*(range4+j) = j;

		for (j=0; j<((ncol+1)*nphase)*(ncol*nphase); j++)
            *(wploc+j) = *(wp+j) / *(dt+i);

		for (j=0; j<ncol; j++) {
			/* Compute value of the polynomial in mesh point */
			for (k=0; k<nphase; k++) {
				*(xtmp+k) = 0;				
				for (l=0; l<(ncol+1); l++)
					*(xtmp+k) = *(xtmp+k) + (*(x+(*(range1+l))*nphase+k)) * (*(wt+j*(ncol+1)+l));
			}
            
			/* Call to Matlab for evaluation of rhs */
			mexCallMATLAB(1,evallhs,3+mxGetNumberOfElements(prhs[2]),evalrhs,"feval");
			frhstmp = mxGetPr(evallhs[0]);
  
			for (k=0; k<nphase; k++) 
				*(frhs+j*nphase+k) = *(frhstmp+k);
       /*     mxDestroyArray(evallhs[0]);*/
            
			/* Call to Matlab for evaluation of jacobian */
            jacrhs[1] = (struct mxArray_tag*) prhs[8];
			mexCallMATLAB(1,jaclhs,5,jacrhs,"cjac");
			frhstmp = mxGetPr(jaclhs[0]);

			/* Store jacobian */
			for (k=0; k<nphase*nphase; k++) 	
				*(jac+k) = *(frhstmp+k);
         /*   mxDestroyArray(jaclhs[0]);*/
            
			/* Call to Matlab for evaluation of jacp */
            jacrhs[1] = (struct mxArray_tag*) prhs[10];
			mexCallMATLAB(1,jaclhs,5,jacrhs,"cjacp");
			frhstmp = mxGetPr(jaclhs[0]);            
            
            for (k=0; k<nfreep; k++)
                for (l=0; l<nphase; l++)
                    *(jacp+k*nphase+l) = *(frhstmp+k*nphase+l);
         /*   mxDestroyArray(jaclhs[0]);*/

			/* temporary sysjac and sysjacp */
			for (k=0; k<nphase; k++) {
				/* sysjac stores kronecker product of jacobian and weights */
			    for (l=nphase; l<(ncol+2)*nphase; l++) {
			        l2 = floor(l/nphase)-1;
			        remm = l % nphase;
    			    *(sysjac + (l-nphase)*nphase*ncol + (*(range4+k))) = *(wt+j*(ncol+1)+l2) * (*(jac+remm*nphase+k));
    			}
            }
            for (l=0; l<nfreep; l++)
                for (k=0; k<nphase; k++)
                    *(sysjacp + l*nphase*ncol + (*(range4+k))) = *(jacp + l*nphase + k);
            
			/* Shift range4 */
			for (k=0; k<nphase; k++) 
				*(range4+k) = *(range4+k) + nphase;
		}
                        
		/* Storage in sparse return matrix */
		/* ------------------------------- */
		
		/* The columns are stored one at a time. Because of the way of storing the matrix, all elements of a column
		must be stored consecutively. Therefore, sometimes some elements will be stored in a temporary matrix. */

		/* Finish computing and store the current (ncol+1)*nphase interval-columns */ 
		for (k=0; k<(ncol+1)*nphase; k++) {
					
			/* Check to see if some elements have been computed previously and were stored temporarily */
			if  ((k<ncol*nphase) || (*(range3+k) > (tps-1)*nphase-1)) {
                
				if (filled) {
					/* Fill in previously computed non-zero elements */
					for (j=0; j<ncol*nphase; j++) {
						if (*(tempmatrix + k*ncol*nphase + j)) {
							*(pr + elementcounter) = *(tempmatrix + k*ncol*nphase + j);
							*(ir + elementcounter) = *(range2 + j) - ncol*nphase;
							elementcounter = elementcounter + 1;
							/* Clear temporary storage matrix */
							*(tempmatrix + k*ncol*nphase + j) = 0;
						}
					}
					if (k == nphase-1)
						/* Reset indicator */
						filled = 0;
				}
                
				/* Do final computations on first component-elements and store the column */
				for (j=0; j<ncol*nphase; j++) 
                    if ((*(wploc+k*(ncol*nphase)+j))-(tmpperiod)*(*(sysjac+k*nphase*ncol+j))) {
		                *(pr + elementcounter) = (*(wploc+k*(ncol*nphase)+j))-(tmpperiod)*(*(sysjac+k*nphase*ncol+j));  
		                *(ir + elementcounter) = *(range2 + j);
		                elementcounter = elementcounter + 1;
		            }				
                
				/* Fill in possible non-zero elements of second and third component */
				if (*(range3+k) < ncoords) {
					if (*(range3+k) < nphase) {
						/* first I-matrix of second component */
    					*(pr + elementcounter) = 1;
						*(ir + elementcounter) = (tps - 1) * nphase + *(range3+k);
        				elementcounter = elementcounter + 1;
    				}   
    				else {
				        if (*(range3+k) > (tps-1)*nphase-1) {
					        /* Second I-matrix of second component */
				            *(pr + elementcounter) = -1;
					        *(ir + elementcounter) = *(range3+k);
						    elementcounter = elementcounter + 1;
						}
					}
				    
				    /* Fill in previously computed element from third component */
				    *(pr + elementcounter) = *(icjac + *(range3+k));
				    *(ir + elementcounter) = ncoords;
					elementcounter = elementcounter + 1;
				}    	
				 
				/* Finish the column */
				*(jc + *(range3+k) + 1) = elementcounter;
			}
            
            else {
				/* These elements are destined for columns which will be assigned other elements later on. 
				So we store these temporarily. */                
				filled = 1;
				for (j=0; j<ncol*nphase; j++)  
                        *(tempmatrix + (k-ncol*nphase)*ncol*nphase + j) = (*(wploc+k*(ncol*nphase)+j))-(tmpperiod)*(*(sysjac+k*nphase*ncol+j));
			}
		}
        
            for (j=0; j<ncol*nphase; j++) {
                 *(Tcol + *(range2+j)) = -(*(frhs+j));
                for (k=0; k<nfreep; k++) {
                    *(freepcols + k*(tps-1)*nphase + *(range2+j)) = -(tmpperiod)*(*(sysjacp + nphase*ncol*k + j));
                }          
            }
        
		/* Shift ranges to next intervals */
		for (j=0; j<ncol+1; j++) {
			*(range1+j) = *(range1+j) + ncol;
			if (j < ncol)
				for (k=0; k<nphase; k++) {
					*(range2+j*nphase+k) = *(range2+j*nphase+k) + ncol*nphase;
					*(range3+j*nphase+k) = *(range3+j*nphase+k) + ncol*nphase;
				}	
			else
				for (k=0; k<nphase; k++)
					*(range3+j*nphase+k) = *(range3+j*nphase+k) + ncol*nphase;
		}	
	}/* end ntst*/
    	
	/* Finally, store the last 2 columns in the sparse return matrix */
    for (i=0; i<nfreep+1; i++) {
        
        for (j=0; j<(tps-1)*nphase; j++) {
            if (i == 0) {
                    /* Store the column from the period */
                if (*(Tcol + j)) {
                    *(pr + elementcounter) = *(Tcol + j);
                    *(ir + elementcounter) = j;
                    elementcounter = elementcounter+1;
                }
            }
            else {
                    /* Store the columns from the free parameters */
                if (*(freepcols + (i-1)*nphase*(tps-1) + j)) {
                    *(pr + elementcounter) = *(freepcols + (i-1)*nphase*(tps-1) + j);
                    *(ir + elementcounter) = j;
                    elementcounter = elementcounter+1;
                }
            }
        }
        
            /* Finish the column */
        *(jc + ncoords+1 + i) = elementcounter;
    }
    
	/* Free all allocated memory */    
	/* ------------------------- */

	mxFree(wploc);
    mxFree(ActiveParams);
	
	mxFree(range1);
	mxFree(range2);
	mxFree(range3);
	mxFree(range4);

	mxFree(jac);
	mxFree(jacp);
	mxFree(icjac);
	mxFree(sysjac);
	mxFree(sysjacp);
	mxFree(frhs);

	mxFree(Tcol);
	mxFree(freepcols);
	mxFree(tempmatrix);    
    
			
	return;
}
