/*
    BVP_LPC_jac.c 
        MEX file corresponding to LPC_jac.m
        Does the evaluation of the jacobian of the BVP of the LPC
        
    calling syntax:
        result = BVP_LPC_jac(lds.func,x,p,T,pars,nc,lds,gds.period,p2)
*/

#include<math.h>
#include<mex.h>
#include<matrix.h>
#include<stdio.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


    /* Declarations */
	/* ------------ */

	double *x,*p,*pars,*pars2, *nc;
	mxArray *thisfield;
	int ntst, ncol, nphase, *ActiveParams, ncoords, nfreep, tps;
	double *upoldp, *mesh, *wt, *wp, *pwi, *LPC_phi, *LPC_psi;

	double *dt, *wploc, T;

	double *pr;
	mwIndex *ir, *jc;
    int ncol_coord;
    
    mxArray *evalrhs[1000], *jacrhs[5], *evalrhs2[1000];
	mxArray *evallhs[1], *jaclhs[1],*evallhs2[1];
	double *xtmp, *xtmp2, *xtmpn;

	int filled, elementcounter, remm;
	int i,j,k,l,l2,m;	/* Indexation variables */

	int *range1, *range2, *range3, *range4;

	double *jac, *bcjac, *sysjac, *zero, *ptmp;
	double *frhs, *frhs2, *frhstmp;

	double *Tcol, *tempmatrix;
	double tmpperiod;

	/* Initializations */
	/* --------------- */
	    
	/* Retrieve parameters. */
	x = mxGetPr(prhs[1]);
	p = mxGetPr(prhs[2]);
    pars2 = mxGetPr(prhs[4]);
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
	T = *(mxGetPr(thisfield));		/* period */
	thisfield = mxGetFieldByNumber(prhs[6],0,34);
	ncol_coord = *(mxGetPr(thisfield));		
    thisfield = mxGetFieldByNumber(prhs[6],0,39);    
	wp = mxGetPr(thisfield);	/* Derivative weights of collocation points */
    /* Kronecker product of the derivative weights and the identity matrix */
	wploc = mxCalloc(mxGetN(thisfield)*mxGetM(thisfield),sizeof(double));
	thisfield = mxGetFieldByNumber(prhs[6],0,41);
	pwi = mxGetPr(thisfield);	/* Extension of weights */
	thisfield = mxGetFieldByNumber(prhs[6],0,60);
	LPC_phi = mxGetPr(thisfield);	/* LPC_phi */
    thisfield = mxGetFieldByNumber(prhs[6],0,61);
	LPC_psi = mxGetPr(thisfield);	/* LPC_psi */
	/* Column numbers of period and free parameters*/
	pars = mxCalloc(nfreep,sizeof(int));
	for (i=0; i<nfreep; i++)
		*(pars+i) = ncoords+i;
    /* Sparse matrix as returnvalue */
	plhs[0] = mxCreateSparse(ncoords+2,ncoords+2,ncoords*ncoords,mxREAL);
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
      	/* Parameters for second rhs-evaluation-call to Matlab */
	evalrhs2[0] = (struct mxArray_tag*) prhs[0];
	evalrhs2[1] = mxCreateDoubleMatrix(1,1,mxREAL);
    zero = mxGetPr(evalrhs2[1]);
    *zero = 0;
	/*evalrhs[2] = mxCreateDoubleMatrix(nphase+ActiveParams,1,mxREAL);*/
    evalrhs2[2] = mxCreateDoubleMatrix(nphase,1,mxREAL);
	xtmp2 = mxGetPr(evalrhs2[2]); 
    for (i=0; i<mxGetNumberOfElements(prhs[2]); i++) {
        evalrhs2[i+3] = mxCreateDoubleMatrix(1,1,mxREAL);
        ptmp = mxGetPr(evalrhs2[i+3]);   
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
	bcjac = mxCalloc(ncoords,sizeof(double));
	frhs = mxCalloc(nphase*ncol,sizeof(double));
    frhs2 = mxCalloc(nphase*ncol,sizeof(double));
	sysjac = mxCalloc(nphase*ncol*(ncol+1)*nphase,sizeof(double));

	Tcol = mxCalloc((tps-1)*nphase+1,sizeof(double));
    
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
        
            /* Call to Matlab for evaluation of rhs */
        for (m=0; m<(ncol+1);m++){
            for (k=0; k<nphase; k++) {     
                *(xtmp2+k) = (*(x+*(range1+m)*nphase+k));
            }
            mexCallMATLAB(1,evallhs2,3+mxGetNumberOfElements(prhs[2]),evalrhs2,"feval");
            frhstmp = mxGetPr(evallhs2[0]);  
            for (k=0; k<nphase; k++){
                *(frhs2+m*nphase+k) = *(frhstmp+k);
            }
          /*  mxDestroyArray(evallhs2[0]);*/
        }
        k=0;
        for (l=0;l<nphase;l++){
            for (m=l;m<(ncol+1)*nphase;m=m+nphase) {             
            /* Compute elements of third component */                
                *(bcjac + *(range3+k)) = *(bcjac + *(range3+k)) + (*(dt+i)) * (*(frhs2+m)) * (*(pwi+m));              
                ++k;  
            }
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
         /*   mxDestroyArray(evallhs[0]);*/

            
			/* Call to Matlab for evaluation of jacobian */
			mexCallMATLAB(1,jaclhs,5,jacrhs,"cjac");
			frhstmp = mxGetPr(jaclhs[0]);

			/* Store jacobian */
			for (k=0; k<nphase*nphase; k++) 	
				*(jac+k) = *(frhstmp+k);
          /*  mxDestroyArray(jaclhs[0]);*/
            
			/* temporary sysjac */
			for (k=0; k<nphase; k++) {
				/* sysjac stores kronecker product of jacobian and weights */
			    for (l=nphase; l<(ncol+2)*nphase; l++) {
			        l2 = floor(l/nphase)-1;
			        remm = l % nphase;
    			    *(sysjac + (l-nphase)*nphase*ncol + (*(range4+k))) = *(wt+j*(ncol+1)+l2) * (*(jac+remm*nphase+k));
    			}
            }
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
				    *(pr + elementcounter) = *(bcjac + *(range3+k));
				    *(ir + elementcounter) = ncoords;
					elementcounter = elementcounter + 1;
				    /* Fill in previously computed element from third component */
				    *(pr + elementcounter) = *(LPC_phi + *(range3+k));                    
				    *(ir + elementcounter) = ncoords+1;
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
        
		/* The last  columns of the jacobian also will be changed at each pass through the loop,
		and therefore will be only effectively stored at the very end. So we store intermediate results temporarily */
            for (j=0; j<ncol*nphase; j++) {
                 *(Tcol + *(range2+j)) = -(*(frhs+j));
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
	}
	/* Finally, store the last 3 columns in the sparse return matrix 
    if (nfreep == 1)*/
    for (j=0; j<(tps-1)*nphase; j++) {
        /* Store the column from the period */
        if (*(Tcol + j)) {
            *(pr + elementcounter) = *(Tcol + j);
            *(ir + elementcounter) = j;
            elementcounter = elementcounter+1;   
        }
    }
    *(pr + elementcounter) = *(LPC_phi+ncoords);        
    *(ir + elementcounter) = ncoords;
    ++elementcounter;
    *(jc + ncoords + 1) = elementcounter;
    for (j=0; j<ncoords+1; j++) {
         /* Store the last column*/
        *(pr + elementcounter) = *(LPC_psi + j);
        *(ir + elementcounter) =  j;
        elementcounter = elementcounter+1;
    }  /* Finish the column */
    *(jc + ncoords + 2) = elementcounter;
	/* Free all allocated memory */    
	/* ------------------------- */

	mxFree(pars);
	mxFree(wploc);
    mxFree(ActiveParams);
	
	mxFree(range1);
	mxFree(range2);
	mxFree(range3);
	mxFree(range4);

	mxFree(jac);
	mxFree(sysjac);
    mxFree(bcjac);
	mxFree(frhs);
    mxFree(frhs2);

	mxFree(Tcol);
	mxFree(tempmatrix); 
    
    /*mxDestroyArray(evalrhs[1]);
	mxDestroyArray(evalrhs[2]);
    for (i=0; i<mxGetNumberOfElements(prhs[2]); i++) {
        mxDestroyArray(evalrhs[3+i]);
    }
    mxDestroyArray(evalrhs2[1]);
	mxDestroyArray(evalrhs2[2]);
    for (i=0; i<mxGetNumberOfElements(prhs[2]); i++) {
        mxDestroyArray(evalrhs2[3+i]);
    }*/

			
	return;
}
