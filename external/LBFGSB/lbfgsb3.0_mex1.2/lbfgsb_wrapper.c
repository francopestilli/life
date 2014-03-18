/*=================================================================
 *
 * lbfgsb_wrapper.c, mex gatewate interface
 * You can call this directly from Matlab,
 * but I strongly suggest using lbfgsb.m, since many
 * variables are modified in-place, and lbfgsb.m will
 * handle them correctly. It also shows you how to use
 * this routine
 *  (so just modify lbfgsb.m to suit your own taste)
 *
 * Written by Stephen Becker, Feb 14 2012
 *
 * This assumes you have installed lbfgsb 3.0
 * (free; see http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html,
 *  direct download is:
 * http://users.eecs.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz )
 *
 * There are also versions 2.1 and 2.4 of the library
 *      For v 2.1, Peter Carbonetto's mex interface works; see
 * http://www.mathworks.com/matlabcentral/fileexchange/15061-matlab-interface-for-l-bfgs-b
 *      and also
 *  http://www.cs.ubc.ca/~pcarbo/lbfgsb-for-matlab.html
 *
 *      For v 2.4 (I don't know where you can find this version though),
 *      use the mex files from here:
 * http://www.cs.toronto.edu/~liam/software.shtml
 *
 *  
 *  You should use version 3.0! It's better.  The old versions are basically
 *  15 years old; changes in version 3.0 are fundamental (see the 2011 "remark" paper
 *  listated at Nocedal's website)
 *
 *
 *
 *
 *  Using a string routine from Peter Carbonetto
 *      It is copyright (c) 2009 Peter Carbonetto
 *      All rights reserved.
 *      Neither Peter Carbonetto nor the University of
 *      British Columbia nor its contributors endorse
 *      this derived work.
 *
 * Inputs (in order):
 *  m       # of memory vectors to keep
 *  x       initial guess or current point, also used to determine the size of the problem (Nx1)
 *  l       list of upper bounds (Nx1)
 *  u       list of lower bounds (Nx1)
 *  nbd     list of which bounds are active (Nx1):
 *              0=neither u nor l, 1 = only l, 2 = both u and l, 3 = only u
 *  f       value of function (1x1)
 *  g       gradient of function (Nx1)
 *
 *  factr   stopping crit: 1e+12 for low accuracy, 1e7 for moderate, 1e1 for high accuracy
 *              (will be multiplied by machine precision)
 *  pgtol   stopping crit for infinity norm of projected gradient
 *  wa      work space array (double)
 *  iwa     work space array (int)
 *  task    string
 *  iprint  how verbose the fortran program should be
 *  csave       some iformation (string, length 60)
 *  lsave[4]    some information (logicals)
 *  isave[44]   some information (integers)
 *  dsave[29]   some information (doubles)
 *
 * Outputs:
 *  f, task, csave, lsave, isave, dsave
 *
 * Warning: the following variables are modified in-place
 *  x, g, wa, iwa
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

/*#include <stdbool.h>*/ /* for using the 'bool' data type */
/* actually, matrix.h defines mxLogical, which is normally 'bool' */

/* matrix.h includes tmwtypes.h, which defines mwSize=size_t
 * unless MX_COMPAT_32 has been defined, in which case mwSize is
 * defined as int.  So assuming mwSize=size_t,
 * then if you pass in an 'int', it will automatically be scaled
 * up if necessary, so there's no problem.
 * You probably do NOT want size variables, like 'm' and 'n',
 * to be size_t in size, because the fortran program will have
 * problems then.
 * */


#ifdef _BLAS64_
#define int_F ptrdiff_t
#else
#define int_F int
#endif

#include <string.h>
#include <limits.h> /* for CHAR_BIT */
#include <assert.h>

/* These are the fortran arguments in order */
#define N_m 0
#define N_x 1
#define N_l 2
#define N_u 3
#define N_nbd 4
#define N_f 5
#define N_g 6
#define N_factr 7
#define N_pgtol 8
#define N_wa 9
#define N_iwa 10
#define N_task 11
#define N_iprint 12
#define N_csave 13
#define N_lsave 14
#define N_isave 15
#define N_dsave 16

/* these depend on the lbfgsb release version */
#define LENGTH_STRING 60
#define LENGTH_LSAVE 4
#define LENGTH_ISAVE 44
#define LENGTH_DSAVE 29

/* the fortran lbfgsb code uses lsave as type LOGICAL
 * which is usually 'int' in C, whereas type LOGICAL*1
 * would be 'bool' in C.
 * BTW, the 'mxLogical' is usually defined as 'bool' */
typedef int fortranLogical;

#ifdef DEBUG
#define debugPrintf mexPrintf
#else
#define debugPrintf fakePrintf
#endif
/* If the 'DEBUG' symbol is undefined, then don't print: */
int fakePrintf(const char *format, ...){
    return 0;
}

/* For using fortran programs in C */
/* On Windows, we generally do not need to append the underscore to the end,
 * whereas on Linux we do.
 * Also, if using a native Windows fortran compiler, we may need
 * to put the function names in all capital letters
 * */
#if defined(NOUNDERSCORE)
    #if defined(UPPERCASE_FORTRAN)
        #define setulb SETULB
    #endif
#else
    #if defined(UPPERCASE_FORTRAN)
        #define setulb SETULB_
    #else
        #define setulb setulb_
    #endif
#endif

/* Declare the L-BFGS-B function */
void setulb( int_F *n, int_F *m, double *x, double *l, double *u,
        int_F *nbd, double *f, double *g, double *factr,
        double *pgtol, double *wa, int_F *iwa,
        char *task, int_F *iprint, char *csave, fortranLogical *lsave,
        int_F *isave, double *dsave );

/* This is taken from the other l-bfgs-b mex interface by Peter Carbonetto. */
/* Copy a C-style string (a null-terminated character array) to a
 * non-C-style string (a simple character array). The length of the
 * destination character array is given by "ndest". If the source is
 * shorter than the destination, the destination is padded with blank
 * characters. 
 * */
void copyCStrToCharArray (const char* source, char* dest, int ndest) {
  int i;
  int nsource = strlen(source);
  /* Only perform the copy if the source can fit into the destination. */
  if (nsource < ndest) {
    strcpy(dest,source);

    /* Fill in the rest of the string with blanks. */
    for (i = nsource; i < ndest; i++)
      dest[i] = ' ';
  }
}

mxLogical isInt( const mxArray *pm ) {
    /* What size 'int' does the fortran program
     * expect ? Not sure... But let's hope that
     * it's the same size as the "int" in C.
     * On my 64-bit computer, CHAR_BIT = 8
     * and sizeof(int) = 4, so it's still 32 bits 
     *
     * CHAR_BIT is from limits.h
     * If using gcc, you can run `gcc -dM -E - < /dev/null | grep CHAR_BIT`
     *  and it should define the symbol __CHAR_BIT__, so this is another way.
     * */
    
    /* debugPrintf("Sizeof(int) is %d\n", sizeof(int) ); */
    switch ( CHAR_BIT * sizeof(int) ) {
        case 16 :
            return mxIsInt16(pm);
        case 32:
            return mxIsInt32(pm);
        case 64:
            return mxIsInt64(pm);
        default:
            mexErrMsgTxt("You have a weird computer that I don't know how to support");
            return false;
    }
}


/* Main mex gateway routine */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] )   { 
    
    int_F   iprint = 1;
    char    task[LENGTH_STRING], csave[LENGTH_STRING];
    char    *task_m, *csave_m;
    double  *isave_m, *lsave_m, *dsave_m; /* temporary pointers */
    int_F   *isave_int;
    mxLogical  *lsave_bool;
    fortranLogical    lsave[LENGTH_LSAVE];
    int_F   n, m, *nbd, *iwa, isave[LENGTH_ISAVE];
    double  *nbd_dbl, *iwa_dbl;
    double  f, factr, pgtol, *x, *l, *u, *g, dsave[LENGTH_DSAVE], *wa;
    int     i;
    mxClassID classID;
    
    int ndim = 2; /* for lcc compiler, must declare these here, not later ... */
    mwSize dims[2] = { LENGTH_ISAVE, 1 };
    /* for copying the data */
    mxLogical   *lsave_out;
    int         *isave_out;
    double      *dsave_out;
    
    /* Parse inputs. Quite boring */
    
    if (nrhs < 5 ) mexErrMsgTxt("Needs at least 5 input arguments");
    m       = (int)*mxGetPr( prhs[N_m] );
    n       = mxGetM( prhs[N_x] );
    if ( mxGetN(prhs[N_x]) != 1 ) mexErrMsgTxt("x must be a column vector");
    if ( mxGetM(prhs[N_l]) != n ) mexErrMsgTxt("l must have same size as x");
    if ( mxGetM(prhs[N_u]) != n ) mexErrMsgTxt("u must have same size as x");
    if ( mxGetM(prhs[N_nbd]) != n ) mexErrMsgTxt("nbd must have same size as x");
    x       = mxGetPr( prhs[N_x] );
    l       = mxGetPr( prhs[N_l] );
    u       = mxGetPr( prhs[N_u] );
    if ( isInt( prhs[N_nbd] ) ) {
        nbd     = (int_F *)mxGetData( prhs[N_nbd] ); 
    } else {
        nbd_dbl = mxGetPr( prhs[N_nbd] );
        nbd     = (int_F *)mxMalloc( n * sizeof(int) );
        assert( nbd != NULL );
        /* convert nbd_dbl (in double format) to integers */
        for (i=0;i<n;i++)
            nbd[i]  = (int)nbd_dbl[i];
    }

    /* f and g, the function and its gradient */
    if ( nrhs < N_f+1 ) 
        f   = 0.0;
    else if  (mxGetNumberOfElements( prhs[N_f] )!=1)
        f   = 0.0;
    else{
        f   = *mxGetPr( prhs[N_f] );
    }
    if (nrhs < N_g-1 ) {
        g   = (double *)mxMalloc( n * sizeof(double) );
        assert( g != NULL );
    }else{
        if ( mxGetM(prhs[N_g]) != n ) mexErrMsgTxt("g must have same size as x");
        g   = mxGetPr( prhs[N_g] );
    }
    
    /* some scalar parameters */
    if ( nrhs < N_factr+1 ) 
        factr   = 1.0e7;
    else if (mxGetNumberOfElements( prhs[N_factr] )!=1)
        factr   = 1.0e7;
    else
        factr   = *mxGetPr( prhs[N_factr] );
    
    if ( nrhs < N_pgtol+1 ) 
        pgtol   = 1.0e-5;
    else if (mxGetNumberOfElements( prhs[N_pgtol] )!=1)
        pgtol   = 1.0e-5;
    else
        pgtol   = *mxGetPr( prhs[N_pgtol] );
    
    
    /* the work arrays 'wa' and 'iwa' */
    if ( nrhs < N_wa+1 ){
        debugPrintf("Allocating memory for wa variable\n"); 
        wa      = (double *)mxMalloc( (2*m*n + 5*n + 11*m*m + 8*m ) * sizeof(double) );
        assert( wa != NULL );
    } else {
        if ( mxGetM(prhs[N_wa]) < (2*m*n + 5*n + 11*m*m + 8*m ) )
            mexErrMsgTxt("work array is too small; make it 2*m*n + 5*n + 11*m*m + 8*m x 1");
        wa      = mxGetPr( prhs[N_wa] );
    }
    
    if ( nrhs < N_iwa+1 ){
        debugPrintf("Allocating memory for iwa variable\n");
        iwa     = (int_F *)mxMalloc( (3*n)*sizeof(int) );
        assert( iwa != NULL );
    }else {
        if ( mxGetM(prhs[N_iwa]) < (3*n) )
            mexErrMsgTxt("i_work array is too small; make it 3*n x 1");
        if (isInt( prhs[N_iwa] )) {
            iwa      = (int_F *)mxGetData( prhs[N_iwa] );
        }else{
            debugPrintf("Converting iwa array to integers\n" );
            iwa_dbl  = mxGetPr( prhs[N_iwa] );
            iwa      = (int_F *)mxMalloc( (3*n)*sizeof(int) );
            for (i=0;i<(3*n);i++){
                iwa[i]  = (int)iwa_dbl[i];
            }
        }
    }
    
    /* the 'task' string */
    if ( nrhs < N_task+1 )
        copyCStrToCharArray("START", task, LENGTH_STRING );
    else{
        task_m = mxArrayToString( prhs[N_task] );
        assert(task_m != NULL);
        copyCStrToCharArray(task_m, task, LENGTH_STRING );
        /* debugPrintf("String2 is: %s, and strlen is %d\n", task_m, strlen(task_m) ); */
    }
    /* debugPrintf("String is: %s, and strlen is %d\n", task, strlen(task) ); */
    

    if ( nrhs < N_iprint+1 )
        iprint  = 1;
    else if (mxGetNumberOfElements( prhs[N_iprint] )!=1)
        iprint  = 1;
    else
        iprint = (int)*(int*)mxGetData( prhs[N_iprint] );
    
    
    
    
    
    /* === Deal with the csave, lsave, isave and dsave variables === */

    if ( nrhs < N_csave+1 ){
        /* do nothing */
    }else {
        csave_m = mxArrayToString( prhs[N_csave] );
        assert(csave_m != NULL);
        copyCStrToCharArray(csave_m, csave, LENGTH_STRING );
    }
    
    if ( nrhs < N_lsave+1 ) {
        /* do nothing, leave at default value */
    }else if (mxGetNumberOfElements( prhs[N_lsave] )!=LENGTH_LSAVE){
        /* do nothing, leave at default value */
        debugPrintf("warning: ignoring the lsave variable, since it has wrong length\n");
    }else{
        /* copy them */
        if (mxIsLogical( prhs[N_lsave] ) ) {
            lsave_bool  = (mxLogical*)mxGetData(prhs[N_lsave]);
            for(i=0;i<LENGTH_LSAVE;i++)
                lsave[i]    = (fortranLogical)lsave_bool[i];
        } else {
            lsave_m     = (double*)mxGetData(prhs[N_lsave]);
            for(i=0;i<LENGTH_LSAVE;i++)
                lsave[i]    = (fortranLogical)lsave_m[i];
        }
    }
           
    
    if ( nrhs < N_isave+1 ) {
        /* do nothing, leave at default value */
    }else if (mxGetNumberOfElements( prhs[N_isave] )!=LENGTH_ISAVE){
        /* do nothing, leave at default value */
        debugPrintf("warning: ignoring the isave variable, since it has wrong length\n");
    }else{
        /* copy them */
        if (mxIsInt32( prhs[N_isave] )) {
            /* debugPrintf("isave input is type int32\n"); */
            isave_int     = (int_F *)mxGetData(prhs[N_isave]);
            for(i=0;i<LENGTH_ISAVE;i++)
                isave[i]    = isave_int[i];
        } else {
            isave_m     = (double*)mxGetData(prhs[N_isave]);
            for(i=0;i<LENGTH_ISAVE;i++)
                isave[i]    = (int)isave_m[i];
        }

    }
    
    if ( nrhs < N_dsave+1 ) {
        /* do nothing, leave at default value */
    }else if (mxGetNumberOfElements( prhs[N_dsave] )!=LENGTH_DSAVE){
        /* do nothing, leave at default value */
    }else{
        /* we don't have to re-cast, but we should copy
         * the data anyhow, since 'dsave' is a stack variable */
        dsave_m     = (double*)mxGetData(prhs[N_dsave]);
        for(i=0;i<LENGTH_DSAVE;i++)
            dsave[i]    = dsave_m[i];
    }
          
    
            
    /* -- Finally, done with parsing inputs. Now, call lbfgsb fortran routine */
    
    /* Be careful! This modifies many variables in-place! 
     * Basically, anything without a '&' before it will be changed in the Matlab
     * workspace */
    
    setulb(&n,&m,x,l,u,nbd,&f,g,&factr,&pgtol,wa,iwa,task,&iprint,
            csave,lsave,isave,dsave);
            
            
    /* don't need to free the memory if I return it in a variable,
     * since then Matlab will handle it automatically */
            
    
    
    /* -- Now, another tedious task: copy variables to output --*/

    if (nlhs != 6 )    mexErrMsgTxt("Should have 6 output arguments");
    

    /* Outputs:
     * f, task, and the "save" variables: csave, lsave, isave, dsave
     * */
    
    
    
    plhs[0] = mxCreateDoubleScalar( f );
    plhs[1] = mxCreateString( task );
    plhs[2] = mxCreateString( csave );  
    plhs[3] = mxCreateLogicalMatrix( LENGTH_LSAVE, 1 );
    
    /* Saving integers... */
    switch ( CHAR_BIT * sizeof(int) ) {
        case 16 :
            classID     = mxINT16_CLASS;
            break;
        case 32:
            classID     = mxINT32_CLASS;
            break;
        case 64:
            classID     = mxINT64_CLASS;
            break;
        default:
            debugPrintf("Value of CHAR_BIT: %d, sizeof(int): %d, together: %d\n",
                    CHAR_BIT, sizeof(int), CHAR_BIT*sizeof(int) );
            mexErrMsgTxt("You have a weird computer that I don't know how to support");
    }
    plhs[4] = mxCreateNumericArray( ndim, dims, classID, mxREAL );
    
    /* and 'double' types are so much easier */
    plhs[5] = mxCreateDoubleMatrix(LENGTH_DSAVE,1, mxREAL);
    
    assert( plhs[0] != NULL );
    assert( plhs[1] != NULL );
    assert( plhs[2] != NULL );
    assert( plhs[3] != NULL );
    assert( plhs[4] != NULL );
    assert( plhs[5] != NULL );
    
    /* copy the data */
    
    lsave_out     = (mxLogical*)mxGetData(plhs[3]);
    for (i=0;i<LENGTH_LSAVE;i++)
        lsave_out[i] = (mxLogical)lsave[i];
    
    isave_out     = (int*)mxGetData(plhs[4]);
    for (i=0;i<LENGTH_ISAVE;i++)
        isave_out[i] = (int)isave[i];
    
    dsave_out     = (double*)mxGetData(plhs[5]);
    for (i=0;i<LENGTH_DSAVE;i++)
        dsave_out[i] = (double)dsave[i];
    
    return;
}
