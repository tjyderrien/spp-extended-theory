/* This is the first of the two examples in the doc of gsl interp.  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

 void wrapper_int_(double y[], double x[], double xint[], double yint[], int *LENGTH, double *STEPSIZE) 

{
  int i;

  double xi, yi;  /*, x[10], y[10];*/
  /* the following only test now
  for (i = 0; i < *LENGTH; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }
  end test */


  {
    /* hansi : use counter for x values */
    int counter;
    counter = 0;
    /* end hansi */

    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_cspline, *LENGTH);

    
    gsl_spline_init (spline, x, y, *LENGTH);
    
    for (xi = x[0]; xi < x[*LENGTH-1]; xi += *STEPSIZE)
      {
	yi = gsl_spline_eval (spline, xi, acc);
	/* instead of printing it out, we now write the new 
	   values to the new arraw xint,yint */
	/*	printf ("%g %g\n", xi, yi);  */
	xint[counter]=xi; 
	yint[counter]=yi; 
	/* 	printf (" %g %g\n",  xint[counter],yint[counter]); */
	/* 	printf ("counter %d \n", counter); 	*/
	counter++;      /*  = counter+1; */
	/* end hansi */
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

}

/*******************************************************************
 ******************************************************************/
/* This is the first of the two examples in the doc of gsl interp.  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

 void wrapper_int_fq_lin_(double y[], double x[], double xint[], double yint[], int *LENGTH, double *STEPSIZE) 

{
  int i;

  double xi, yi;  /*, x[10], y[10];*/
  /* the following only test now
  for (i = 0; i < *LENGTH; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }
  end test */


  {
    /* hansi : use counter for x values */
    int counter;
    counter = 0;
    /* end hansi */

    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_linear, *LENGTH);

    //       printf("linear\n");

    /*	orig: 
	gsl_spline *spline 
	= gsl_spline_alloc (gsl_interp_cspline, *LENGTH);
    */
    
    gsl_spline_init (spline, x, y, *LENGTH);
    
    for (xi = x[0]; xi < x[*LENGTH-1]; xi += *STEPSIZE)
      {
	yi = gsl_spline_eval (spline, xi, acc);
	/* instead of printing it out, we now write the new 
	   values to the new arraw xint,yint */
	/*	printf ("%g %g\n", xi, yi);  */
	xint[counter]=xi; 
	yint[counter]=yi; 
	/* 	printf (" %g %g\n",  xint[counter],yint[counter]); */
	/* 	printf ("counter %d \n", counter); 	*/
	counter++;      /*  = counter+1; */
	/* end hansi */
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

}


/*******************************************************************
 ******************************************************************/
/* This is the first of the two examples in the doc of gsl interp.  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

 void wrapper_int_fq_pol_(double y[], double x[], double xint[], double yint[], int *LENGTH, double *STEPSIZE) 

{
  int i;

  double xi, yi;  /*, x[10], y[10];*/
  /* the following only test now
  for (i = 0; i < *LENGTH; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }
  end test */


  {
    /* hansi : use counter for x values */
    int counter;
    counter = 0;
    /* end hansi */

    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_polynomial, *LENGTH);

    //  printf("polynomial\n");

    /*	orig: 
	gsl_spline *spline 
	= gsl_spline_alloc (gsl_interp_cspline, *LENGTH);
    */
    
    gsl_spline_init (spline, x, y, *LENGTH);
    
    for (xi = x[0]; xi < x[*LENGTH-1]; xi += *STEPSIZE)
      {
	yi = gsl_spline_eval (spline, xi, acc);
	/* instead of printing it out, we now write the new 
	   values to the new arraw xint,yint */
	/*	printf ("%g %g\n", xi, yi);  */
	xint[counter]=xi; 
	yint[counter]=yi; 
	/* 	printf (" %g %g\n",  xint[counter],yint[counter]); */
	/* 	printf ("counter %d \n", counter); 	*/
	counter++;      /*  = counter+1; */
	/* end hansi */
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

}




/*******************************************************************
 ******************************************************************/
/* This is the first of the two examples in the doc of gsl interp.  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

 void wrapper_int_fq_cspline_(double y[], double x[], double xint[], double yint[], int *LENGTH, double *STEPSIZE) 

{
  int i;

  double xi, yi;  /*, x[10], y[10];*/
  /* the following only test now
  for (i = 0; i < *LENGTH; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }
  end test */


  {
    /* hansi : use counter for x values */
    int counter;
    counter = 0;
    /* end hansi */

    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_cspline, *LENGTH);
    //    printf("cspline\n");

    /*	orig: 
	gsl_spline *spline 
	= gsl_spline_alloc (gsl_interp_cspline, *LENGTH);
    */
    
    gsl_spline_init (spline, x, y, *LENGTH);
    
    for (xi = x[0]; xi < x[*LENGTH-1]; xi += *STEPSIZE)
      {
	yi = gsl_spline_eval (spline, xi, acc);
	/* instead of printing it out, we now write the new 
	   values to the new arraw xint,yint */
	/*	printf ("%g %g\n", xi, yi);  */
	xint[counter]=xi; 
	yint[counter]=yi; 
	/* 	printf (" %g %g\n",  xint[counter],yint[counter]); */
	/* 	printf ("counter %d \n", counter); 	*/
	counter++;      /*  = counter+1; */
	/* end hansi */
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

}




/*******************************************************************
 ******************************************************************/
/* This is the first of the two examples in the doc of gsl interp.  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

 void wrapper_int_fq_akima_(double y[], double x[], double xint[], double yint[], int *LENGTH, double *STEPSIZE) 

{
  int i;

  double xi, yi;  /*, x[10], y[10];*/
  /* the following only test now
  for (i = 0; i < *LENGTH; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }
  end test */


  {
    /* hansi : use counter for x values */
    int counter;
    counter = 0;
    /* end hansi */

    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_akima, *LENGTH);

    //    printf("akima spline\n");

    /*	orig: 
	gsl_spline *spline 
	= gsl_spline_alloc (gsl_interp_cspline, *LENGTH);
    */
    
    gsl_spline_init (spline, x, y, *LENGTH);
    
    for (xi = x[0]; xi < x[*LENGTH-1]; xi += *STEPSIZE)
      {
	yi = gsl_spline_eval (spline, xi, acc);
	/* instead of printing it out, we now write the new 
	   values to the new arraw xint,yint */
	/*	printf ("%g %g\n", xi, yi);  */
	xint[counter]=xi; 
	yint[counter]=yi; 
	/* 	printf (" %g %g\n",  xint[counter],yint[counter]); */
	/* 	printf ("counter %d \n", counter); 	*/
	counter++;      /*  = counter+1; */
	/* end hansi */
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

}

















































/*******************************************************************
 ******************************************************************/

/* The only difference in this wrapper is that we replaced 
   gsl_spline_eval by gsl_spline_eval_deriv */

void wrapper_differentiate_(double y[], double x[], double xint[], double yint[], int *LENGTH, double *STEPSIZE) 

{
  int i;

  double xi, yi;  /*, x[10], y[10];*/
  /* the following only test now
  for (i = 0; i < *LENGTH; i++)
    {
      printf ("%g %g\n", x[i], y[i]);
    }
  end test */


  {
    /* hansi : use counter for x values */
    int counter;
    counter = 0;
    /* end hansi */

    gsl_interp_accel *acc 
      = gsl_interp_accel_alloc ();
    gsl_spline *spline 
      = gsl_spline_alloc (gsl_interp_cspline, *LENGTH);
    
    gsl_spline_init (spline, x, y, *LENGTH);
    
    for (xi = x[0]; xi < x[*LENGTH-1]; xi += *STEPSIZE)
      {
	yi = gsl_spline_eval_deriv (spline, xi, acc);
	/* instead of printing it out, we now write the new 
	   values to the new arraw xint,yint */
	/*	printf ("%g %g\n", xi, yi);  */
	xint[counter]=xi; 
	yint[counter]=yi; 
	/*	printf (" %g %g\n",  xint[counter],yint[counter]); */
	/*	printf ("counter %d \n", counter); 	*/
	 counter++;
	/* end hansi */
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }

}


