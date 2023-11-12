/*
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * Author:  Hong-Yan ZHANG
 * Contact: hongyan@hainnu.edu.cn
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "KuiperSolver.h"

void CheckParameter(double alpha, unsigned n, unsigned k)
{
	unsigned UtpIsGood = 1;
	if(alpha < 0.0 || alpha > 1.0){
		printf("Input error, the alpha is not the upper tail probability\n");
        UtpIsGood = 0;
	}
	unsigned CapacityIsGood = 1;
	if(n < 6){
		printf("Input error, the sample capacity n should be n >= 6");
		CapacityIsGood = 0;
	}

	unsigned OrderIsGood = 1;
	if(k == 0){
		printf("Input error, the order k should be k >= 1");
		OrderIsGood = 0;
	}
	
	if(UtpIsGood*CapacityIsGood*OrderIsGood == 0) exit(0);	
}



/***************** Solving the Fixed-Point *******************/
double FixedPointSolver(
	Ptr2FunUpdate T, 
	Ptr2Fun f, 
	Ptr2FunDist d, // the implementation is the function Distance 
    double precision, 
    double guess, 
    double alpha, 
    unsigned n,
    unsigned k)
{
	double improve = T(f, guess, alpha, n, k);
	while( d(guess, improve) >= precision)
	{
		guess = improve;
		improve = T(f, guess, alpha, n, k);		
	}
	return improve;	
}

/***************** Calculating the distance ******************/ 
double Distance(double x, double y)
{
	double dist = fabs(x - y);
	return dist;	
}


/*********** Auxiliary procedures **********************/
double FunA0(unsigned n, unsigned k)
{
	double A0 = -1.0;             // k = 1
	if(k > 1)  A0 += 1.0/(18.0*n);  // k = 2, 3      
    if(k > 3)  A0 += - 1.0/(648.0*n*n); // k = 4, 5, ...
    return A0;	
}


double FunA1(double c, unsigned n, unsigned k)
{
	double  A1 = (8*c*c -2) - 8*(4*pow(c,3)-3*c)/(3*sqrt(n)); // k = 1
    if( k > 1){  // k = 2, 3, ...
		A1 += (64*pow(c,4) - 100*pow(c,2) + 13)/(9*n);  
	}
	if( k > 2){ // k = 3, 4, ...
		A1 += -32*(8*pow(c,5)-22*pow(c,3)+9*c)/(81*pow(n,1.5));
	}
	if( k > 3){ // k = 4, 5, ...
		A1 += (1024*pow(c,6)-4496*pow(c,4)+3864*c*c-363)/(972*n*n);
	}
	if( k > 4){  // k = 5, 6, ...
		A1 += -32*(512*pow(c,7)-3376*pow(c,5)+5080*pow(c,3)
		      -1485*c)/(3645*pow(n,2.5));
	}
	return A1;
}

double FunA2(double c, unsigned n, unsigned k)
{
	double  A2  = (32*c*c -2) - 32*(16*pow(c,3) - 3*c)/(3*sqrt(n)); // k = 1
	if( k > 1){ // k = 2, 3, ...
		A2 += (4096*pow(c,4)- 1552*c*c + 49)/(9*n);
	}
	if( k > 2){ // k = 3, 4, ...
		A2 += -64*(1024*pow(c,5)-656*pow(c,3)+63*c)/(81*pow(n,1.5));
	}
	if( k > 3){ // k = 4, 5, ...
		A2 += (1048576*pow(c,6) - 1024256*pow(c,4)
		       + 199776*c*c - 2403)/(972*n*n);
	}
	if( k > 4){ // k = 5, 6, ...
		A2 += - 32*(2097152*pow(c,7) - 2919424*pow(c,5)
		      + 964480*pow(c,3) - 63540*c)/(3645*pow(n, 2.5));
	} 
	return A2;	
}

/*******  for V_n test, direct iterative method  ********/
double FunFctm(double c, double alpha, unsigned n, unsigned k)
{
	double A0 = FunA0(n, k);
	double A1 = FunA1(c, n, k);
	double A2 = FunA2(c, n, k);
	double y = sqrt((log(A1 + A2*exp(-6*c*c)) - log(alpha-1-A0))/2);
	return y;	
}

/*******  for V_n test, Newton iterative method  ********/
double FunFnlm(double c, double alpha, unsigned n, unsigned k)
{
	double A0 = FunA0(n, k);
    double A1 = FunA1(c, n, k);
    double A2 = FunA2(c, n, k);
    double y = 2*c*c + log(alpha - 1 - A0) - log(A1 + A2*exp(-6*c*c));
    return y;
}


double FunFste(double c, double alpha, unsigned n, unsigned k)
{
	/* Note that it is for Stephens' modified Kuiper Test 
	 * the value returned  does not depend on the n and k.
	 * The n and k is used to keep the interface of the function
	 * so that FunFste can be referenced by the pointer f which
	 * is an argument (function object) of the fixed-point solver.
	 */
    double y = (-2+8*c*c)*exp(-2*c*c) - alpha;
	//double y = alpha*exp(2*c*c) -(8*c*c-2);
	return y;	
}


double GenInitValue(Ptr2Fun f, double a, double b, double h, 
                    double alpha, unsigned n, unsigned k)
{
        double guess = (b+a)/2; /* set the midian point as the guess */ 
        double Delta = fabs(b-a);
        while( Delta > h){ // h = 0.05 is good enough
                if(f(guess, alpha, n, k)*f(a, alpha, n, k) > 0){
        	   	a = guess;
        	}else{
        		b = guess;          	
        	}
        	guess = (a+b)/2;
        	Delta = fabs(b-a)/2;
        }
    return guess;	
}

/************* Updating Mapping *********************/
double UpdateMethodDirect(Ptr2Fun f_ctm, double c, double alpha, unsigned n, unsigned k)
{
	double  u = f_ctm(c, alpha, n, k);  // f_ctm is a contractive map;
	return  u;	
}
	
double UpdateMethodNewton(Ptr2Fun f_nlm, double c, double alpha, unsigned n, unsigned k)
{
  double h = 1e-5;
  double slope = (f_nlm(c+h, alpha, n, k) - f_nlm(c, alpha, n, k))/h;
  double c_new = c - f_nlm(c, alpha, n, k)/slope;
  return c_new;
}


/******************* Sovling the Kuiper's Pair ********************/
KuiperPair KuiperPairSolver(double alpha, unsigned n, int method, unsigned k)
{
	Ptr2FunUpdate T;
	Ptr2Fun f;
	
	
	if(method == 1) {// Direct iterative for Kuiper's $V_n$-test /
		T = UpdateMethodDirect;
		f = FunFctm; // $f_{ctm}$ 
	}else{           // Newton's iterative for Kuiper's $V_n$-test 
		T = UpdateMethodNewton;
		f = FunFnlm; // $f_{nlm}$
	}
	
	double  epsilon = 1e-5;
	Ptr2FunDist d = Distance;	
	double a = 0.6, b = 3.0, h = 0.05;
	double guess = GenInitValue(FunFnlm,a, b, h, alpha, n, k);
    double c = FixedPointSolver(T, f, d, epsilon, guess, alpha, n, k);
    double v = c/sqrt(n);	
	KuiperPair sol = {c, v};	
	return sol;
}

/*********************** Attention, please! *********************
* KuiperUTQ(x, n, k) is equivalent to KuiperLTQ(1-x, n, k)
*****************************************************************/
double KuiperLTQ(double lower_tail_prob, unsigned n, unsigned k)
{
	if(lower_tail_prob <= 1e-4) {return 0.0;}
	
	const int method = 2;
	double alpha = 1.0 - lower_tail_prob;// Upper tail probability	
	double a = 0.6, b = 3.0, h = 0.05;
	double guess = GenInitValue(FunFnlm, a, b, h, alpha, n, k);
	KuiperPair sol = KuiperPairSolver(alpha,n,method,k);
    return sol.v;	
}

/*********************** Attention, please! *********************
* KuiperUTQ(x, n, k) is equivalent to KuiperLTQ(1-x, n, k)
*****************************************************************/
double KuiperUTQ(double upper_tail_prob, unsigned n, unsigned k)
{
	if(upper_tail_prob >= 1.0 - 1e-4) {return 0.0;}
	
	Ptr2FunUpdate  T = UpdateMethodNewton;
	Ptr2Fun        f = FunFnlm;
	double       eps = 1e-5;
	Ptr2FunDist    d = Distance;	
	double a = 0.6, b = 3.0, h = 0.05;
	double guess = GenInitValue(FunFnlm,a, b, h, upper_tail_prob, n, k);	
	double c = FixedPointSolver(T,f,d,eps,guess,upper_tail_prob,n, k);
	double utq = c/sqrt(n); // upper tail quantile
	return utq;	
}

double KuiperInvCDF(double prob, unsigned n, unsigned k)
{
	// double x = KuiperLTQ(prob, n, k) is OK
	double x = KuiperUTQ(1.0 - prob, n, k); 	
    return x;	
}


double ModiKuiperUTQ(double alpha)
{
	if(alpha >= 1.0 - 1e-4) {return 0.0;}
	
	Ptr2FunUpdate  T = UpdateMethodNewton;
	Ptr2Fun        f = FunFste;
	double       eps = 1e-5;
	Ptr2FunDist    d = Distance;	
	//double a = 0.8, b = 3.0, h = 0.05;
	//double guess = GenInitValue(f,a, b, h, alpha, 6, 1);	
	double guess = 1.8; 
	double c = FixedPointSolver(T, f, d, eps, guess, alpha, 6, 1);
	return c;	
}


/* This is an equivalent implementation of ModiKuiperUTQ above.
 * Note that here INTEGER = 100000000 is a replacement of infty
 * for the limit and the overflow of n*n may happen in computing
 * the value of A0, A1 and A2 with C/C++/... programs.
 */ 

/*
double ModiKuiperUTQ(double alpha)
{
	if(alpha >= 1.0 - 1e-4) {return 0.0;}
	
	Ptr2FunUpdate  T = UpdateMethodNewton;
	Ptr2Fun        f = FunFnlm;
	double       eps = 1e-5;
	Ptr2FunDist    d = Distance;	
	double a = 0.6, b = 3.0, h = 0.1;
	unsigned INTEGER = 10000000; // for replacing the infty
	double guess = GenInitValue(FunFnlm,a, b, h, alpha, INTEGER, 1);	
	double c = FixedPointSolver(T,f,d,eps,guess, alpha, INTEGER, 1);
	return c;		
}
*/



/*********************** Attention, please! *********************
* ModiKuiperUTQ(x) = ModiKuiperLTQ(1-x)
*****************************************************************/
double ModiKuiperLTQ(double ltp)
{
    double v = ModiKuiperUTQ(1-ltp);
	return v;	
}





