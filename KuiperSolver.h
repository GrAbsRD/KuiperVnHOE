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
#include <stdbool.h>

typedef struct{
	double  c; // critical value,  \alpha = \Pr(\sqrt{n}V_n>c}
	double  v; // quantile for upper tail probability, v = c/sqrt{n}
}KuiperPair;

void CheckParameter(double alpha, unsigned n, unsigned k);


typedef double (*Ptr2Fun)(double, double, unsigned, unsigned k); // for f(c, alpha, n, k) 
typedef double (*Ptr2FunUpdate)(Ptr2Fun f, double guess, double alpha, unsigned n, unsigned k); 
typedef double (*Ptr2FunDist)(double, double);

double FunA0(unsigned n, unsigned k);
double FunA1(double c, unsigned n, unsigned k);
double FunA2(double c, unsigned n, unsigned k);
double FunFnlm(double c, double alpha, unsigned n, unsigned k);
double FunFctm(double c, double alpha, unsigned n, unsigned k);
double FunFste(double c, double alpha, unsigned n, unsigned k);

/* generate the initial value with bisection method */
double GenInitValue(
	Ptr2Fun f, 
	double a, 
	double b, 
	double h, 
	double alpha, 
	unsigned n, 
	unsigned k);


double UpdateMethodDirect(Ptr2Fun f, double c, double alpha, unsigned n, unsigned k);
double UpdateMethodNewton(Ptr2Fun f, double c, double alpha, unsigned n, unsigned k);

double Distance(double x, double y);
KuiperPair KuiperPairSolver(double alpha, unsigned n, int method, unsigned k);

/*********************** Attention, please! *********************
* KuiperUTQ(x, n, k) is equivalent to KuiperLTQ(1-x, n, k)
*****************************************************************/
double KuiperUTQ(double upper_tail_prob, unsigned n, unsigned k);  
double KuiperLTQ(double lower_tail_prob, unsigned n, unsigned k); 

/*********************** Attention, please! *********************
* KuiperInvCDF(x, n) is equivalent to KuiperLTQ(x, n) 
*****************************************************************/
double KuiperInvCDF(double prob, unsigned n, unsigned k); // prob is in [0,1)


/*********************** Attention, please! *********************
* ModiKuiperUTQ(x) is equivalent to ModiKuiperLTQ(1-x)
*****************************************************************/
double ModiKuiperUTQ(double alpha);
double ModiKuiperLTQ(double alpha);








