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

#include <stdio.h>
#include <stdlib.h>

#include "KuiperSolver.h"

/*********** Compile+Link+Run at Unix/Linux terminal *****************
 * Compile+Link:   $ gcc KuiperPairSolver.c main.c -o Kuiper -lm
 *           or    $ g++ KuiperPairSolver.c main.c -o Kuiper
 * Run:            $ ./Kuiper 0.05 1 
 *                 $ ./Kuiper 0.05 2 
 *                 $ ./Kuiper 0.05 1  >> test.txt
 *                 $ ./Kuiper 0.05 2  >> test.txt
 * 
 * *******************************************************************/
 
void PrintArray1d(unsigned a[], unsigned size, const char* arrayname)
{
   printf("%s = [", arrayname);
   for(int i = 0; i < size; i++){
	   printf("%6d", a[i]);
	   if(i < size -1)
			printf(",");
	   else
			printf("];\n");
   }
}

void PrintUtqSeq(double UTQ[], unsigned size, const char* arrayname)
{
   printf("%s = [", arrayname);
   for(int i = 0; i < size; i++){
	   printf("%.4f", UTQ[i]);
	   if(i < size -1)
			printf(", ");
	   else
			printf("];\n");
   }
}

void PrintKuiperPair(KuiperPair sol){
	printf("(%.4f, %.4f)", sol.c, sol.v);
}

int main(int argc, char* argv[])
{
	double   alpha = atof(argv[1]); // 0.05
	unsigned method= atoi(argv[2]); // 1 for Direct and 2 for Newton iterative method

	// CheckParameter(alpha, n, k);	
	
	printf("V_n test, ");
	if( method == 1){
		printf("Directive Iterative Method\n");
	}else{
		printf("Newton's Iterative Method\n");
	}
	
    printf("alpha = %.4f \n", alpha);
	
	
	/* Attention, please!
	 *  For small n, the alpha (probability of type I) should not be
	 *   too small. Otherwise, it is meaningless since we can not get
	 *   precise conclusion with small samples in the sense of science
	 *   and statistics.
	 *  For large n, the alpha can be small since we can get enough 
	 *  knowledge and empirical information from the samples. 
	 *   For example, if n < 10, we should set alpha >= 0.01.
	 *   Moreover, if n >= 60, alpha = 0.001 is OK.
	 */ 
	unsigned kArray[] = {1, 2, 3, 4, 5};
	unsigned nArray[] = {6,  7,   8,    9, 10, 15,  20, 25,  30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85,
		                 90, 95, 100, 180, 200, 250, 300, 400, 500, 550, 560, 600, 700, 800,
		                 900, 1000, 2000, 3000, 4000, 5000,
		                 1000000};
	//unsigned nArray[] = {200, 300, 400, 500, 600, 1000, 5000, 8000, 100000};
	const int N = sizeof(nArray)/sizeof(nArray[0]);
	const int K = sizeof(kArray)/sizeof(kArray[0]);
	KuiperPair Mat[N][K];
	for(int i = 0; i < N; i++){
		for(int j = 0; j < K; j++){
			Mat[i][j] = KuiperPairSolver(
							alpha, 
							nArray[i],  /* n */
							method, 
							kArray[j]   /* k */
							);

		}
	}
	
	printf("\n\n");
	PrintArray1d(kArray, K, "order k");
	for(int i = 0; i < N; i++){
		printf("Sample capacity n = %-6d: ", nArray[i]);
		for(int j = 0; j < K; j++){
			PrintKuiperPair(Mat[i][j]);
			if(j < K -1){printf(", "); }else{ printf("\n");}
		}		
	}
	
	
	PrintArray1d(kArray, K, "order k");
	printf("LaTeX code for table\n\n");
	printf("\\begin{table}[htbp]\n");
	printf("\\centering\n");
	printf("\\caption{Kuiper pair $\\mpair{c^\\alpha_n(k)}{v^\\alpha_n(k)}$");
	printf(" for order $1\\le k \\le 5$ and upper tail probability $\\alpha = %.4f$}\n",alpha);
	printf("\\begin{tabular}{|c|ccccc|}\n");
			printf("\\hline\n");
	printf("\\diagbox{$n$}{$k$}  & $1$ & $2$ &  $3$ &  $4$ &  $5$ \\\\ \n");
	printf("\\hline\n");
	for(int i = 0; i < N; i++){
		printf("$%-6d$ & ", nArray[i]);
		for(int j = 0; j < K; j++){
			printf("$");
			PrintKuiperPair(Mat[i][j]);
			printf("$");
			if(j < K -1){printf("& "); }else{ printf("\\\\ \n");}
		}		
	}
	printf("\\hline\n");
	printf("\\end{tabular}\n");
	printf("\\end{table}\n");
	
     
	/*
	double alpha_list[10] = {0.693, 0.528, 0.377, 0.252, 0.158,
	                         0.093, 0.052, 0.027, 0.0135, 0.0063};
	printf("\n\nTABLE I for the  $V_n$-test by Kuiper 1962: Acomparison:\n");
	KuiperPair sol;
	for(int r= 0; r < K; r++){
	    printf("Order k = %d\n", kArray[r]);
	    printf("critical value c,   quantile v,         alpha \n");
		for(int i = 0; i < 10; i++){
			sol = KuiperPairSolver(alpha_list[i], 10, method, kArray[r]); // n = 10
			printf("\t %.2f, \t\t %.4f, \t %.4f\n", 
		          sol.c, sol.v, alpha_list[i]);
		} 
		printf("\n");
    }
	*/
	
	/*	
	printf("\nKuiper's Upper Tail Quantile and Lower Tail Quantile in V_n test, ");
	printf("with element (utq(a), ltq(1-a))\n");
	const int M1 = 4;
	const int N1 = 7;
	double AlphaUtq[M1]= {0.10, 0.05, 0.02, 0.01};
	double AlphaLtq[M1]= {0.90, 0.95, 0.98, 0.99};
	double utq[M1][N1];
	double ltq[M1][N1];
	unsigned capacity[N1] = {6, 7, 8, 9, 10, 100, 180};
	printf("n = ");
	for(int j = 0; j < N1; j++) printf("%d, ", capacity[j]); 
	printf("\n");
	printf("For alpha_utq + alpha_ltq = 1.0, we have the property utq(a) = ltq(1-a): \n");
	for(int r = 0; r < K; r++){ 
		printf(" k = %d\n", kArray[r]);
		for(int i = 0; i < M1; i++){
			printf("a = %.4f, 1-a = %.4f: ", AlphaUtq[i], AlphaLtq[i]);
			for(int j = 0; j < N1; j++){
				utq[i][j]  = KuiperUTQ(AlphaUtq[i], capacity[j], kArray[r]);
				ltq[i][j]  = KuiperLTQ(AlphaLtq[i], capacity[j], kArray[r]);
				printf("(%.4f, %.4f)", utq[i][j], ltq[i][j]);
				if(j <= N1-2) printf(", ");
			}
			printf("\n");		
		}
		printf("\n\n");
	}

	*/
	
	double Alpha[] = {0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.020,  0.03, 0.04, 0.05};
	unsigned size = sizeof(Alpha)/sizeof(Alpha[0]);
	double UTQste[size];
	for(int i = 0; i < size; i++){
		UTQste[i] = ModiKuiperUTQ(Alpha[i]);
	}
	printf("\n\n");
	printf("Upper tail quantile for Stephens' modified Kuiper statistic Tn:\n");
    PrintUtqSeq(Alpha, size, "utp alpha           ");
    PrintUtqSeq(UTQste,size, "modified Kuiper UTQ ");     
	return 0;
}

