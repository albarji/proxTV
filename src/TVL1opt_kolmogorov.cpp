/*
    Extracted from the source code in http://pub.ist.ac.at/~vnk/papers/TV_TREE.html
    with the implementation of the method in the paper

    Total Variation on a Tree
    by Vladimir Kolmogorov, Thomas Pock and Michal Rolínek. 

    Minor modifications by Álvaro Barbero.
*/
#include <stdio.h>
#include <stdlib.h>
#include "TVopt.h"

/********************************************
    Data Structures
********************************************/
struct Breakpoint_a1
{
	double lambda;
	int _slope; // "_" means that it's not the true slope
};


///////////////////////////////////////////////////////////////////////////////////////////
//                    1. Convex quadratic unary terms + TV regularizer                   //
///////////////////////////////////////////////////////////////////////////////////////////

//       unary terms f_i(z):        = 1/2 z*z - b_i z     with a_i>0
//       pairwise terms f_{ij}(z):  = w_{ij} |z|

// Input: a, b: arrays of size n
//        w: array of size n-1
// Output: solution (this array of size n must be allocated by the user)
//
// Complexity: O(n)


void SolveTVConvexQuadratic_a1(int n, double* b, double* w, double* solution)
{
	if (n <= 1)
	{
		if (n == 1) solution[0] = b[0];
		return;
	}

	int i;
	Breakpoint_a1* P = new Breakpoint_a1[2*n-1];
	Breakpoint_a1* L = &P[n-1];
	Breakpoint_a1* R = L+1;
	double* lambda_min = solution;
	double* lambda_max = new double[n-1];

	int A = 1;


	// message at node 1
	(L-1)->_slope = -A;
	L->lambda = lambda_min[0] = (-w[0] + b[0]) / A;
	L->_slope = 0;
	R->lambda = lambda_max[0] = ( w[0] + b[0]) / A;
	R->_slope = -A;

	for (i=1; ; i++)
	{
		A ++;

		int slope;
		double W_prev = (i>0) ? w[i-1] : 0;
		double W = (i < n-1) ? w[i] : 0;
		double message_min = -W_prev + 1*L->lambda - b[i];
		double message_max =  W_prev + 1*R->lambda - b[i];

		// lower bound
		slope = 1;
		while ( message_min < -W )
		{
			slope = L->_slope + A;
			L ++;
			if (L > R) break;
			message_min += (L->lambda - (L-1)->lambda) * slope;
		}

		if (i == n-1)
		{
			if (L > R) L --;
			solution[i] = L->lambda - message_min / slope;
			break;
		}

		L --;
		(L-1)->_slope = -A;

		if (L == R)
		{
			R ++;

			R->_slope = -A;
			R->lambda = lambda_max[i] = L->lambda - (message_max - W) / 1;
			L->lambda = lambda_min[i] = L->lambda - (message_max + W) / 1;

			continue;
		}

		L->lambda = lambda_min[i] = (L+1)->lambda - (W + message_min) / slope;

		// upper bound
		slope = 1;
		while ( message_max > W )
		{
			R --;
			slope = R->_slope + A;
			message_max -= ((R+1)->lambda - R->lambda) * slope;
			if (R == L) break;
		}

		R ++;
		R->_slope = -A;
		R->lambda = lambda_max[i] = (R-1)->lambda + (W - message_max) / slope;
	}

	for (i-- ; i>=0; i--)
	{
		if (solution[i+1] < lambda_min[i]) solution[i] = lambda_min[i];
		else if (solution[i+1] > lambda_max[i]) solution[i] = lambda_max[i];
		else solution[i] = solution[i+1];
	}

	delete [] P;
	delete [] lambda_max;
}

// non-weighted version
void SolveTVConvexQuadratic_a1_nw(int n, double* b, double w, double* solution)
{
	if (n <= 1)
	{
		if (n == 1) solution[0] = b[0];
		return;
	}

	int i;
	Breakpoint_a1* P = new Breakpoint_a1[2*n-1];
	Breakpoint_a1* L = &P[n-1];
	Breakpoint_a1* R = L+1;
	double* lambda_min = solution;
	double* lambda_max = new double[n-1];

	int A = 1;


	// message at node 1
	(L-1)->_slope = -A;
	L->lambda = lambda_min[0] = (-w + b[0]) / A;
	L->_slope = 0;
	R->lambda = lambda_max[0] = ( w + b[0]) / A;
	R->_slope = -A;

	for (i=1; ; i++)
	{
		A ++;

		int slope;
		double W_prev = (i>0) ? w : 0;
		double W = (i < n-1) ? w : 0;
		double message_min = -W_prev + 1*L->lambda - b[i];
		double message_max =  W_prev + 1*R->lambda - b[i];

		// lower bound
		slope = 1;
		while ( message_min < -W )
		{
			slope = L->_slope + A;
			L ++;
			if (L > R) break;
			message_min += (L->lambda - (L-1)->lambda) * slope;
		}

		if (i == n-1)
		{
			if (L > R) L --;
			solution[i] = L->lambda - message_min / slope;
			break;
		}

		L --;
		(L-1)->_slope = -A;

		if (L == R)
		{
			R ++;

			R->_slope = -A;
			R->lambda = lambda_max[i] = L->lambda - (message_max - W) / 1;
			L->lambda = lambda_min[i] = L->lambda - (message_max + W) / 1;

			continue;
		}

		L->lambda = lambda_min[i] = (L+1)->lambda - (W + message_min) / slope;

		// upper bound
		slope = 1;
		while ( message_max > W )
		{
			R --;
			slope = R->_slope + A;
			message_max -= ((R+1)->lambda - R->lambda) * slope;
			if (R == L) break;
		}

		R ++;
		R->_slope = -A;
		R->lambda = lambda_max[i] = (R-1)->lambda + (W - message_max) / slope;
	}

	for (i-- ; i>=0; i--)
	{
		if (solution[i+1] < lambda_min[i]) solution[i] = lambda_min[i];
		else if (solution[i+1] > lambda_max[i]) solution[i] = lambda_max[i];
		else solution[i] = solution[i+1];
	}

	delete [] P;
	delete [] lambda_max;
}
