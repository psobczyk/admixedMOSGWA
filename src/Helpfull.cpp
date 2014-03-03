/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
 *										*
 *	This program is free software; you can redistribute it and/or modify	*
 *	it under the terms of the GNU General Public License as published by	*
 *	the Free Software Foundation; either version 3 of the License, or	*
 *	(at your option) any later version.					*
 *										*
 *	This program is distributed in the hope that it will be useful,		*
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of		*
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.			*
 *	See the GNU General Public License for more details.			*
 ********************************************************************************/

#include "Helpfull.hpp"
#include <iostream>

using namespace std;

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// helpful functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/** computes 3 Vectors for the Regression accounting for Population 1,2,3
* @param vecPos: which of the 3 Vectors
* @param individum: the Indivium < Sample number
*/
int populationVector ( const int vecPos, const int individum ) {
	switch ( vecPos ) {
		case 1:
			if ( individum < 45 ) return 0; //CHB
			if ( individum < 90 ) return 0; //JPT
			if ( individum < 150 ) return -1;	//CEU
			if ( individum < 210 ) return 1;	//YRI
			break;
		case 2:
			if ( individum < 45 ) return 1; //CHB
			if ( individum < 90 ) return 0; //JPT
			if ( individum < 150 ) return -1;	//CEU
			if ( individum < 210 ) return 0;	//YRI
			break;
		case 3:
			if ( individum < 45 ) return 0; //CHB
			if ( individum < 90 ) return 1; //JPT
			if ( individum < 150 ) return -1;	//CEU
			if ( individum < 210 ) return 0;	//YRI
			break;
		default: return 0;	// error
	}
	return 0;	// error
}

/** convert string to interger */
int str2int ( const string &str ) {
	stringstream ss(str);
	int n;
	ss >> n;	// TODO<BB>: here and in similar spots: this is not fault-tolerant.
	return n;
}


/** convert interger to string */
string int2str ( const int n ) {
	stringstream ss;
	ss << n;
	return ss.str();
}


/** convert integer to string */
string double2str ( const double x ) {
	stringstream ss;
	ss << x;
	return ss.str();
}


/** convert interger to a formated string (of a given width with filling in the character fillwith) */
string int2strPadWith ( const int n, const int width, const char fillwith = ' ' ) {
	stringstream ss;
	ss << setfill(fillwith) << setw (width) << n;
	return ss.str();
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// mathematical functions
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



int factorial ( int i ) {	// TODO<BB> is integer calculation enough i.e. i less than about 15?
	if ( 0 > i ) throw;	// TODO<BB> handle exceptions
	int prod;
	for ( prod = 1; 1 < i; --i ) {
		prod *= i;	// TODO<BB> make faster using precalculated table
	}
	return prod;
}

double log_factorial ( int i ) {
	if ( 0 > i ) throw;	// TODO<BB> handle exceptions
	double sum;
	for ( sum = 0.0; 1 < i; --i ) {
		sum += log( i );	// TODO<BB> make faster using precalculated table
	}
	return sum;
}

/** computes logistic firth-regression for design matrix x and target 
 plain C code: logistf package for R from Georg Heinze
 All Parameters are needed to be set!!! */
bool logistffit(
				// Input
				// should be k and n 
				const int		*in_k,		// # of rows of X (# of variables) BB: I'd rather call that "columns"
									// BB: Mind that col_fit[] is of length k, not n. Some confusion here.
				const int		*in_n,		// # of columns of X (# of samples) BB: I'd rather call that "rows"
				//no need for pointer here
				//double*
				const gsl_matrix* X,			// design matrix ( n * k )
				const int*         y,			// target vector ( n )
				// Output
				double*			beta,		// regression coefficients (vector k)
				double*			var,		// var matrix ( k * k )
				double*			pi,		// probabilities for target ( n )
				double*			H,			// diag of H matrix vector (n)
				double*			out_loglik, // log-likelihood of the model
				int*			iter,		// # of main iterations;
				int*			evals,		// # of evalations of evaluatios
				double*			lchange,	// the change in the log-likelihood between steps
				double*			ret_max_U_star,
				double*			ret_max_delta,
				// Optional Input
				const double*		weight,		// the weights
				const double*		offset,		// offset
				const int*		firth,		//  use firth-regression
				const int*		col_fit,	// a "boolean" vector indication which colums to use
				const double*		init,		// initial values for beta
				const int*		eval_llh,	// only evaluate likelihood
				// Control Parameters
				const int*		maxit,
				const int*		maxhs,
				const int*		maxstep,
				const double*	lconv,
				const double*	gconv,
				const double*	xconv
			)
{
	int k= *in_k;
	int n= *in_n;
	int i, j,l, ll; //Loop variables
	int s; // needed for computing determinates and inverses of matrices
	int used_col; // the number of columns in col_fit

	// Initialise Matrix X
	//gsl_matrix_view X = gsl_matrix_view_array(x, n, k);

	////////////////////////////////////////////////////////////////////////////
	// Initialize Offset
	gsl_vector_const_view offsetV = gsl_vector_const_view_array(offset, n);

	////////////////////////////////////////////////////////////////////////////
	//// Initialise col_fit
	//// col_fit is boolean vector, indicating if the i-th column is used (true) or not (false)
	//// used_col counts the number of used columns

	used_col=0;
	for (i=0; i < k; i++)
	{
		if (col_fit[i])
		{
			used_col++;
		}
	}
/*used col kann schon außerhalb definiert werden und col_fit ist ja auch außerhalb bekannt.
 * Intern sind diese Dinge constante Größen!
 */
	//////////////////////////////////////////////////////////////////////////////
	//// Initialise Beta

	memcpy(beta, init, sizeof(double)*k);
	gsl_vector_view betaV = gsl_vector_view_array (beta, k);



	//////////////////////////////////////////////////////////////////////////////
	*lchange = 5.0;
	*iter=0;


	// set pi
	gsl_vector * helpVn = gsl_vector_alloc( n );
	gsl_vector_memcpy(helpVn, &offsetV.vector); // helpVn = offsetV
	// TODO<BB>: dubious: here, whole of X is used here, not only the columns marked in col_fit
	gsl_blas_dgemv(CblasNoTrans, -1.0, X, &betaV.vector, -1.0,helpVn);  //helpVn = -1(X*betaV) + (-1)offset

	// reviewer-remark<BB>: mind the rather extensive code-duplication between here and the inner for-loop
	for (i=0; i < n; i++)
	{
		pi[i]=1.0/( 1  + exp (gsl_vector_get(helpVn, i) ) ); //0.5
		//cerr<<"y["<<i<<"]="<<y[i]<<endl;
	}


	// compute log-likelihood
	double loglik=0;
	for (i=0; i< n; i++)
	{
		loglik+= weight[i]*( y[i]*log(pi[i]) + (1-y[i])*log(1-pi[i]) );
		//cerr<<"l="<<weight[i]*( y[i]*log(pi[i]) + (1-y[i])*log(1-pi[i]) )<<endl;//1

	}

	// Initialise Matrices, Vectors and Permutaitons
	gsl_matrix * XW2 = gsl_matrix_alloc (k, n);
	gsl_vector * helpVk = gsl_vector_alloc( k);
	gsl_matrix * Fisher = gsl_matrix_alloc (k, k);
	gsl_permutation * perm_k = gsl_permutation_alloc (k); // needed for LU_decomp
	gsl_matrix * covs = gsl_matrix_alloc (k, k);
	gsl_vector * UStar = gsl_vector_alloc( k );
	gsl_vector * delta = gsl_vector_alloc( k );


	if (*firth)
	{
		// XW2 = X^T W^(1/2) =	crossprod(x, diag(weight * pi * (1 - pi))^0.5)
		for (i=0; i<n; i++)
		{
			gsl_matrix_get_row(helpVk,X, i);//helpVk= X(i,:)
			gsl_blas_dscal(sqrt( weight[i]*pi[i]*( 1-pi[i] ) ),helpVk); //y=\sqrt(w[i].*pi[i].*(1-pi[i]))*helpVk
			gsl_matrix_set_col(XW2, i, helpVk);
		}

		//Fisher Matrix: X^T W X
		//es wird Fisher=+1*XW2*XW2'+0*Fisher
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0, XW2, XW2,0.0, Fisher);

		// to compute determinate, we need an LU-decompositon of Fisher
		s=0;
		gsl_linalg_LU_decomp(Fisher, perm_k, &s); //
		loglik+=  0.5*gsl_linalg_LU_lndet(Fisher);//
		//cerr<<loglik<<endl;

	}

	//Loop;
	*evals = 1;
	int stop=0;
	double loglik_old;


	double mx;

	////////////////////////////////////////////////////////////////////////////
	while (!stop)
	{
		loglik_old = loglik;
		// XW2 = X^T W^(1/2) =	crossprod(x, diag(weight * pi * (1 - pi))^0.5)
		for (i=0; i<n; i++)
		{
			gsl_matrix_get_row(helpVk,X, i);
			gsl_blas_dscal(sqrt( weight[i]*pi[i]*( 1-pi[i] ) ),helpVk);
			gsl_matrix_set_col(XW2, i, helpVk);
		}

		//Fisher Matrix: X^T W X
		gsl_blas_dgemm(CblasNoTrans,CblasTrans,1, XW2, XW2,0, Fisher);

		int s=0;
		gsl_linalg_LU_decomp(Fisher, perm_k, &s);

		// check if matrix is invertible, that is, all diagonal elements are non-zero
		i=0;
		while( i < k)
		{
			if ( fabs( gsl_matrix_get (Fisher, i, i)) <1e-15)
			{
				// Free Memory
				gsl_vector_free (helpVn);
				gsl_matrix_free (XW2);
				gsl_vector_free (helpVk);
				gsl_matrix_free(Fisher);
				gsl_permutation_free(perm_k);
				gsl_matrix_free(covs);
				gsl_vector_free (UStar);
				gsl_vector_free (delta);
				cerr << "Fisher Matrix is not invertible" << endl;
				return false;
			}
			i++;
		}
		// covs =  (X^T W X)^-1
		gsl_linalg_LU_invert(Fisher, perm_k, covs);

		// Compute the diagonal-Elements of (XW2)^T * covs * XW2
		// reviewer-remark<BB>: H is used (unless as return value) only if *firth is true. Put in the if below?
		// reviewer-remark<BB>: Also covs and Fisher are used (unless as return value) only if *firth is true.
		for (i=0; i<n; i++)
		{
			H[i]=0;

			for(j=0; j<k; j++)
			{
				for(l=0; l<k; l++)
				{
					H[i]+=gsl_matrix_get(XW2, j, i)*gsl_matrix_get(covs, j, l)*gsl_matrix_get(XW2, l, i);//covs nur hier gebraucht
				}
			}
		}


		if (*firth)
		{
			for(i=0; i< n; i++)
			{
				gsl_vector_set(helpVn, i, weight[i]*(y[i] - pi[i])+H[i]*(0.5-pi[i])); //nur hier wird H ausgelesen, also nicht im nicht firth Fall
			}
		}
		else
		{
			for(i=0; i< n; i++)
			{
				gsl_vector_set(helpVn, i, weight[i]*(y[i] - pi[i]) );
			}
		}

		// Set UStar
		gsl_blas_dgemv(CblasTrans, 1.0, X, helpVn, 0.0, UStar);

		// allocate XX_covs with 0;
		gsl_matrix * XX_covs = gsl_matrix_calloc (k, k);

		if ( !*eval_llh )
		{
			double diag_element;
			gsl_matrix * XX_XW2 = gsl_matrix_alloc (used_col, n);

			// XX_XW2 = crossprod(x[, col.fit, drop = FALSE], diag(weight * pi * (1 - pi))^0.5)
			// reviewer-remark<BB>: the below code transposes the selected columns
			for (i =0; i < n; i++)
			{
				l=0;
				diag_element = sqrt( weight[i]*pi[i]*(1-pi[i]) );
				for (j = 0; j < k; j++)
				{
					if(col_fit[j]==1)
					{       //gsl vector mal skalar sollte hier gehen
						gsl_matrix_set(XX_XW2, l, i,  diag_element * gsl_matrix_get(X, i, j) );
						l++;
					}
				}
			}

			// Allocate local matrices. dependent on eval_llh, so no earlier allocation possible
			gsl_matrix * XX_Fisher = gsl_matrix_alloc (used_col, used_col);
			gsl_matrix * XX_Fisher_inv = gsl_matrix_alloc (used_col, used_col);
			gsl_permutation * perm_used_col = gsl_permutation_alloc (used_col); // needed for LU_decomp

			//XX_Fisher = (XX_XW2)^T XX_XW2
			//reviewer-remark<BB>: below code transposes the other way round: XX_XW2 (XX_XW2)^T = Xsel^T W Xsel
			gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0, XX_XW2, XX_XW2,0.0, XX_Fisher);

			// invert XX_Fisher, destroys XX_Fisher

			int s=0;
			gsl_linalg_LU_decomp(XX_Fisher, perm_used_col, &s); //


			// check if matrix is invertible, that is, all diagonal elements are non-zero
			i=0;
			while( i < used_col)
			{
				if ( 0.0  == gsl_matrix_get (XX_Fisher, i, i) )
				{

					// Free Memory

					gsl_vector_free (helpVn);
					gsl_matrix_free (XW2);
					gsl_vector_free (helpVk);
					gsl_matrix_free(Fisher);
					gsl_permutation_free(perm_k);
					gsl_matrix_free(covs);
					gsl_vector_free (UStar);
					gsl_vector_free (delta);
					gsl_matrix_free(XX_Fisher);
					gsl_matrix_free(XX_Fisher_inv);
					gsl_matrix_free(XX_XW2);
					gsl_permutation_free(perm_used_col);

					cerr << "error XX_Fish"<<endl;
					return false;
				}
				i++;
			}

			// XX_Fisher_inv= (X^T W X)^-1
			gsl_linalg_LU_invert(XX_Fisher, perm_used_col, XX_Fisher_inv);

			l=0;
			for (i=0; i < k; i++)
			{
				if (col_fit[i]==1)
				{
					ll=0;
					for (j=0; j < k; j++)
					{
						if (col_fit[j]==1)
						{
							gsl_matrix_set(XX_covs, i, j, gsl_matrix_get(XX_Fisher_inv, l, ll));
							ll++;
						}
					}
					l++;
				}
			}

			// Deallocate Memory
			gsl_matrix_free(XX_Fisher);
			gsl_matrix_free(XX_Fisher_inv);
			gsl_matrix_free(XX_XW2);
			gsl_permutation_free(perm_used_col);
		}


		gsl_blas_dgemv(CblasNoTrans, 1.0, XX_covs, UStar, 0.0, delta);

		// return XX_covs as var
		// reviewer-remark<BB>: this return data seems never to be used in MOSGWA
		for (i=0; i < k; i++)
		{
			for (j=0; j < k; j++)
			{
				var[i*k + j]=gsl_matrix_get(XX_covs,i,j);
			}

		}

		// deallocate XX_covs
		gsl_matrix_free(XX_covs);


		// set mx = max(abs(delta))/maxstep
		mx=0;
		for (i = 0; i < k; i++)
		{
			if ( fabs(gsl_vector_get( delta, i ))/(*maxstep) > mx )
			{
				mx = fabs(gsl_vector_get( delta, i ))/(*maxstep);
			}
		}


		if (mx > 1)
		{
				gsl_blas_dscal((1.0/mx), delta);
		}
		(*evals)++;

		if (*maxit > 0)
		{

			(*iter)++;

			// beta= beta + delta
			gsl_blas_daxpy(1.0, delta, &betaV.vector);

			int halfs; //will not be used here only within the loop
			// Half-Steps
			for (halfs=1; halfs <= *maxhs; halfs++)
			{

				// set pi
				gsl_vector_memcpy(helpVn, &offsetV.vector); // helpVn = offsetV
				gsl_blas_dgemv(CblasNoTrans, -1, X, &betaV.vector, -1,helpVn);	//helpVn = -1(X*betaV) + (-1)offset
				for (i=0; i < n; i++)
				{
					pi[i]=1/(1+exp (gsl_vector_get(helpVn, i) ) );//is set to 0.5 in our case
				}

				// compute log-likelihood
				loglik=0;
				for (i=0; i< n; i++)
				{
					loglik+= weight[i]*( y[i]*log(pi[i]) + (1-y[i])*log(1-pi[i]) );
				}

				if (*firth)
				{
					// XW2 = X^T (W ^ 1/2)
					for (i=0; i<n; i++)
					{
						gsl_matrix_get_row(helpVk, X, i);
						gsl_blas_dscal(sqrt( weight[i]*pi[i]*( 1-pi[i] ) ),helpVk);
						gsl_matrix_set_col(XW2, i, helpVk);
					}

					//Fisher Matrix: X^T W X
					gsl_blas_dgemm(CblasNoTrans,CblasTrans,1, XW2, XW2,0, Fisher);

					s=0;
					gsl_linalg_LU_decomp(Fisher, perm_k, &s); //
					loglik= loglik + 0.5*gsl_linalg_LU_lndet(Fisher);// gsl_linalg_LU_lndet computes ln(|det(A)|)
				}

				(*evals)++;
				*lchange = loglik-loglik_old;


				if (loglik > loglik_old)
				{
					break; //	end Half-Step-loop
				}

				for (i=0; i < k; i++) // beta = beta - delta*2^(-halfs)
				{
					beta[i] = beta[i] - ( gsl_vector_get(delta, i) * pow( 2, -halfs ) );
					// TODO<BB>: use daxpy and replace pow( 2, . ) by quick table lookup.
					// and better re-start from original betas than back-correcting beta+t*delta
				}

			}
		}

		//////////////////////////////////////////////////////////////////////////////
		// Stop-Condition:
		if (*iter == *maxit) // maximum Nummer of Steps reached
		{
				stop=1;
		}
		else //or convergence-criterias met
		{
			i=0; // max(|delta|< xconv => |delta_i| <  xconv for all i AND |UStar_i| < gonv for all i
			while( i<k && fabs(gsl_vector_get(delta,i))<= *xconv && fabs(gsl_vector_get(UStar,i))< *gconv  )
			{i++;}

			if (i==k && *lchange < *lconv)
			{
				stop=1;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	// Output-Variables

	*out_loglik=loglik;
	*ret_max_U_star=fabs(gsl_vector_get(UStar,0));
	for (i=0; i<k; i++)
	{
		if (fabs(gsl_vector_get(UStar,i)) > *ret_max_U_star )
		{
			*ret_max_U_star = fabs(gsl_vector_get(UStar,i));
		}
	}
	// ret_max_delta = max(|delta_i|), mx= max(|delta_i|)/maxstep
	*ret_max_delta= mx*(*maxstep);


	// Deallocate Memory
	gsl_vector_free (helpVn);
	gsl_matrix_free (XW2);
	gsl_vector_free (helpVk);
	gsl_matrix_free(Fisher);
	gsl_permutation_free(perm_k);
	gsl_matrix_free(covs);
	gsl_vector_free (UStar);
	gsl_vector_free (delta);

	return true;
}
