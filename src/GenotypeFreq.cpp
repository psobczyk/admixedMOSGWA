/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2015, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#include "GenotypeFreq.hpp"
#include "logging/Logger.hpp"
#include <gsl/gsl_cdf.h>	// needed to compute p-values for teststatistics

using namespace linalg;
using namespace logging;

GenotypeFreq::GenotypeFreq(
	const MData & mData,
	const size_t snp
) : n_( mData.getIdvNo() ) {
	r0_ = r1_ = r2_ = s0_ = s1_ = s2_ = 0u;

	AutoVector
		genotypes( n_ ),
		phenotypes( n_ );

	mData.getXcolumn( snp, genotypes );
	mData.getY( phenotypes );

	if ( parameter->case_value == parameter->control_value ) {
		logger->error(
			"Case and control individuals have the same phenotype value %f==%f",
			parameter->case_value,
			parameter->control_value
		);
	}

	if ( parameter->case_value == parameter->missing_phenotype_value ) {
		logger->error(
			"Case individuals have the same value as missing phenotype %f==%f",
			parameter->case_value,
			parameter->missing_phenotype_value
		);
	}

	if ( parameter->control_value == parameter->missing_phenotype_value ) {
		logger->error(
			"Control individuals have the same value as missing phenotype %f==%f",
			parameter->control_value,
			parameter->missing_phenotype_value
		);
	}

	// TODO<BB>: some values are the same for all SNPs, so do not calculate twice.
	for ( size_t idv = 0; idv < n_; ++idv ) {
		const double y = phenotypes.get( idv );
		if ( y == parameter->case_value ) {
			switch( static_cast<int>( genotypes.get( idv ) ) ) {
				case -1: ++r0_; break;
				case 0: ++r1_; break;
				case 1: ++r2_; break;
			}
		} else if ( y == parameter->control_value ) {
			switch( static_cast<int>( genotypes.get( idv ) ) ) {
				case -1: ++s0_; break;
				case 0: ++s1_; break;
				case 1: ++s2_; break;
			}
		} else {
			logger->error(
				"Phenotype[%u] %f is neither case nor control. Statistical tests will be inaccurate.",
				idv,
				y
			);
		}
	}

	n0_ = r0_ + s0_;
	n1_ = r1_ + s1_;
	n2_ = r2_ + s2_;
	r_ = r0_ + r1_ + r2_;
	s_ = s0_ + s1_ + s2_;

	if ( 0 == r_ ) {
		logger->warning( "No individuals in the case group." );
	}
	if ( 0 == s_ ) {
		logger->warning( "No individuals in the control group." );
	}
}

/** @see GenotypeFreq::ChiSquare() */
double GenotypeFreq::chiSquareSummand (
	const double obs,
	const double totalRow,
	const double totalCol,
	const double total
) const {
	if ( 0 < totalRow && 0 < totalCol ) {
		double numerator = totalRow * totalCol / total;
		double diff = obs - numerator;
		return
			diff * diff / numerator;
	} else {
		return 0.0;
	}
}

/**
	\f[ \chi^2 = \sum_{i=0}^2 [ (r_i - r*n_i/n)^2 /(r*n_i/n) + (s_i - s*n_i/n)^2 /(s*n_i/n)] \f]
	n > 0, check if r*n_i, s*n_i > 0 -> compute summand, else summand = 0
	chiSquareSummand computes this.
	Check correctness: http://www.people.ku.edu/~preacher/chisq/chisq.htm
*/
double GenotypeFreq::calculateChiSquare () const {
	double chi2 = chiSquareSummand(r0_, r_, n0_, n_)
		+ chiSquareSummand(s0_, s_, n0_, n_)
		+ chiSquareSummand(r1_, r_, n1_, n_)
		+ chiSquareSummand(s1_, s_, n1_, n_)
		+ chiSquareSummand(r2_, r_, n2_, n_)
		+ chiSquareSummand(s2_, s_, n2_, n_);
	// return p-value
	// return gsl_cdf_chisq_Q( chi2, 2 ); // Q(x) = \int_{x}^{+\infty} p(x') dx'
	return exp( chi2 * -0.5 );	// BB: in case of 2 degrees of freedom, chisquare greatly simplifies
}

/** Cochran-Armitage trend Test (CATT),
* according to "Robust Test for Single-marker Analysis in Case-Control Genetic Association Studies",
* Qizhai Li et al.
*/
double GenotypeFreq::calculateCATT () const {
	double x[3];
	double theta[3];

	// TODO<BB>: Why are all x[0] == 0? Efficiency!
	switch( parameter->genetic_model ) {	// 1 = Recessive, 2 = Additive, 3 = Dominant
		case 1: // 1 means Recessive
			x[0] = 0.0; x[1] = 0.0; x[2] = 1.0;
			break;
		case 3: // 3 means Dominant
			x[0] = 0.0; x[1] = 1.0; x[2] = 1.0;
			break;
		default:
			// 2 means Additive, default Model
			x[0] = 0.0; x[1] = 0.5; x[2] = 1.0;
			break;
	}

	// Caution! all values are integers,
	theta[0] = (double) n0_ / n_;
	theta[1] = (double) n1_ / n_;
	theta[2] = (double) n2_ / n_;

	const double SumNom
		= x[0] * ( (double) r0_ / r_ - (double) s0_ / s_)
		+ x[1] * ( (double) r1_ / r_ - (double) s1_ / s_)
		+ x[2] * ( (double) r2_ / r_ - (double) s2_ / s_);
	const double SumDenS
		= x[0] * x[0] * theta[0]
		+ x[1] * x[1] * theta[1]
		+ x[2] * x[2] * theta[2];
	const double SSumDen
		= x[0] * theta[0]
		+ x[1] * theta[1]
		+ x[2] * theta[2];
	const double SSumDen2 = SSumDen * SSumDen;

	// TODO<BB> simplify sqrt expression.
	// As to "C-like" casts: http://cs.txstate.edu/~br02/cs1428/SupportFiles/Programming/TypeCasting.htm
	const double Zx = SumNom * sqrt(
		n_ / (
			(
				static_cast<double>( n_ ) / r_
				+
				static_cast<double>( n_ ) / s_
			) * ( SumDenS - SSumDen2 )
		)
	);
	if ( Zx != Zx ) {	// checks if NaN is returned in case of 0/0; TODO<BB> check this idiom
		                // look at: http://stackoverflow.com/questions/570669/checking-if-a-double-or-float-is-nan-in-c
		return 1.0;
	}

	// return p-value
	return 2 * gsl_cdf_ugaussian_P(-fabs(Zx));
}
