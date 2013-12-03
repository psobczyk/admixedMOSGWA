#include "GenotypeFreq.hpp"
#include <gsl/gsl_cdf.h>	// needed to compute p-values for teststatistics

using namespace linalg;

GenotypeFreq::GenotypeFreq( const MData & mData, const int snp ) {
	// case
	r0_=0;
	r1_=0;
	r2_=0;

	// control
	s0_=0;
	s1_=0;
	s2_=0;

	n_ = mData.getIdvNo();
	r_ = mData.getCaseNo();//warning when 0!!
	s_ = mData.getContNo();//warning when 0!!
        if(0==r_){ cerr<<"Warning CaseNo()==0\n";exit (1);}
        if(0==s_){ cerr<<"Warning ContNo()==0\n"; exit(1);}
	// every individual is checked whether it is case or control, and the genotype is counted
	const Vector genotypes = const_cast<MData&>( mData ).getXcolumn( snp );
	const Vector phenotypes = const_cast<MData&>( mData ).getY();
	for ( int idv = 0; idv < n_; ++idv ) {
		if ( parameter.case_value == phenotypes.get( idv ) ) {	// case_value should be 1 normally
			switch( static_cast<int>( genotypes.get( idv ) ) ) {
				case -1: ++r0_; break;
				case 0: ++r1_; break;
				case 1: ++r2_; break;
			}
		} else {	// control
			switch( static_cast<int>( genotypes.get( idv ) ) ) {
				case -1: ++s0_; break;
				case 0: ++s1_; break;
				case 1: ++s2_; break;
			}
		}
	}

	n0_=r0_+s0_;
	n1_=r1_+s1_;
	n2_=r2_+s2_;
}

/** @see GenotypeFreq::ChiSquare() */
double GenotypeFreq::chiSquareSummand (
	const double obs,
	const double totalRow,
	const double totalCol,
	const double total
) const {
	if ( 0 < totalRow && 0 < totalCol ) {
/*		cout<<"obs "<<obs<<endl;
		cout<<"totalRow "<<totalRow<<endl;
		cout<<"totalCol "<<totalCol<<endl;
		cout<<"total "<<total<<endl;*/
		double numerator = totalRow * totalCol / total;
//		cout<<"numerator "<<numerator<<endl;
		double diff = obs - numerator;
//		cout<<"diff "<<diff<<endl;
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
/*
	cout<<"r0_ "<<r0_<<endl;
	cout<<"r1_ "<<r1_<<endl;
	cout<<"r2_ "<<r2_<<endl;
	cout<<"s0_ "<<s0_<<endl;
	cout<<"s1_ "<<s1_<<endl;
	cout<<"s2_ "<<s2_<<endl;
	cout<<"n0_ "<<n0_<<endl;
	cout<<"n1_ "<<n1_<<endl;
	cout<<"n2_ "<<n2_<<endl;
	cout<<"r_ "<<r_<<endl;
	cout<<"s_ "<<s_<<endl;
	cout<<"n_ "<<n_<<endl;

	cout<<"chiSquareSummand(r0_, r_, n0_, n_) "<<chiSquareSummand(r0_, r_, n0_, n_)<<endl;
	cout<<"chiSquareSummand(s0_, s_, n0_, n_) "<<chiSquareSummand(s0_, s_, n0_, n_)<<endl;
	cout<<"chiSquareSummand(r1_, r_, n1_, n_) "<<chiSquareSummand(r1_, r_, n1_, n_)<<endl;
	cout<<"chiSquareSummand(s1_, s_, n1_, n_) "<<chiSquareSummand(s1_, s_, n1_, n_)<<endl;
	cout<<"chiSquareSummand(r2_, r_, n2_, n_) "<<chiSquareSummand(r2_, r_, n2_, n_)<<endl;
	cout<<"chiSquareSummand(s2_, s_, n2_, n_) "<<chiSquareSummand(s2_, s_, n2_, n_)<<endl;
	cout<<"calculateChiSquare"<<chi2<<endl;
*/	
	// return p-value
	return gsl_cdf_chisq_Q( chi2, 2 ); // Q(x) = \int_{x}^{+\infty} p(x') dx'
}

/** Cochran-Armitage trend Test (CATT),
* according to "Robust Test for Single-marker Analysis in Case-Control Genetic Association Studies",
* Qizhai Li et al.
*/
double GenotypeFreq::calculateCATT () const {
	double x[3];
	double theta[3];

	// TODO<BB>: Why are all x[0] == 0? Efficiency!
	switch( parameter.genetic_model ) {	// 1 = Recessive, 2 = Additive, 3 = Dominant
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

	// TODO<BB> check whether double() function is equivalent typecast
	// TODO<BB> simplify sqrt expression.
	const double Zx
		= ( sqrt( double(n_) ) * SumNom )
		/
		sqrt(
			( (double) n_ / r_ + (double) n_ / s_ )
			*
			( SumDenS - SSumDen2 )
		);
	if ( Zx != Zx ) {	// checks if NaN is returned in case of 0/0; TODO<BB> check this idiom
		return 1.0;
	}

	// return p-value
	return gsl_cdf_ugaussian_P(-fabs(Zx)) + gsl_cdf_ugaussian_Q(fabs(Zx));
}
