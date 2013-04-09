#ifndef PACKAGE_IMPUTATION_HPP
#define PACKAGE_IMPUTATION_HPP

#include "../MData.hpp"

/** Imputation algorithms. These try to find plausible data to replace missing entries in {@link MData}. */
namespace imputation {

	/** Perform the imputation. */
	void imputate ( MData& mData );

}

#endif	/* PACKAGE_IMPUTATION_HPP */
