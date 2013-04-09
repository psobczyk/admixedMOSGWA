#ifndef IO_INPUT_HPP
#define IO_INPUT_HPP

#include "../SNP.hpp"
#include "../Individual.hpp"
#include "../linalg/Vector.hpp"

/** Provides input/output exchange of regression data.
* @author Bernhard Bodenstorfer
*/
namespace io {

	/** Interface to provide the regression data input. */
	struct Input {

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () const = 0;

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () const = 0;

		/** Retrieve the data for the given SNP. */
		virtual SNP getSnp ( const size_t snpIndex ) = 0;

		/** Retrieve the data for the given individual. */
		virtual Individual getIndividual ( const size_t individualIndex ) = 0;

		/** Retrieve the {@link countIndividuals} sized vector of genotype information for the given SNP. */
		virtual linalg::Vector getGenotypeVector ( const size_t snpIndex ) = 0;

		/** Retrieve the {@link countIndividuals} sized vector of phenotype information. */
		virtual linalg::Vector getPhenotypeVector () = 0;

		/** Declare access to be finished, release all resources. */
		virtual ~Input ();
	};

}

#endif	/* IO_INPUT_HPP */
