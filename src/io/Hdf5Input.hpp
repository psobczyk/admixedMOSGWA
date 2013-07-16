/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Bernhard Bodenstorfer.				*
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

#ifndef IO_HDF5INPUT_HPP
#define IO_HDF5INPUT_HPP

#include <string>
#include <hdf5.h>
#include "Input.hpp"
#include "../Exception.hpp"
#include "../linalg/AutoVector.hpp"
#include "../linalg/AutoMatrix.hpp"

namespace io {

	/** Reads input data from HDF5 file format. */
	class Hdf5Input : public Input {

		/** Paths of the relevant objects in the HDF5 file. */
		static const char
			* const genotypeMatrixPath,
			* const individualListPath,
			* const snpListPath,
			* const phenotypeVectorPath;

		/** The path of the opened HDF5-file, used for error reports. */
		const std::string hdf5filename;

		/** HDF5 identifier of the opened HDF5-file. */
		const hid_t hdf5file;

		/** Genome data matrix.
		* It is stored as transposed matrix in order to optimise memory access.
		* The vectors holding all individuals' data for one SNP should be in a contiguous piece of memory
		* so that the whole vector fits into cache.
		* It is these vectors which are added to the regression algorithm in the search phase.
		*/
		linalg::AutoMatrix genotypeMatrixTransposed;

		/** Phenotype data vector. */
		linalg::AutoVector phenotypeVector;

		protected:

		/** Determine length of a one-dimensional array. */
		size_t countDimensions ( const char * const objectPath ) const;

		/** Retrieve a string from the given one-dimensional array of strings.
		* @param objectPath must point to a one-dimensional array of strings.
		*/
		std::string getString ( const char * const objectPath, const size_t index ) const;

		public:

		/** Set up the HDF5 file reading. */
		Hdf5Input ( const char* const filename );

		/** Return the number of SNPs in the data. */
		virtual size_t countSnps () const;

		/** Return the number of individuals in the data. */
		virtual size_t countIndividuals () const;

		/** Retrieve the data for the given SNP. */
		virtual SNP getSnp ( const size_t snpIndex );

		/** Retrieve the data for the given individual. */
		virtual Individual getIndividual ( const size_t individualIndex );

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypesIntoVector ( const size_t snpIndex, linalg::Vector& vector );

		/** Copy the {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypesIntoVector ( linalg::Vector& vector );

		/** Declare access to be finished, release all resources. */
		virtual ~Hdf5Input ();

		private:	// HDF5 helper functions

		hid_t h5fOpen ( const char* const filename ) throw ( Exception );
		void h5fClose () throw ( Exception );

		hid_t h5dOpen ( const char * const objectPath ) const throw ( Exception );
		hid_t h5dType ( const hid_t datasetId, const char * const objectPath ) const throw ( Exception );
		hid_t h5dSpace ( const hid_t datasetId, const char * const objectPath ) const throw ( Exception );
		void h5dReadAll ( const hid_t datasetId, const char * const objectPath, double *buffer ) const throw ( Exception );
		void h5dClose ( const hid_t datasetId, const char * const objectPath ) const throw ( Exception );

		hid_t h5tCopy ( const hid_t datatypeId, const char * const typeDescription ) const throw ( Exception );
		size_t h5tSize ( const hid_t datatypeId, const char * const objectPath ) const throw ( Exception );
		bool h5tIsVarString ( const hid_t datatypeId, const char * const objectPath ) const throw ( Exception );
		void h5tClose ( const hid_t datatypeId, const char * const objectPath ) const throw ( Exception );

		int h5sDims ( const hid_t dataspaceId, const char * const objectPath ) const throw ( Exception );
		void h5sDimSizes ( const hid_t dataspaceId, const char * const objectPath, hsize_t *sizes ) const throw ( Exception );
		void h5sClose ( const hid_t dataspaceId, const char * const objectPath ) const throw ( Exception );

	};

}

#endif	/* IO_HDF5INPUT_HPP */
