/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2014, Bernhard Bodenstorfer.				*
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

#ifndef IO_PLINKINPUT_HPP
#define IO_PLINKINPUT_HPP

#include "InputAdapter.hpp"
#include "../linalg/AutoMatrix.hpp"
#include <fstream>
#include <map>

namespace io {

	/** Reads input data from PLink BIM, FAM and BED formatted files. */
	class PlinkInput : public InputAdapter {

		/** File name for genotype data. */
		std::string genotypeFilename;

		/** Genotype data file. */
		std::ifstream genotypeFile;

		/** Points to genotype array start within the genonotype data file. */
		std::streampos genotypeArrayStart;

		/** Genotype data order. */
		enum { IDV_MAJOUR, SNP_MAJOUR } genotypeArrayTransposition;

		/** Phenotype data vectors for all traits.
		* It is stored as transposed matrix in order to optimise memory access.
		* The vectors holding all individuals' data for one SNP should be in a contiguous piece of memory
		* so that the whole vector fits into cache.
		* It is these vectors which are added to the regression algorithm in the search phase.
		*/
		linalg::AutoMatrix phenotypeMatrixTransposed;

		/** Covariate matrix.
		* Similar to {@link #phenotypeMatrixTransposed},
		* it is stored as transposed matrix in order to optimise memory access.
		*/
		linalg::AutoMatrix covariateMatrixTransposed;

		/** Reads covariate and additional phenotype (yvm) files. */
		void readExtraFile (
			const char * const topic,
			const std::string& filename,
			const std::map<const std::string,size_t>& idvIndex,
			std::vector<std::string>& names,
			linalg::AutoMatrix& matrixTransposed
		);

		public:

		/** Set up the Plink file reading. */
		PlinkInput ( const char* const filename );

		/** Copy the {@link countIndividuals} sized vector of genotype information for the given SNP into the given vector. */
		virtual void retrieveGenotypeVector ( const size_t snpIndex, linalg::Vector& vector );

		/** Copy the {@link countIndividuals} sized vector of phenotype information into the given vector. */
		virtual void retrievePhenotypeVector ( const size_t traitIndex, linalg::Vector& vector );

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void retrieveCovariateVector ( const size_t covIndex, linalg::Vector& vector );

		/** Declare access to be finished, release all resources. */
		virtual ~PlinkInput ();
	};

}

#endif	/* IO_PLINKINPUT_HPP */
