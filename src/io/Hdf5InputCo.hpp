/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2013, Bernhard Bodenstorfer.					*
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

#ifndef IO_HDF5INPUTCO_HPP
#define IO_HDF5INPUTCO_HPP

#include "InputCo.hpp"
#include "Hdf5Input.hpp"

namespace io {

	/** Reads input data from HDF5 file format with covariates. */
	class Hdf5InputCo : virtual public Hdf5Input, virtual public InputCo {

		protected:

		/** Paths of the relevant objects in the HDF5 file. */
		static const char
			* const covariateListPath,
			* const covariateMatrixPath;

		/** Gets covariate names into MOSGWA. */
		Hdf5StringList covariateNames;

		/** Gets covariate matrix into MOSGWA. */
		Hdf5DoubleTable covariatesTransposed;

		public:

		/** Set up the HDF5 file reading. */
		Hdf5InputCo ( const char * const filename );

		/** Return the number of covariate vectors in the data. */
		virtual size_t countCovariateVectors ();

		/** Get the name of the given covariate. */
		virtual std::string getCovariateName ( const size_t covIndex );

		/** Copy the {@link countIndividuals} sized vector of covariate information for the given covariate into the given vector. */
		virtual void retrieveCovariatesIntoVector ( const size_t covIndex, linalg::Vector& vector );

		/** Declare access to be finished, release all resources. */
		virtual ~Hdf5InputCo ();
	};

}

#endif	/* IO_HDF5INPUTCO_HPP */
