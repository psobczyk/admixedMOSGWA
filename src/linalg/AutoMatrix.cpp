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

#include "AutoMatrix.hpp"
#include <string.h>
#include <assert.h>
#include "package.hpp"

using namespace std;

namespace linalg {

	size_t AutoMatrix::calculateSize ( const size_t rows, const size_t cols ) {
		if ( 0 == rows || 0 == cols ) return 0;	// guard against division by zero in assert below
		const size_t required = rows * cols * sizeof( double );
		assert( required / rows / cols == sizeof( double ) );	// Guard against integer overflow
		return required;
	}

	AutoMatrix::AutoMatrix (
		const size_t rows,
		const size_t cols
	) : Matrix(
		rows,
		cols,
		// if not here, size would be initialised after parent class portion
		static_cast<double*>( malloc( size = calculateSize( rows, cols ) ) ),
		cols
	) {}

	AutoMatrix::AutoMatrix (
		const size_t rows,
		const size_t cols,
		const size_t allocateRows,
		const size_t allocateCols
	) : Matrix(
		rows,
		cols,
		// if not here, size would be initialised after parent class portion
		static_cast<double*>( malloc( size = calculateSize( allocateRows, allocateCols ) ) ),
		allocateCols
	) {}

	AutoMatrix::AutoMatrix ( const AutoMatrix& original ) : Matrix(
		original.countRows(),
		original.countColumns(),
		static_cast<double*>( malloc( original.size ) ),
		original.matrix.tda
	), size( original.size ) {
		const size_t
			rows = countRows(),
			cols = countColumns();
		if ( cols == matrix.tda ) {
			// copy one big lump
			memcpy( matrix.data, original.matrix.data, calculateSize( rows, cols ) );
		} else {
			// copy row by row
			const size_t blockSize = cols * sizeof( double );
			const double *source = original.matrix.data;
			double *target = matrix.data;
			for ( size_t row = 0; row < rows; ++row ) {
				memcpy( target, source, blockSize );
				source += original.matrix.tda;
				target += matrix.tda;
			}
		}
	}

	AutoMatrix& AutoMatrix::operator= ( const AutoMatrix& original ) {
		if ( this != &original ) {
			const size_t
				rows = original.countRows(),
				cols = original.countColumns();
			const size_t used = calculateSize( rows, cols );
			double* newArray;
			if ( used <= size ) {
				newArray = matrix.data;
			} else {
				newArray = static_cast<double*>( malloc( used ) );
				if ( NULL == newArray ) {
					throw bad_alloc();
				}
				free( matrix.data );
				matrix.data = newArray;
				size = used;
			}
			gslInit(
				rows,
				cols,
				newArray,
				calculateSize( rows, original.matrix.tda ) <= size ? original.matrix.tda : cols
				// try to copy tda
			);
			// the following copy is like in the copy constructor
			if ( cols == original.matrix.tda ) {
				// copy one big lump
				memcpy( matrix.data, original.matrix.data, used );
			} else {
				// copy row by row
				const size_t blockSize = cols * sizeof( double );
				const double *source = original.matrix.data;
				double *target = matrix.data;
				for ( size_t row = 0; row < rows; ++row ) {
					memcpy( target, source, blockSize );
					source += original.matrix.tda;
					target += matrix.tda;
				}
			}
		}
		return *this;
	}

	void AutoMatrix::exactSize ( const size_t rows, const size_t cols ) {
		const size_t required = calculateSize( rows, cols );
		double *newArray;
		if ( required == size ) {	// Shortcut: no malloc necessary
			if ( cols <= matrix.tda ) {	// move towards begin
				const size_t currentRows = countRows();
				const size_t blockSize = cols * sizeof( double );
				const double *source = matrix.data;
				double *target = matrix.data;
				for ( size_t row = 0; ++row < currentRows; ) {
					// For row = 0, no copy is needed, so the pre-increment of row and pointers is what we want
					memmove( target += cols, source += matrix.tda, blockSize );
				}
			} else {	// move towards end requires opposite loop direction (new cols bigger, so new rows smaller)
				const size_t currentCols = countColumns();
				const size_t blockSize = currentCols * sizeof( double );
				const double *source = &matrix.data[ matrix.tda * rows ];
				double *target = &matrix.data[ rows * cols ];
				// Guard against unsigned underflow
				if ( 0 < rows ) {
					for ( size_t row = rows; --row >= 1; ) {
						memmove( target -= cols, source -= matrix.tda, blockSize );
					}
				}
			}
			newArray = matrix.data;
		} else {
			newArray = static_cast<double*>( malloc( required ) );
			assert( NULL != newArray );
			const size_t copyRows = min( rows, countRows() );
			const size_t blockSize = min( cols, countColumns() ) * sizeof( double );
			for ( size_t row = 0; row < copyRows; ++row ) {
				memcpy( &newArray[ cols * row ], &matrix.data[ matrix.tda * row ], blockSize );
			}
			free( matrix.data );
			size = required;
		}
		gslInit( rows, cols, newArray, cols );	// update GSL struct variables
	}

	/** The method guarantees that the physical dimensions will not shrink.
	* Moreover, if growth is necessary, it will be in geometrically increasing steps.
	*/
	void AutoMatrix::upSize ( const size_t rows, const size_t cols ) {
		const size_t
			physicalRows = rows <= countRows() ? countRows() : upperPowerOf2( rows ),
			physicalCols = cols <= matrix.tda ? matrix.tda : upperPowerOf2( cols );
		exactSize( physicalRows, physicalCols );
		gslInit( rows, cols, matrix.data, matrix.tda );
	}

	Vector AutoMatrix::newRow () {
		const size_t row = countRows();
		upSize( row + 1, countColumns() );
		return rowVector( row );
	}

	Vector AutoMatrix::newColumn () {
		const size_t col = countColumns();
		upSize( countRows(), col + 1 );
		return columnVector( col );
	}

	AutoMatrix::~AutoMatrix () {
		free( matrix.data );
		matrix.data = NULL;
	}

}
