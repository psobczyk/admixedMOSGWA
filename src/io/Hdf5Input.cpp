#include "Hdf5Input.hpp"
#include "../Exception.hpp"
#include <cstdlib>
#include <cassert>
#include <vector>

using namespace std;
using namespace linalg;

namespace io {

	const char
		* const Hdf5Input::genotypeMatrixPath = "/genome_matrix",
		* const Hdf5Input::individualListPath = "/individuals",
		* const Hdf5Input::snpListPath = "/single_nucleotide_polymorphisms",
		* const Hdf5Input::phenotypeVectorPath = "/phenotypes";

	Hdf5Input::Hdf5Input ( const char* const filename )
		:
		hdf5filename( filename ),
		hdf5file( h5fOpen( filename ) ),
		genotypeMatrixTransposed( 0, 0 ),
		phenotypeVector( 0 )
	{
		{
			// Read genotype matrix
			const hid_t genotypeDataset = h5dOpen( genotypeMatrixPath );
			const hid_t genotypeType = h5dType( genotypeDataset, genotypeMatrixPath );
			const hid_t genotypeSpace = h5dSpace( genotypeDataset, genotypeMatrixPath );
			const int genotypeDim = h5sDims( genotypeSpace, genotypeMatrixPath );
			if( 2 != genotypeDim ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\" expected 2 dimensions but got %d.",
					filename,
					genotypeMatrixPath,
					genotypeDim
				);
			}

			hsize_t dim[2];
			h5sDimSizes( genotypeSpace, genotypeMatrixPath, dim );
			const size_t size = dim[0] * dim[1];
			assert( 0 == dim[0] || dim[1] == size/dim[0] ); // guard against integer overflow
			vector<double> array( size );
			h5dReadAll( genotypeDataset, genotypeMatrixPath, array.data() );
			genotypeMatrixTransposed.exactSize( dim[0], dim[1] );
			genotypeMatrixTransposed.fill( array.data() );

			h5sClose( genotypeSpace, genotypeMatrixPath );
			h5tClose( genotypeType, genotypeMatrixPath );
			h5dClose( genotypeDataset, genotypeMatrixPath );
		}

		{
			// Read phenotype vector
			const hid_t phenotypeDataset = h5dOpen( phenotypeVectorPath );
			const hid_t phenotypeType = h5dType( phenotypeDataset, phenotypeVectorPath );
			const hid_t phenotypeSpace = h5dSpace( phenotypeDataset, phenotypeVectorPath );
			const int phenotypeDim = h5sDims( phenotypeSpace, phenotypeVectorPath );
			if( 1 != phenotypeDim ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\" expected 1 dimension but got %d.",
					filename,
					phenotypeVectorPath,
					phenotypeDim
				);
			}
			hsize_t dim[1];
			h5sDimSizes( phenotypeSpace, phenotypeVectorPath, dim );
			if ( countIndividuals() != dim[0] ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\" second (i.e. storage minor) dimension size is %d,"
					" dataset \"%s\" dimension size is %d, but they should be equal.",
					filename,
					genotypeMatrixPath, countIndividuals(),
					phenotypeVectorPath, dim[0]
				);
			}
			vector<double> array( dim[0] );
			h5dReadAll( phenotypeDataset, phenotypeVectorPath, array.data() );
			phenotypeVector.exactSize( dim[0] );
			phenotypeVector.fill( array.data() );
			h5sClose( phenotypeSpace, phenotypeVectorPath );
			h5tClose( phenotypeType, phenotypeVectorPath );
			h5dClose( phenotypeDataset, phenotypeVectorPath );
		}

		// Test for input data consistency
		if ( countSnps() != countDimensions( snpListPath ) ) {
			throw Exception(
				"HDF5 input file \"%s\" lists %l SNPs in dataset \"%s\", but would require %l in dataset \"%s\".",
				hdf5filename.c_str(),
				countDimensions( snpListPath ),
				snpListPath,
				countSnps(),
				genotypeMatrixPath
			);
		}
		if( countIndividuals() != countDimensions( individualListPath ) ) {
			throw Exception(
				"HDF5 input file \"%s\" lists %l individuals in dataset \"%s\", but would require %l in dataset \"%s\".",
				hdf5filename.c_str(),
				countDimensions( individualListPath ),
				individualListPath,
				countIndividuals(),
				genotypeMatrixPath
			);
		}
	}

	size_t Hdf5Input::countDimensions ( const char * const objectPath ) const {
		const hid_t datasetId = h5dOpen( objectPath );
		const hid_t dataspaceId = h5dSpace( datasetId, objectPath );
		const int dimensions = h5sDims( dataspaceId, objectPath );
		if( 1 != dimensions ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" expected 1 dimension but got %d.",
				hdf5filename.c_str(),
				objectPath,
				dimensions
			);
		}
		hsize_t dim[1];
		h5sDimSizes( dataspaceId, objectPath, dim );

		h5sClose( dataspaceId, objectPath );
		h5dClose( datasetId, objectPath );

		assert( dim[0] == (size_t) dim[0] );	// Guard against overflowingly large sizes

		return dim[0];
	}

	string Hdf5Input::getString ( const char * const objectPath, const size_t index ) const {
		const hid_t datasetId = h5dOpen( objectPath );

		const hid_t datatypeId = h5dType( datasetId, objectPath );
		const size_t datatypeSize = h5tSize( datatypeId, objectPath ) + 1;	// + 1 for trailing \000
		const bool isVarString = h5tIsVarString( datatypeId, objectPath );
		const hid_t memType = h5tCopy( H5T_C_S1, "H5T_C_S1" );
		if ( 0 > H5Tset_size( memType, isVarString ? H5T_VARIABLE : datatypeSize ) ) {
			throw isVarString
				? Exception( "HDF5 set string size %l failed.", datatypeSize )
				: Exception( "HDF5 set string size variable failed." );
		}

		const hid_t dataspaceId = h5dSpace( datasetId, objectPath );
		const int dimensions = h5sDims( dataspaceId, objectPath );
		if( 1 != dimensions ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" expected 1 dimension but got %d.",
				hdf5filename.c_str(),
				objectPath,
				dimensions
			);
		}
		hsize_t dim[1];
		h5sDimSizes( dataspaceId, objectPath, dim );
		assert( index < dim[0] );	// requested index in range?
		const hsize_t coordinates[1] = { index };
		if ( 0 > H5Sselect_elements( dataspaceId, H5S_SELECT_SET, 1, coordinates ) ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" select element %l failed.",
				hdf5filename.c_str(),
				objectPath,
				index
			);
		}
		const hid_t memSpace = H5Screate( H5S_SCALAR );
		if ( 0 > memSpace ) {
			throw Exception( "HDF5 create scalar space (for memory) failed." );
		}

		string value;
		if ( isVarString ) {
			char* buffer[1];
			if ( 0 > H5Dread( datasetId, memType, memSpace, dataspaceId, H5P_DEFAULT, buffer ) ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\" read variable length string[%l] failed.",
					hdf5filename.c_str(),
					objectPath,
					index
				);
			}
			value = string( buffer[0] );
			if ( 0 > H5Dvlen_reclaim( memType, memSpace, H5P_DEFAULT, buffer ) ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\" read variable length string[%l] failed.",
					hdf5filename.c_str(),
					objectPath,
					index
				);
			}
		} else {
			vector<char> buffer( datatypeSize );
			if ( 0 > H5Dread( datasetId, memType, memSpace, dataspaceId, H5P_DEFAULT, buffer.data() ) ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\" read fixed length string[%l] failed.",
					hdf5filename.c_str(),
					objectPath,
					index
				);
			}
			value = string( buffer.data(), 0, datatypeSize );
		}
		// TODO: handle erddror return values
		// TODO: make sure a trailing \000 is written in all cases.

		h5sClose( memSpace, "[memory]" );
		h5tClose( memType, "[memory]" );
		h5sClose( dataspaceId, objectPath );
		h5tClose( datatypeId, objectPath );
		h5dClose( datasetId, objectPath );

		return value;
	}

	size_t Hdf5Input::countSnps () const {
		return genotypeMatrixTransposed.countRows();
	}

	size_t Hdf5Input::countIndividuals () const {
		return genotypeMatrixTransposed.countColumns();
	}

	SNP Hdf5Input::getSnp ( const size_t snpIndex ) {
		const string snpId = getString( snpListPath, snpIndex );
		const size_t idLength = snpId.length();
		size_t
			i,
			chromosomeIdLength = 0,
			positionStringStart = 0;
		for ( i = 0; i < idLength; ++i ) {
			const char c = snpId[i];
			if ( '_' == c && 0 == positionStringStart ) {
				chromosomeIdLength = i;
				positionStringStart = i + 1;
			} else if ( c < '0' || '9' < c ) {
				throw Exception(
					"HDF5 input file \"%s\" dataset \"%s\""
					" SNP[%d] has bad character in name \"%s\" at position %l;"
					" expecting decimals indicating chromosome and position,"
					" separated by a single underscore.",
					hdf5filename.c_str(),
					snpListPath,
					snpIndex,
					snpId.c_str(),
					i
				);
			}
		}
		const string chromosomeId( snpId, 0, chromosomeIdLength );
		const unsigned long position = strtoul( snpId.c_str() + positionStringStart, NULL, 10 );
		return SNP( chromosomeId, snpId, 0.0, position, 0, 0 );
	}

	Individual Hdf5Input::getIndividual ( const size_t individualIndex ) {
		const string individualId = getString( individualListPath, individualIndex );

		const Individual individual(
			"",
			individualId,
			"",
			"",
			Individual::MISSING,
			getPhenotypeVector().get( individualIndex )
		);

		return individual;
	}

	Vector Hdf5Input::getGenotypeVector ( const size_t snpIndex ) {
		return genotypeMatrixTransposed.rowVector( snpIndex );
	}

	Vector Hdf5Input::getPhenotypeVector () {
		return phenotypeVector;
	}

	Hdf5Input::~Hdf5Input () {
		h5fClose();
	}

	hid_t Hdf5Input::h5fOpen ( const char * const filename ) throw ( Exception ) {
		const hid_t fileId = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
		if ( 0 > fileId ) {
			throw Exception(
				"Open HDF5 input file \"%s\" failed.",
				filename
			);
		}
		return fileId;
	}

	void Hdf5Input::h5fClose () throw ( Exception ) {
		if ( 0 > H5Fclose( hdf5file ) ) {
			throw Exception(
				"HDF5 input file \"%s\" close failed.",
				hdf5filename.c_str()
			);
		}
	}

	hid_t Hdf5Input::h5dOpen ( const char * const objectPath ) const throw ( Exception ) {
		const hid_t datasetId = H5Dopen2( hdf5file, objectPath, H5P_DEFAULT );
		if ( 0 > datasetId ) {
			throw Exception(
				"HDF5 input file \"%s\" does not contain dataset \"%s\".",
				hdf5filename.c_str(),
				objectPath
			);
		}
		return datasetId;
	}

	hid_t Hdf5Input::h5dType ( const hid_t datasetId, const char * const objectPath ) const throw ( Exception ) {
		const hid_t datatypeId = H5Dget_type( datasetId );
		if ( 0 > datatypeId ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" failed to determine datatype.",
				hdf5filename.c_str(),
				objectPath
			);
		}
		return datatypeId;
	}

	hid_t Hdf5Input::h5dSpace ( const hid_t datasetId, const char * const objectPath ) const throw ( Exception ) {
		const hid_t dataspaceId = H5Dget_space( datasetId );
		if ( 0 > dataspaceId ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" failed to determine dataspace.",
				hdf5filename.c_str(),
				objectPath
			);
		}
		return dataspaceId;
	}

	void Hdf5Input::h5dReadAll ( const hid_t datasetId, const char * const objectPath, double *buffer ) const throw ( Exception ) {
		if ( 0 > H5Dread( datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer ) ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" read failed.",
				hdf5filename.c_str(),
				objectPath
			);
		}
	}

	void Hdf5Input::h5dClose ( const hid_t datasetId, const char * const objectPath ) const throw ( Exception ) {
		if ( 0 > H5Dclose( datasetId ) ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" close failed.",
				hdf5filename.c_str(),
				objectPath
			);
		}
	}

	hid_t Hdf5Input::h5tCopy ( const hid_t datatypeId, const char * const typeDescription ) const throw ( Exception ) {
		const hid_t copyId = H5Tcopy( H5T_C_S1 );
		if ( 0 > copyId ) {
			throw Exception(
				"HDF5 input file \"%s\" type \"%s\" copy failed.",
				hdf5filename.c_str(),
				typeDescription
			);
		}
		return copyId;
	}

	size_t Hdf5Input::h5tSize ( const hid_t datatypeId, const char * const objectPath ) const throw ( Exception ) {
		const size_t datatypeSize = H5Tget_size( datatypeId );
		if ( 0 >= datatypeSize ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" failed to get datatype size.",
				hdf5filename.c_str(),
				objectPath
			);
		}
		return datatypeSize;
	}

	bool Hdf5Input::h5tIsVarString ( const hid_t datatypeId, const char * const objectPath ) const throw ( Exception ) {
		const htri_t isVarString = H5Tis_variable_str( datatypeId );
		if ( 0 > isVarString ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" failed to determine whether datatype is a variable length string.",
				hdf5filename.c_str(),
				objectPath
			);
		}
		return isVarString;
	}

	void Hdf5Input::h5tClose ( const hid_t datatypeId, const char * const objectPath ) const throw ( Exception ) {
		if ( 0 > H5Tclose( datatypeId ) ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" datatype close failed.",
				hdf5filename.c_str(),
				objectPath
			);
		}
	}

	int Hdf5Input::h5sDims ( const hid_t dataspaceId, const char * const objectPath ) const throw ( Exception ) {
		const int dimensions = H5Sget_simple_extent_ndims( dataspaceId );
		if ( 0 > dimensions ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" failed to determine dataspace dimension.",
				hdf5filename.c_str(),
				objectPath
			);
		}
		return dimensions;
	}

	void Hdf5Input::h5sDimSizes ( const hid_t dataspaceId, const char * const objectPath, hsize_t *sizes ) const throw ( Exception ) {
		if ( 0 > H5Sget_simple_extent_dims( dataspaceId, sizes, NULL ) ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" failed to determine dataspace dimension sizes.",
				hdf5filename.c_str(),
				objectPath
			);
		}
	}

	void Hdf5Input::h5sClose ( const hid_t dataspaceId, const char * const objectPath ) const throw ( Exception ) {
		if ( 0 > H5Sclose( dataspaceId ) ) {
			throw Exception(
				"HDF5 input file \"%s\" dataset \"%s\" dataspace close failed.",
				hdf5filename.c_str(),
				objectPath
			);
		}
	}

}
