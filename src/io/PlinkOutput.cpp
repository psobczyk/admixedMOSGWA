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

#include "PlinkOutput.hpp"
#include "PlinkConstants.hpp"
#include "../Exception.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

using namespace std;
using namespace linalg;
using namespace io::PlinkConstants;

namespace io {

	const char PlinkOutput::fillBytes[] = {
		0, 0x1, 0x5, 0x15, 0x55
	};

	PlinkOutput::PlinkOutput (
		const char * const filenameTrunk,
		const size_t individualCount,
		const size_t snpCount,
		const size_t covariateCount,
		const size_t traitCount
	)
	:
		OutputAdapter( individualCount, snpCount, covariateCount, traitCount ),
		filenameTrunk( filenameTrunk ),
		genotypeStream( ( this->filenameTrunk + genotypeMatrixExtension ).c_str(), ios::binary ),
		covariateMatrixTransposed( covariateCount, individualCount ),
		phenotypeMatrixTransposed( traitCount, individualCount )
	{
		genotypeStream.write( bedFileMagic, sizeof( bedFileMagic ) );
		genotypeOrigin = genotypeStream.tellp();
		for ( size_t snp = 0; snp < snpCount; ++snp ) {
			for ( size_t idv = 0; idv < individualCount; ) {
				if ( idv+4 <= individualCount ) {
					genotypeStream.write( fillBytes + 4, 1 );
					idv += 4;
				} else {
					// fillBytes[0] is never used, but given for clarity
					genotypeStream.write( fillBytes + individualCount - idv, 1 );
					break;
				}
			}
		}
	}

	void PlinkOutput::setSnps ( const SNP * snps ) {
		const size_t snpCount = countSnps();
		ofstream snpStream( ( this->filenameTrunk + snpListExtension ).c_str() );
		for ( size_t snpIndex = 0; snpIndex < snpCount; ++snpIndex ) {
			const SNP& snp = snps[snpIndex];
			snpStream
				<< snp.getChromosome()
				<< " "
				<< snp.getSnpId()
				<< "\t"
				<< snp.getGeneticDistance()
				<< "\t"
				<< snp.getBasePairPosition()
				<< "\t"
				<< snp.getAllele1()
				<< "\t"
				<< snp.getAllele2()
				<< endl;
		}
		snpStream.close();
	}

	/** This method is not thread-safe, since the stream is a shared resource. */
	void PlinkOutput::storeGenotypeVector ( const size_t snpIndex, const Vector& v ) {
		const size_t
			idvCount = countIndividuals(),
			bytePerSnp = ( idvCount >> 2 ) + ( 0x3 & idvCount ? 1 : 0 ),
			position = snpIndex * bytePerSnp;
		assert( 0 == idvCount || snpIndex == position / bytePerSnp );	// Guard against overflow
		assert( v.countDimensions() == idvCount );
		assert( snpIndex < countSnps() );
		streampos snpOrigin = genotypeOrigin + streamoff( position );
		genotypeStream.seekp( snpOrigin );	// this implies that the method is not thread-safe
		unsigned int count = 0;
		char accumulator = 0;
		for ( size_t idv = 0; idv < idvCount; ++idv ) {
			const double x = v.get( idv );
			unsigned int pattern;
			for (
				pattern = 0;
				pattern < sizeof( genotypeTranslation ) / sizeof( genotypeTranslation[0] );
				++pattern
			) {
				if ( genotypeTranslation[pattern] == x ) {
					break;
				}
			}
			if ( sizeof( genotypeTranslation ) / sizeof( genotypeTranslation[0] ) <= pattern ) {
				throw Exception(
					"Cannot convert genotype for individual %u and SNP %u"
					" value of %f to BED file."
					" Only values -1, 0, +1 and NaN are allowed.",
					idv,
					snpIndex,
					x
				);
			}
			accumulator |= pattern << count;
			if ( 8 <= ( count+= 2 ) ) {
				genotypeStream.write( &accumulator, sizeof( accumulator ) );
				count = accumulator = 0;
			}
		}
		if ( count ) {
			genotypeStream.write( &accumulator, sizeof( accumulator ) );
		}
	}

	void PlinkOutput::storeCovariateVector ( const size_t covIndex, const Vector& v ) {
		Vector covariateVector = covariateMatrixTransposed.rowVector( covIndex );
		covariateVector.copy( v );
	}

	void PlinkOutput::storePhenotypeVector ( const size_t traitIndex, const Vector& v ) {
		Vector phenotypeVector = phenotypeMatrixTransposed.rowVector( traitIndex );
		phenotypeVector.copy( v );
	}

	PlinkOutput::~PlinkOutput () {
		genotypeStream.close();
		const string
			famPath = this->filenameTrunk + individualListExtension,
			covPath = this->filenameTrunk + covariateMatrixExtension,
			yvmPath = this->filenameTrunk + phenotypeMatrixExtension;
		const size_t
			idvCount = countIndividuals(),
			covCount = countCovariates(),
			traitCount = countTraits();
		const Individual * individuals = getIndividuals();
		{
			ofstream famStream( famPath.c_str() );
			for ( size_t idv = 0; idv < idvCount; ++idv ) {
				const Individual& individual = individuals[idv];
				famStream
					<< individual.getFamilyID()
					<< " "
					<< individual.getIndividualID()
					<< "\t"
					<< individual.getPaternalID()
					<< "\t"
					<< individual.getMaternalID()
					<< "\t"
					<< individual.getSexCode();
				if ( 0 < traitCount ) {
					famStream
						<< "\t"
						<< phenotypeMatrixTransposed.get( 0, idv );
				}
				famStream << endl;
			}
			famStream.close();
		}
		if ( 0 < covCount ) {
			const string * covariates = getCovariates();
			ofstream covStream( covPath.c_str() );
			covStream << "FID IID";
			for ( size_t cov = 0; cov < covCount; ++cov ) {
				covStream
					<< "\t"
					<< covariates[ cov ].c_str();
			}
			covStream << endl;
			for ( size_t idv = 0; idv < idvCount; ++idv ) {
				const Individual& individual = individuals[idv];
				covStream
					<< individual.getFamilyID()
					<< " "
					<< individual.getIndividualID();
				const Vector covVector = covariateMatrixTransposed.columnVector( idv );
				for ( size_t cov = 0; cov < covCount; ++cov ) {
					covStream
						<< "\t"
						<< covVector.get( cov );
				}
				covStream << endl;
			}
			covStream.close();
		} else {
			// remove disturbing covariates file, if there is any
			remove( covPath.c_str() );
			// TODO: exception on existence here
		}
		if ( 1 < traitCount ) {		// because trait[0] is treated in .fam file
			const string * traits = getTraits();
			ofstream yvmStream( yvmPath.c_str() );
			yvmStream << "FID IID";
			for ( size_t trait = 1; trait < traitCount; ++trait ) {
				yvmStream
					<< "\t"
					<< traits[ trait ].c_str();
			}
			yvmStream << endl;
			for ( size_t idv = 0; idv < idvCount; ++idv ) {
				const Individual& individual = individuals[idv];
				yvmStream
					<< individual.getFamilyID()
					<< " "
					<< individual.getIndividualID();
				const Vector traitVector = phenotypeMatrixTransposed.columnVector( idv );
				for ( size_t trait = 1; trait < traitCount; ++trait ) {
					yvmStream
						<< "\t"
						<< traitVector.get( trait );
				}
				yvmStream << endl;
			}
			yvmStream.close();
		} else {
			// remove disturbing covariates file, if there is any
			remove( yvmPath.c_str() );
			// TODO: exception on existence here
		}
	}

}
