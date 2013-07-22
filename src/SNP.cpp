#include "SNP.hpp"
#include "Exception.hpp"
#include <ostream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

SNP::SNP (
	const unsigned int chromosomeId,
	const string &snpId,
	const double geneticDistance,
	const int basePairPosition,
	const char allele1,
	const char allele2
) :
	chromosomeId( chromosomeId ),
	snpId( snpId ),
	geneticDistance( geneticDistance ),
	basePairPosition( basePairPosition ),
	allele1( allele1 ),
	allele2( allele2 )
{
}

SNP::SNP (
	const string &chromosomeId,
	const string &snpId,
	const double geneticDistance,
	const int basePairPosition,
	const char allele1,
	const char allele2
) :
	snpId( snpId ),
	geneticDistance( geneticDistance ),
	basePairPosition( basePairPosition ),
	allele1( allele1 ),
	allele2( allele2 )
{
	const char * const chrStr = chromosomeId.c_str();
	const int chr = atoi( chrStr );
	if ( 0 == chr && 0 != strncmp( "0", chrStr, 2 ) ) {
		//check if chromosome is a special case
		if ( 0 == strncmp( "X", chrStr, 2 ) || 0 == strncmp( "x", chrStr, 2 ) ) {
			this->chromosomeId = 23;
		} else if ( 0 == strncmp( "Y", chrStr, 2 ) || 0 == strncmp( "y", chrStr, 2 ) ) {
			this->chromosomeId = 24;
		} else if ( 0 == strncmp( "XY", chrStr, 3 ) || 0 == strncmp( "xy", chrStr, 3 ) ) {
			this->chromosomeId = 25;
		} else if (
			0 == strncmp( "MT", chrStr, 3 ) || 0 == strncmp( "mt", chrStr, 3 )
			||
			0 == strncmp( "M", chrStr, 3 ) || 0 == strncmp( "m", chrStr, 3 )
		) {
			this->chromosomeId = 26;
		} else {
			// printLOG( "Warning: Invalid Chromosome-Name \"" + chromosome + "\" detected: using 0 for unplaced" );
			throw Exception( "Invalid chromosome-name \"%s\".", chrStr );
			this->chromosomeId = 0;
		}
	} else if ( 255 >= chr ) {
		this->chromosomeId = chr;
	} else {
		// printLOG( "Warning: Invalid Chromosome-Name \"" + chromosome + "\" detected: using 0 for unplaced" );
		throw Exception( "Invalid chromosome-name \"%s\".", chrStr );
		this->chromosomeId = 0;
	}
}

// TODO<BB>: factor out single marker test value from SNP class.
void SNP::setSingleMarkerTest ( const double singleMarkerTest ) {
	this->singleMarkerTest = singleMarkerTest;
}

// TODO<BB>: Rename to getChromosomeId.
int SNP::getChromosome () const {
	return chromosomeId;
}

std::string SNP::getSnpId () const {
	return snpId;
}

double SNP::getGeneticDistance () const {
	return geneticDistance;
}

int SNP::getBasePairPosition () const {
	return basePairPosition;
}

char SNP::getAllele1 () const {
	return allele1;
}

char SNP::getAllele2 () const {
	return allele2;
}

double SNP::getSingleMarkerTest () const {
	return singleMarkerTest;
}

ostream& operator<< ( ostream& s, const SNP& snp ) {
	return s
		<< "SNP( \""
		<< snp.getChromosome()
		<< "\",\""
		<< snp.getSnpId()
		<< "\",\""
		<< snp.getGeneticDistance()
		<< "\",\""
		<< snp.getBasePairPosition()
		<< "\",\""
		<< snp.getAllele1()
		<< "\",\""
		<< snp.getAllele2()
		<< "\",\""
		<< snp.getSingleMarkerTest()
		<< ")";
}
