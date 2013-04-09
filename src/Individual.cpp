#include "Individual.hpp"

using namespace std;

Individual::Individual (
	const string& familyId,
	const string& individualId,
	const string& paternalId,
	const string& maternalId,
	const Sex sex,
	const double phenotype
) :
	familyId( familyId ),
	individualId( individualId ),
	paternalId( paternalId ),
	maternalId( maternalId ),
	sex( sex ),
	phenotype( phenotype )
{
}

// TODO<BB>: separate out phenotype from individual.
// TODO<BB>: return const string& from all string getters.
// TODO<BB>: write Id instead ID and leave away the _.

string Individual::getFamilyID () const {
	return familyId;
}

string Individual::getIndividualID () const {
	return individualId;
}

string Individual::getPaternalID () const {
	return paternalId;
}

string Individual::getMaternalID () const {
	return maternalId;
}

Individual::Sex Individual::getSexCode () const {
	return sex;
}

double Individual::getPhenotype () const {
	return phenotype;
}

void Individual::setPhenotype ( const double phenotype ) {
	this->phenotype = phenotype;
}

ostream& operator<< ( ostream& s, const Individual& individual ) {
	s
		<< "Individual( \""
		<< individual.getFamilyID()
		<< "\",\""
		<< individual.getIndividualID()
		<< "\",\""
		<< individual.getPaternalID()
		<< "\",\""
		<< individual.getMaternalID()
		<< "\",";
	switch ( individual.getSexCode() ) {
		case Individual::MALE: s << "MALE"; break;
		case Individual::FEMALE: s << "FEMALE"; break;
		default: s << "MISSING";
	}
	s
		<< ","
		<< individual.getPhenotype()
		<< ")";
	return s;
}
