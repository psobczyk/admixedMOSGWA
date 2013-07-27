/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2011–2013, Erich Dolejsi, Bernhard Bodenstorfer.		*
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

#include "Individual.hpp"

using namespace std;

Individual::Individual (
	const string& familyId,
	const string& individualId,
	const string& paternalId,
	const string& maternalId,
	const Sex sex
) :
	familyId( familyId ),
	individualId( individualId ),
	paternalId( paternalId ),
	maternalId( maternalId ),
	sex( sex )
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
		<< ")";
	return s;
}
