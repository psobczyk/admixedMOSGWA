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

#ifndef INDIVIDUAL_HPP
#define INDIVIDUAL_HPP

#include <string>
#include <ostream>

/** A biological individual with certain properties. */
class Individual {

public:

	/** Encodes biological sex or a missing value. */
	enum Sex {
		MISSING = 0,
		MALE = 1,
		FEMALE = 2
	};

private:

	/** Identifies the individual's family. */
	std::string familyId;

	/** Within the family, uniquely identifies the individual. */
	std::string individualId;

	/** Identifies the individual's father. */
	std::string paternalId;

	/** Identifies the individual's mother. */
	std::string maternalId;

	/** Specifies the individual's biological sex. */
	Sex sex;

public:

	/** Construct an Individual. */
	Individual(
		const std::string& familyId,
		const std::string& individualId,
		const std::string& paternalId,
		const std::string& maternalId,
		const Sex sex
	);

	// Getters
	std::string getFamilyID () const;
	std::string getIndividualID () const;
	std::string getPaternalID () const;
	std::string getMaternalID () const;
	Sex getSexCode () const;
};

/** Ouput a {@link Vector}. */
std::ostream& operator<< ( std::ostream& s, const Individual& individual );

#endif
