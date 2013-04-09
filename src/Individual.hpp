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

	/** Specifies the individual's affectedness. */
	double phenotype;

public:

	/** Construct an Individual. */
	Individual(
		const std::string& familyId,
		const std::string& individualId,
		const std::string& paternalId,
		const std::string& maternalId,
		const Sex sex,
		const double phenotype
	);

	// Getters
	std::string getFamilyID () const;
	std::string getIndividualID () const;
	std::string getPaternalID () const;
	std::string getMaternalID () const;
	Sex getSexCode () const;
	double getPhenotype () const;

	/** Modify the stored phenotype value. */
	void setPhenotype ( const double phenotype );
};

/** Ouput a {@link Vector}. */
std::ostream& operator<< ( std::ostream& s, const Individual& individual );

#endif
