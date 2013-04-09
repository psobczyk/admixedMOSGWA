#ifndef ONESNPALLIND_HPP
#define ONESNPALLIND_HPP

#include <string>

#include "Parameter.hpp"
#include "Helpfull.hpp"
#include "Log.hpp"

/** Stores a cross section of the major mode genotypes for one SNP and for all individuals */
class OneSnpAllInd {

public:
	vector<bool> one_;
	vector<bool> two_;
};
#endif
