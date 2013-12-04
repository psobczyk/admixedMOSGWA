/********************************************************************************
 *	This file is part of the MOSGWA program code.				*
 *	Copyright ©2012–2013, Artur Gola, Erich Dolejsi, Bernhard Bodenstorfer.	*
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

#ifndef POOL_HPP
#define POOL_HPP

#include <vector>
#include <iostream>
#include "types.hpp"

class PoolItem
{
  friend std::ostream& operator << ( std::ostream& out, const PoolItem& p );

  /** the postition - no of the SNPs in the Model */
  std::vector<snp_index_t> snps;

  /** model selection criterion of this model */
  long double msc;

public: 

  PoolItem( std::vector<snp_index_t> &snp, double msc );

  PoolItem( const PoolItem &p );

  PoolItem& operator= ( const PoolItem &p );

  ~PoolItem() {}

  bool operator < ( const PoolItem &p ) const;

  unsigned int getModelSize () const { return snps.size(); }

  std::vector<snp_index_t> getPoolItem () const { return snps; }
  
  long double getMsc() const {return msc;}
  
  bool findSns(snp_index_t s);
  
};

#endif
