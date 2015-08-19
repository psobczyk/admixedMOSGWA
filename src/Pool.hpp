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
#include <string>

class PoolItem
{
  friend std::ostream& operator << ( std::ostream& out, const PoolItem& p );

  /** the postition - no of the SNPs in the Model */
  std::vector<size_t> snps;

  /** model selection criterion of this model */
  long double msc;     // model seection criterion
  long double heriatability;    // Heritability of the model
  unsigned int id;
  char char_id;

public: 

  PoolItem( std::vector<size_t> &snp, long double msc, long double h, unsigned int id, char char_id );

  PoolItem( const PoolItem &p );

  PoolItem& operator= ( const PoolItem &p );

  ~PoolItem() {}

  bool operator < ( const PoolItem &p ) const;

  unsigned int getModelSize () const { return snps.size(); }

  std::vector<size_t> getPoolItem () const { return snps; }
  
  long double getMsc() const {return msc;}
  long double getHeritability() const {return heriatability;}
  
  bool findSns( size_t s );
  unsigned int getID() const {return id;}
};

#endif
