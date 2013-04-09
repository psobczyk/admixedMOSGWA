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
