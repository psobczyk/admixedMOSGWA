#ifndef SORT_VEC_HPP
#define SORT_VEC_HPP

#include <vector>
#include <algorithm>

using namespace std;

/** A sorted list of (int id, double value) pairs.
* Access by index, not by iterators.
*/
class SortVec {

	private:

	/** One pair to sort */
	struct SortItem {

		/** Used for SNP id */
		const int id;

		/** Used for P-value */
		const double value;

		/** Construct a pair */
		SortItem ( const int id, const double value );
	};

	friend bool order_function ( const SortVec::SortItem* i, const SortVec::SortItem* j );

friend bool order_function2 ( const SortVec::SortItem* i, const SortVec::SortItem* j );
	/** Vector of sortable pairs */
	vector <SortItem*> list_;

	/** Forbid use of the (default) copy constructor. */
	SortVec ( const SortVec& original );

	/** Forbid use of the (default) assignment operator. */
	SortVec& operator= ( const SortVec& original );

	protected:

	/** Clear all entries. */
	void clear ();

	public:

	/** Default constructor */
	SortVec ();

	/** Constructor with a given size */
	SortVec ( const int n );

	/** Constructor with arrays containing position and values, so that  ids[i], values[i] are a pair */
	SortVec ( const int n, const int ids[], const double values[],bool bigger=true );

	/** Destructor */
	~SortVec ();

	/** Clear and overwrite SortVec with arrays */
	void fillVec ( const int n, const int ids[], const double values[], bool bigger=true);

	/** Get position (or id_) for the k-th smallest value */
	int getId ( const int k ) const;

	/** Get the k-th smallest value */
	double getValue ( const int k ) const;
};

#endif	/* SORT_VEC_HPP */
