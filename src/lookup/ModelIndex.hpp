#ifndef MODEL_INDEX_HPP
#define MODEL_INDEX_HPP

#include <set>
#include <vector>
#include <iostream>
#include "../types.hpp"

namespace lookup {

	/** Represents a sparse selection of {@link SNP}s.
	* Basically, it is a subset of position-indices of SNPs of MData.
	* The objects of class <code>ModelIndex</code> are immutable
	* so that they can be used as safe map keys.
	* However, subclasses such as {@link MutableModelIndex} are not immutable.
	* The memory footpring is rather minimalistic (assuming a sparse selection)
	* so that many objects can be stored in memory.
	* Note: An alternative implementation of ModelIndex optimised for larger model sizes
	* could use a Patricia tree (radix tree) or trie to centrally administrate the multi-indices in use.
	*/
	class ModelIndex {

		/** Shared storage location for the SNP positions.
		* The private class <code>ReferenceCountedShared</code> is immutable by design.
		* It exploits immutability to easily avoid expensive copying.
		*/
		class ReferenceCountedShared {

			/** Forbid empty constructor. */
			ReferenceCountedShared ();

			/** Forbid copy constructor. */
			ReferenceCountedShared ( const ReferenceCountedShared& original );

			/** Forbid assignment. */
			ReferenceCountedShared& operator= ( const ReferenceCountedShared& original );

			public:

			/** Reference counter, i.e. number of ModelIndex objects refering to this struct */
			size_t referenceCount;

			/** Size of the model, i.e. number of SNPs */
			const size_t size;

			/** Ascendingly sorted list of SNP positions.
			* The array is of size {@link size}.
			*/
			const snp_index_t * snps;

			/** Construct an instance.
			* The obvious assumption is that upon construction, the reference count will be 1.
			*/
			ReferenceCountedShared ( const size_t size, const snp_index_t *snps );

			/** Destruct the instance, when nobody refers to it any more. */
			~ReferenceCountedShared ();
		} * snpStruct;

		/** Intialise with the indices from a set.
		* WARNING: This method is meant to be called as part of a constructor.
		* It must not be called on an already constructed object!
		* Otherwise a memory leak will occur.
		* This factoring out of common constructor code is necessary
		* because C++ constructors cannot yet delegate to each other.
		*/
		void init ( const std::set<snp_index_t> &snpSet );

		protected:

		/** Protected assignment operator.
		* Objects of this class should be immutable to make them safe map keys.
		* However, mutable subclasses may allow assignment and, hence, use this operator.
		*/
		ModelIndex& operator= ( const ModelIndex& original );

		/** Make this and that share the same internal index array.
		* The assumption is that this and that must be equal
		* in the sense of {@link operator==( const ModelIndex& that )}.
		* Otherwise, the promised immutability is broken.
		* In this case, the method is semantically <code>const</code>.
		*/
		void shareInternals ( ModelIndex& that );

		public:

		/** Iterator to retrieve the indices present in a multi-index of class ModelIndex. */
		typedef const snp_index_t* const_iterator;

		/** Construct an empty model index. */
		ModelIndex ();

		/** Construct with the indices from a set. */
		ModelIndex ( const std::set<snp_index_t> &snpSet );

		/** Convenience constructor to get indices from a vector,
		* sort them and make them unique.
		* This is meant to facilitate the move
		* from legacy vectors of SNP-indices
		* to modern ModelIndex.
		*/
		ModelIndex ( const std::vector<snp_index_t> &snpVec );

		/** Copy constructor. */
		ModelIndex ( const ModelIndex& original );

		/** Copy except for one SNP constructor.
		* This reminds of the prototyping pattern:
		* All SNPs are copied, except the given one,
		* which is added, if it is not member of the original and removed if it is so.
		* @param original to be taken as prototype to be copied except for one SNP.
		* @param snp to be added or left out when copying,
		* depending on whether it is not or is contained in the original.
		*/
		ModelIndex ( const ModelIndex& original, const snp_index_t snp );

		/** Destruct and potentially free (depending on reference count) referenced data. */
		~ModelIndex ();

		/** Compare with another.
		* This method <em>semantically</em> deals with const objects;
		* behind the scenes, however, if equality is detected,
		* both objects are made share a single reference-counted internal object.
		* @see #shareInternals( ModelIndex& that )
		* @returns -1 if this is less than that, 0 if it is equal and +1 if it is greater.
		*/
		int compare ( const ModelIndex& that ) const;

		/** Comparison for equality.
		* @see #compare( const ModelIndex& that )
		*/
		bool operator== ( const ModelIndex& that ) const;

		/** Comparison for unequality.
		* @see #compare( const ModelIndex& that )
		*/
		bool operator!= ( const ModelIndex& that ) const;

		/** An less relation consistent with equality and the other order-related operators.
		* @see #compare( const ModelIndex& that )
		*/
		bool operator< ( const ModelIndex& that ) const;

		/** An less-or-equal relation consistent with equality and the other order-related operators.
		* @see #compare( const ModelIndex& that )
		*/
		bool operator<= ( const ModelIndex& that ) const;

		/** An greater relation consistent with equality and the other order-related operators.
		* @see #compare( const ModelIndex& that )
		*/
		bool operator> ( const ModelIndex& that ) const;

		/** An greater-or-equal relation consistent with equality and the other order-related operators.
		* @see #compare( const ModelIndex& that )
		*/
		bool operator>= ( const ModelIndex& that ) const;

		/** Whether this is a subset of that.
		* This means this contains all indices of that, and potentially more.
		*/
		bool isSubsetOf ( const ModelIndex& that ) const;

		/** Whether this is a superset of that.
		* This means all indices of this are also contained in that, and potentially more.
		*/
		bool isSupersetOf ( const ModelIndex& that ) const;

		/** Get the number of SNPs. */
		size_t size () const;

		/** Quasi iterate the SNP positions in ascending order. */
		const_iterator begin () const;

		/** The last iterated position. */
		const_iterator end () const;

		/** Tests whether this ModelIndex contains a SNP, represented by its position. */
		bool contains ( const snp_index_t snp ) const;
	};

	/** Ouput a {@link ModelIndex}. */
	std::ostream& operator<< ( std::ostream& s, const ModelIndex& m );

}

#endif	/* MODEL_INDEX_HPP */
