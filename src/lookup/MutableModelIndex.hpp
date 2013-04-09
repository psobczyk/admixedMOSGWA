#ifndef MUTABLE_MODEL_INDEX_HPP
#define MUTABLE_MODEL_INDEX_HPP

#include "ModelIndex.hpp"

namespace lookup {

	/** A mutable {@link ModelIndex}. */
	class MutableModelIndex : public ModelIndex {

		public:

		/** Copy constructor, to avoid the default. */
		MutableModelIndex ( const MutableModelIndex& original );

		/** Copy from {@link ModelIndex} constructor. */
		MutableModelIndex ( const ModelIndex& original );

		/** Assignment operator, to avoid the default. */
		MutableModelIndex& operator= ( const MutableModelIndex& original );

		/** Assignment from {@link ModelIndex} operator. */
		MutableModelIndex& operator= ( const ModelIndex& original );
	};

}

#endif	/* MUTABLE_MODEL_INDEX_HPP */
