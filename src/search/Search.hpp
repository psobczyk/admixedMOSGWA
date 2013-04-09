#ifndef SEARCH_SEARCH_HPP
#define SEARCH_SEARCH_HPP

#include "../io/Input.hpp"
#include "../lookup/ResultStore.hpp"

/** Search algorithms. These try to find optimal {@link Model}s for given {@link Input} data.
* @author Bernhard Bodenstorfer
*/
namespace search {

	/** Abstract base class for search algorithms. */
	class Search {

		protected:

		/** The source of data with regard to which to search. */
		io::Input& input;

		/** Already calculated model selection criteria. */
		lookup::ResultStore& resultStore;

		/** Construct with reference to data and result store. */
		Search ( io::Input& input, lookup::ResultStore& resultStore );

		public:

		/** Perform the search. */
		virtual void run () = 0;

		/** Destruct and potentially free resources. */
		virtual ~Search ();
	};

}

#endif	/* SEARCH_SEARCH_HPP */
