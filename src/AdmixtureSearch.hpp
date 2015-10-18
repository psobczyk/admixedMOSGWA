

#ifndef ADMIXTURESEARCH_HPP
#define ADMIXTURESEARCH_HPP

#include "AMdata.hpp"
#include "AModel.hpp"
#include "search/Search.hpp"

/** For a given AMdata object, AModel object and MarginalTests object
AdmixtureSearch performes greedy search over the space of all possible 
models.

*/
class AdmixtureSearch  : public search::Search {
private:
  /** Data needed by the model below. */
  AMdata data;
  
  /** The best model yet found. */
  AModel model;

  //* Marginal tests */
  MarginalTests tests;

public:
  /** Set up the search environment. */
  AdmixtureSearch ();
  
  /** Run the search. */
  virtual void run ();
  
  /** Retrieve the winning model.
   * This should be called only after {@link run}.
   */
  virtual const AModel* result ();
  

};

#endif

