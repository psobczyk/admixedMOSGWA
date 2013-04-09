#ifndef MOSGWA_TYPES_H
#define MOSGWA_TYPES_H

#include <stdint.h>

/** The index identifying an {@link SNP}. */
typedef uint32_t snp_index_t;

/** Maximum number of SNPs in a {@link Model}.
* Models will typically be much smaller.
*/
#define MAX_MODEL_SIZE	4294967295U

#endif	/* MOSGWA_TYPES_H */
