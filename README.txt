Welcome to MOSGWA

MOSGWA is a tool meant to be operated from the Linux command line.
It is planned to make it compile on other operating systems.

In the top directory, you see:
README.txt	This overview information
INSTALL.txt	Build instructions
CHANGES.txt	Change log
COPYING.txt	GNU General Public License, which applies to this software
CMakeLists.txt	Top level configuration file for build with the cmake tool
src		Contains the C++ source and header files and a suitable makefile

Installation:
Follow the steps described in INSTALL.txt.

Running:
MOSGWA takes its configuration from files given on the command line.
You run MOSGWA with the command syntax:

MOSGWA config_file_name[s]

Config files look similar to Windows INI-files.
They determine files used, and any parameters for the search strategy,
in cases when the default values are not deemed optimal.
The following is an example:

[input]
plink_files = "random"
[data]
trait_index = 0
[output]
files = "random_out"
[single_marker]
test = cochran_armitage
[model_preselection]
mBIC_expected_causal_SNPs = 25
[model_selection]
selection_criterium = mBIC2
regression_type = firth
fast_multi_forward = false

You see the sections of the file headed by section headings, which are enclosed in square brackets [].
Within each section, the names of the parameters are unique.
You set parameters with an equals sign.

MOSGWA currently uses four types of parameters.
boolean (true or false)
integer (e.g. 0, 1, 2)
floating point (e.g. -9.3e3)
string (e.g. "random_out")

The [input] section must specify where to read the data from.

plink_files = "random"

specifies that the input format is plink's binary format, and the files to read are:
random.bim			contains information about SNPs
random.fam			contains information about individuals including the phenotype for one trait
random.bed			contains the fact table of genotypes
random.cov	if existing	contains additional covariates if there are any
random.yvm	if existing	contains phenotypes for additional traits if there are any

Concerning the file formats see:
* http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
* http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

The [data] section

trait_index = 0

specifies that the phenotype for the first trait should be taken. In plink format, it is contained in the file with suffix .fam.

The [output] section specifies where log- and other output files will be written with the option

files = "random_out"

This string will be used as prefix for the output filenames.

For fine-tuning model selection, the section [model_selection] contains options

expected_causal_snps_MBIC	lets the first step with the relaxed selection criterium consider models of up to about the given size.

Further useful options:

[log]
level			choice		DEBUG, INFO (default), WARNING, ERROR: how much to log

[input]
cache_limit		integer		genotype data for how many SNPs may be cached in memory, default 1024

[output]
correlation_threshold	double		restricts the listing of closely correlated SNPs, default 0.999

[single_marker]
test			choice		chi_square, cochran_armitage ... which one to use

[model_preselection]
mBIC_expected_causal_SNPs	integer		parameter for first round of model search (with mBIC), which determines the starting point for the actual model search

[model_selection]
regression_type		choice		linear (default), firth: to use for calculating the model selection criterium
selection_criterium	choice		mBIC or mBIC2 (default): criterium to use in second round of model selection
mBIC_expected_causal_SNPs	integer		parameter for second round of model search, when mBIC is chosen, irrelevant for mBIC2
EBIC_gamma		double		parameter in the EBIC criterium, defaults to 1 - log( #individuals ) / ( 2 * log( #SNPs ) )
search_strategy		choice		greedy (default), memetic_algorithm
maximalModelSize	integer		limits the search to models of size up to the given; saves time
PValueBorder		size_type	only so many SNPs are considered in multi-forward steps, ranked by p-value
forward_step_max	integer		bounds the numer of SNPs in the forward step from the empty model
fast_multi_forward	integer		bounds the numer of SNPs in the forward step from nonempty models
nSNPKriterium		integer		useful for running with a subset of top-ranking SNPs: the original number #SNPs of SNPs, to be used by the selection criteria; not used if 0, which is the default

Upon successful run, you will find (assuming output filename prefix "random_out") files with the names

random_out_IT.txt			results from individual SNP tests
random_outYvecout			states the phenotype vector used
random_out.mod				describes the chosen model
random_out.log				log file from the search
random_out_0the_result_Corr.txt		information about SNPs which are highly correlated to those in the model
random_out0the_resultCorr.h5		similar, but in HDF5 format
