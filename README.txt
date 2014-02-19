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
cochran_armitage = true 
[model_selection]
nSNPKriterium = 5000
expected_causal_snps_MBIC = 25
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
• http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#bed
• http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

The [data] section

trait_index = 0

specifies that the phenotype for the first trait should be taken. In plink format, it is contained in the file with suffix .fam.

The [output] section specifies where log- and other output files will be written with the option

files = "random_out"

This string will be used as prefix for the output filenames.

For fine-tuning model selection, the section [model_selection] contains options

expected_causal_snps_MBIC	lets the first step with the relaxed selection criterium consider models of up to about the given size.

Further useful options:

[single_marker]
chi_square		boolean		whether to use χ²-test in the single marker phase
cochran_armitage	boolean		whether to use Cochran-Armitage test in the single marker phase

[model_selection]
maximalModelSize	integer		limits the search to models of size up to the given; saves time
PValueBorder		integer		only so many SNPs are considered in multi-forward steps, ranked by p-value
forward_step_max	integer		bounds the numer of SNPs in the forward step from the empty model
fast_multi_forward	integer		bounds the numer of SNPs in the forward step from nonempty models
nSNPKriterium		integer		useful for running with a subset of top-ranking SNPs: the original number of SNPs, to be used by mBIC 
expected_causal_snps_MBIC	integer	indication to the initial search with mBIC (not mBIC2!) how many SNPs to expect. The decisive second run with mBIC2 does not need this parameter.

Upon successful run, you will find (assuming output filename prefix "random_out") files with the names

random_out_IT.txt			results from individual SNP tests
random_outYvecout			states the phenotype vector used
random_out.mod				describes the chosen model
random_out.log				log file from the search
random_out_0the_result_Corr.txt		information about SNPs which are highly correlated to those in the model
random_out0the_resultCorr.h5		similar, but in HDF5 format
