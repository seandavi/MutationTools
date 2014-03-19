## MutationTools

This is a small collection of tools for working with variants derived from DNA
sequencing, with a particular focus on cancer mutation.

## Installation

There are a couple of choices.  Installation directly from github is possible using the `devtools` package:

    library(devtools)
	install_github('MutationTools','seandavi')

Or, install from my software repository using `biocLite`:

    library(BiocInstaller)
	biocLite('MutationTools',siteRepos='http://watson.nci.nih.gov/~sdavis/software/R',type='source')
