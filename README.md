**Warning**: Repository and package names are subject to change.

**Warning**: This is heavily work in progress.


## repository contents ##
* ``bio`` reads fasta and translates genetic code
* ``matrix`` is a wrapper for some GNU Scientific Library functions (not used currently)
* ``optimize`` is the MCMC & downhill simplex implementation
* ``tree`` is tree manipulation library
* ``tlh`` is a simple tree likelihood calculation (JC69)
* ``cdlh`` is MCMC sampler for branchsite model
* ``norm`` sampler for multiple normal distribution model

## cdlh ##
* ``M0.go`` — M0 model
* ``branch_site.go`` — branch site model
* ``cdlh.go`` — main file (mcmc sampler itself)
* ``codon_frequency.go`` — F0, F3X4
* ``codon_sequences.go`` — codon alignment class
* ``ematrix.go`` — matrix class which remembers its eigen decomposition
* ``model.go`` — tree + alignment model base class

### cdlh tests ###
* ``likelihood_test.go`` — likelihood test (compare with codeml)
* ``mcmc_test.go`` — MCMC benchmark
* ``mcmcpar_test.go`` — test that likelihood is consistent during chain evaluation

## optimize ##
* ``mh.go`` — simple metropolis hastings implementation
* ``parameter.go`` — float64 parameter class
* ``adaptive.go`` — adaptive parameter class
* ``prior.go`` — prior functions
* ``proposal.go`` — proposal functions
