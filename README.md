**Warning**: Repository and package names are subject to change.

**Warning**: This is heavily work in progress.


## repository contents ##
* ``bio`` reads fasta and translates genetic code
* ``matrix`` is a wrapper for some GNU Scientific Library functions (not used currently)
* ``optimize`` is the MCMC & downhill simplex implementation
* ``tree`` is tree manipulation library
* ``godon`` is MCMC sampler/maximum likelihood for M0 and branchsite model
* ``norm`` sampler for multiple normal distribution model

## cmodel ##
* ``M0.go`` — M0 model
* ``branch_site.go`` — branch site model
* ``codon_frequency.go`` — F0, F3X4
* ``codon_sequences.go`` — codon alignment class
* ``ematrix.go`` — matrix class which remembers its eigen decomposition
* ``matrix.go`` — codon transition matrix
* ``model.go`` — tree + alignment model base class
* ``tools.go`` — misc helper functions

### cmodel tests ###
* ``likelihood_test.go`` — likelihood test (compare with codeml)
* ``mcmc_test.go`` — MCMC benchmark
* ``mcmcpar_test.go`` — test that likelihood is consistent during chain evaluation

## optimize ##
* ``mh.go`` — simple metropolis hastings implementation
* ``parameter.go`` — float64 parameter class
* ``adaptive.go`` — adaptive parameter class
* ``prior.go`` — prior functions
* ``proposal.go`` — proposal functions
* ``optimizer.go`` — Optimizer and Optimizable intefaces
* ``simplex.go`` — simplex method
* ``utils.go`` — helper functions
