**Warning**: Repository and package names are subject to change.

**Warning**: This is heavily work in progress.


## repository contents ##

* ``bio`` reads fasta and translates genetic code
* ``matrix`` is a wrapper for some GNU Scientific Library functions (not used currently)
* ``tlh`` is a simple tree likelihood calculation (JC69)
* ``cdlh`` is MCMC sampler for branchsite model

## cdlh content ##

* ``M0.go`` — M0 model
* ``branch_site.go`` — branch site model
* ``cdlh.go`` — main file (mcmc sampler itself)
* ``codon_frequency.go`` — F0, F3X4
* ``codon_sequences.go`` — codon alignment class
* ``ematrix.go`` — matrix class which remembers its eigen decomposition
* ``likelihood_test.go`` — likelihood test (compare with codeml)
* ``mcmc.go`` — MCMC implementation
* ``mcmc_test.go`` — MCMC benchmark
* ``model.go`` — tree + alignment model base class
* ``tools.go`` — misc tools
* ``tools_test.go`` — test for tools