[![Build Status](https://drone.io/bitbucket.org/Davydov/godon/status.png)](https://drone.io/bitbucket.org/Davydov/godon/latest?branch=master)

Godon is codon models software written in go.

Godon currently supports
[M0](http://mbe.oxfordjournals.org/content/11/5/725.abstract) and
[Branch-Site](http:/mbe.oxfordjournals.org/content/22/12/2472) models
of positive selection.

Godon supports likelihood optimization
([L-BFGS-B](https://en.wikipedia.org/wiki/Limited-memory_BFGS#L-BFGS-B),
[downhill simplex](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
&
[simulated annealing](https://en.wikipedia.org/wiki/Simulated_annealing))
as well as MCMC
([Metropolis-Hastings algorithm](https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm)).

Godon supports state aggregation (options `-observed` and `-fixed`,
the last is applied only to the absolutely conserved positions). See
[the paper](http://biorxiv.org/content/early/2015/12/24/035063) for
the details.

**Warning**: Godon is currently in an early beta stage.

## Examples

Perform likelihood maximization using L-BFGS-B optimizer for the
Branch-Site model without ooptimizing the branch lengths (use only a
single CPU).
```
#!bash
$ godon -cpu 1 -nobrlen -method lbfgsb -model BS ali.fst tree.nwk
```

Run MCMC using M0 model.
```
#!bash
$ godon -method mh -model M0 ali.fst tree.nwk
```


## repository contents ##
* ``bio`` reads fasta and translates genetic code
* ``matrix`` is a wrapper for some GNU Scientific Library functions
  (not used currently)
* ``optimize`` is the MCMC & downhill simplex implementation
* ``tree`` is tree manipulation library
* ``tlh`` is a simple tree likelihood calculation (JC69)
* ``godon`` is MCMC sampler/maximum likelihood for M0 and branchsite
  model
* ``norm`` sampler for multiple normal distribution model

### cmodel ###
* ``M0.go`` — M0 model
* ``branch_site.go`` — branch site model
* ``codon_frequency.go`` — F0, F3X4
* ``codon_sequences.go`` — codon alignment class
* ``ematrix.go`` — matrix class which remembers its eigen
  decomposition
* ``matrix.go`` — codon transition matrix
* ``model.go`` — tree + alignment model base class
* ``tools.go`` — misc helper functions

#### cmodel tests ####
* ``likelihood_test.go`` — likelihood test (compare with codeml)
* ``mcmc_test.go`` — MCMC benchmark
* ``mcmcpar_test.go`` — test that likelihood is consistent during
  chain evaluation

### optimize ###
* ``mh.go`` — metropolis hastings & simulated annealing
  implementations
* ``parameter.go`` — float64 parameter class
* ``adaptive.go`` — adaptive parameter class
* ``prior.go`` — prior functions
* ``proposal.go`` — proposal functions
* ``optimizer.go`` — Optimizer and Optimizable intefaces
* ``simplex.go`` — simplex method
* ``utils.go`` — helper functions
