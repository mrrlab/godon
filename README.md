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

Godon supports state aggregation (option `--aggregate`). See
[the paper](http://biorxiv.org/content/early/2015/12/24/035063) for
the details. For the paper we used v0.5 (39bf774). Since then
likelihood computations code were substantially changed.

**Warning**: Godon is currently in a beta stage.

## Installation

The software was tested only on Linux so far.

Requirements:

* Go (preferably v1.7 or later)
* Git
* C and Fortran compilers
* [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt)
* BLAS (e.g. [OpenBLAS](http://www.openblas.net/))
* [Gonum BLAS](https://github.com/gonum/blas/) C-bindings

Once you have got all of that you can run:

```
$ bin/install.sh
```

### Full installation on Ubuntu 16.04

1. Install Go v1.7 or later. You can start by installing Go v1.6 and
   then updating using [godeb](https://github.com/niemeyer/godeb).
2. Install dependencies:
   `sudo apt-get install git libnlopt-dev libopenblas-dev build-essentials gfortran`
3. Install Gonum BLAS:
   `CGO_LDFLAGS="-lopenblas" go install github.com/gonum/blas/cgo`
4. (Optional) If your Go is older than v1.7 install
   [go-lbfsg](https://github.com/idavydov/go-lbfgsb).
5. Install godon:
   `bin/install.sh`

### Mac OS X installation

1. Make sure you have C compiler, build tools and gfortran.
2. Install Go (1.7 or later).
3. Install NLopt.
4. Get Godon source code with `go get -d bitbucket.org/Davydov/godon/godon`.
5. Install BLAS Cgo: `CGO_LDFLAGS="-framework Accelerate" go install github.com/gonum/blas/cgo`.
6. Install godon. You cannot really build a static binary for Mac OS
   X, so you need to use a `-d` option. Make sure you use correct
   paths to nlopt library and include files, as well as a correct path
   to the fortran library (`-lgfortran`).

    ```
    CGO_CFLAGS="-I/path/to/nlopt/include" CGO_LDFLAGS="-L/usr/local/Cellar/gcc/6.2.0/lib/gcc/6 -L/path/to/nlopt/lib" $GOPATH/src/bitbucket.org/Davydov/godon/bin/install.sh -d
	```

## Examples

Perform likelihood maximization using L-BFGS-B optimizer for the
Branch-Site model without optimizing the branch lengths (use only a
single CPU).
```
#!bash
$ godon -p 1 -n BS ali.fst tree.nwk
```

Run MCMC using M0 model with the downhill simplex optimization.
```
#!bash
$ godon -m mh M0 ali.fst tree.nwk
```


## repository contents ##
* ``bin`` installation script
* ``bio`` reads fasta and translates genetic code
* ``cmodel`` codon models
* ``codon`` working with codon and transition matrices
* ``godon`` is MCMC sampler/maximum likelihood for M0 and branchsite
  model
* ``misc`` various utilities
* ``optimize`` is the MCMC & downhill simplex and other algorithms implementation
* ``paml`` various functions ported from [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html)
* ``tree`` is tree manipulation library

### codon ###
* ``codon_frequency.go`` — F0, F3X4
* ``codon_sequences.go`` — codon alignment class
* ``ematrix.go`` — matrix class which remembers its eigen
  decomposition
* ``matrix.go`` — transition matrix routines

### cmodel ###
* ``aggregation.go`` — codon aggregation code
* ``branch_site.go`` — branch site model
* ``M0.go`` — M0 model
* ``model.go`` — tree + alignment model base class
* ``tools.go`` — misc helper functions

#### cmodel tests ####
* ``likelihood_test.go`` — likelihood test (compare with codeml)
* ``mcmc_test.go`` — MCMC benchmark
* ``mcmcpar_test.go`` — test that likelihood is consistent during
  chain evaluation

### optimize ###
* ``adaptive.go`` — adaptive parameter class
* ``lbfgsb.go`` — L-BFGS-B optimizer
* ``mh.go`` — metropolis hastings & simulated annealing
  implementations
* ``nlopt_callback.go`` — NLopt callback wrapper
* ``nlopt.go`` — NLopt wrapper
* ``optimizer.go`` — Optimizer and Optimizable intefaces
* ``parameter.go`` — float64 parameter class
* ``prior.go`` — prior functions
* ``proposal.go`` — proposal functions
* ``simplex.go`` — simplex method
* ``utils.go`` — helper functions

### misc ###
* ``brexp`` exports branch lengths and node labels in various
formats
* ``brmatch`` matches branch labels between two trees
* ``norm`` is a sampler for multiple normal distributions model
