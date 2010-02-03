#lang scribble/manual

@(require (for-label scheme)
          (for-label "denest.ss")
          (for-label "mcmc.ss")
          (for-label scheme/flonum))

@title{Markov-Chain Monte Carlo With Density-Tree Evidence}

The @filepath{mcmc-denest.plt} package provides code for running
@link["http://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo"]{Markov-Chain
Monte Carlo} calculations in PLT Scheme.  The @filepath["denest.ss"]
module contains an algorithm of my own design based (roughly) on a
@link["http://en.wikipedia.org/wiki/Kd-tree"]{kD-tree} datastructure
that allows to approximate the Bayesian evidence from a set of samples
output from an MCMC.  (Traditionally, a major weakness of the MCMC
approach to Bayesian parameter estimation has been the difficulty of
computing the evidence from the output of a MCMC.  The evidence is
required to compute the Bayes ratio for model selection.  This
difficulty has led to approaches such as the
@link["http://www.mrao.cam.ac.uk/software/multinest/"]{MultiNest}
algorithm, or the evidence approximation techniques in
@link["http://arxiv.org/abs/0911.2150v1"]{arXiv:0911.2150}.)

The following sections document the procedures exported from the
package, and attempt to give a feeling for the algorithms involved.
The @secref{License} section describes the license under which this
code is released: the GPL, v3.

@section{Markov-Chain Monte Carlo}

@defmodule[(planet wmfarr/mcmc-denest:1:0/mcmc)]

A Markov-Chain Monte Carlo (MCMC) is a sequence of samples (a "chain")
chosen in such a way that the probability of appearance of a
particular sample point in the chain asymptotically approaches a given
probability distribution.  Often, the probability distribution that is
being approximated is a Bayesian estimate of the "posterior
probability" of some set of parameters in a model given a set of
measurements in an experiment.  In this case, examining the chain can
provide estimates of the quantities of interest.  For example, to
compute the marginalized probability distribution of a particular
parameter, it is only necessary to form a histogram of the number of
times various values of that parameter appear in the chain.

The Markov property of the chain from an MCMC means that it is
possible to generate subsequent states in the chain from only
knowledge of the current state.  Therefore, the
@scheme[make-mcmc-sampler] procedure produces a generator thunk that,
when applied to the current state of the chan, produces the next
value; in this way, chains of arbitrary length can be produced, and
processed, "on-line" without requiring arbitrary amounts of memory to
store past chain values.

The @scheme[make-mcmc-sampler] procedure uses the
@link["http://en.wikipedia.org/wiki/Metropolis-Hastings_algorithm"]{Metropolis-Hastings
algorithm} to choose the samples that go into the chain.  This
algorithm requires a probabilistic choice between the current sample
and a proposed sample; the value of the
@scheme[current-pseudo-random-generator] is used to generate the
random numbers required for this probabilistic choice.  PLT Scheme
uses the MRG32k3a pseudo-random number generation algorithm from
L'Eucyer, which is highly regarded in this field, and should be
adequate for this purpose.

@defproc[(make-mcmc-sampler (propose (-> any/c any))
                            (log-posterior (-> any/c inexact-real?))
                            (log-jump-probability (-> any/c any/c inexact-real?)))
         (-> (cons/c any/c inexact-real?) (cons/c any/c inexact-real?))]{

Produces a generator procedure that, when applied to the current
sample in the chain, returns the next sample.  A chain sample is a
@scheme[cons] of an arbitrary value representing the chain state and
the log-posterior of that value: @scheme[(cons state (log-posterior
state))].

The chain will asymptotically sample the distribution on samples
expressed by the exponential of the @scheme[log-posterior] procedure
(but note that this procedure need only produce a value that is the
log of a value @emph{proportional} to the posterior probability---it
need not be normalized).  The log of the posterior is used to avoid
underflow in the event of small posteriors.

The @scheme[propose] argument proposes a new state from the current
one (this function is sometimes called the "jump-proposal").  The
@scheme[log-jump-probability] procedure takes a source state and a
target state, and computes the log-probability---again up to an
arbitrary normalization factor---that the @scheme[propose] procedure
will propose a transition from the source to the target.  Because of
the details of the sampling, if @scheme[proposal] is symmetric---that
is, if
@schemeblock[
(= (log-jump-probability source target) 
   (log-jump-probability target source))]
then the sampling algorithm is insensitive to the return value of
@scheme[log-jump-probability], and @scheme[(lambda args 0.0)] is an
acceptable value for this procedure.

}
                      
@defproc[(in-mcmc (propose (-> any/c any))
                  (log-posterior (-> any/c inexact-real?))
                  (log-jump-probability (-> any/c any/c inexact-real?))
                  (start-state any/c))
         sequence?]{

Produces a sequence of MCMC samples (recall that a sample is
@scheme[(cons state (log-posterior state))]) starting with
@scheme[start-state].  This sequence will never terminate, so
@scheme[in-mcmc] must be used with another sequence that terminates,
for example
@schemeblock[
(in-parallel (in-range N-samples)
             (in-mcmc propose 
                      log-posterior 
                      log-jump-probability 
                      start-state))]

}

@section[#:tag "Denest"]{Density Estimator}

@defmodule[(planet wmfarr/mcmc-denest:1:0/denest)]

A k-D-tree-like datastructure allows to partition R^n adaptively
around a set of points, assigning each point a (hopefully small)
volume.  This allows computing quantities like the Bayesian evidence
from MCMC samples: evidence = <sum over points: posterior * volume>.
The adaptive nature of the k-D tree means that points in regions of
low density (i.e. low-posterior points) have large volumes (but
hopefully contribute little to the integral), while points in regions
of dense sampling (i.e. high-posterior points) are assigned small
volumes, allowing for more accurate integration.  The result has a
flavor of importance sampling.  To compute the evidence integral from
a tree constructed from MCMC samples, use
@schemeblock[
(fold
 (lambda (tree sum)
   (let ((volume (volume tree))
         (obj (object-density-tree-obj tree)))
     (let ((obj (if (pair? obj) (car obj) obj))) ; The tree can store
                                                 ; a list of objects
                                                 ; with the same
                                                 ; coordinates
       (fl+ sum (fl* vol (exp (cdr obj))))))
   0.0
   tree))]
(Recall that MCMC samples are @scheme[(cons state log-posterior)].)

@deftogether[
((defstruct empty-density-tree
   ((low flvector?)
    (high flvector?)))
 (defstruct object-density-tree
   ((low flvector?)
    (high flvector?)
    (obj any/c)))
 (defstruct cell-density-tree
   ((low flvector?)
    (high flvector?)
    (left-tree density-tree?)
    (right-tree density-tree?)))
 (defproc (density-tree? (obj any/c)) boolean?))]{

The three different types of density tree.  A density tree is either

@itemlist[
@item{Empty.}

@item{A container for a single object.}

@item{A cell that is split along one dimension into a left- and right-density-tree.}
]

Every density tree describes a region of R^n given by the points,
@scheme[x], such that @scheme[low] <= @scheme[x] < @scheme[high].
Objects whose coordinates lie in this region are contained in the tree
or one of its sub-trees.

The predicate @scheme[density-tree?] tests whether an object
represents a density tree.

}

@defproc[(volume (t density-tree?)) inexact-real?]{

Computes the volume of the region described by the tree.

}

@deftogether[
((defproc (objects->density-tree 
           (->coords (-> any/c flvector?))
           (combiner (-> (listof any/c) any))
           (objects (listof any/c)))
   density-tree?)
 (defproc (objects->density-tree/bounds
           (->coords (-> any/c flvector))
           (combiner (-> (listof any/c) any))
           (objects (listof any/c))
           (low flvector?)
           (high flvector?))
   density-tree?))]{

Constructs a density tree from a list of objects, using
@scheme[->coords] to find the coordinates of each object, and
@scheme[combiner] to combine a list of objects with the same
coordinates into a single object to place into an
@scheme[object-density-tree].  The
@scheme[objects->density-tree/bounds] procedure takes a @scheme[low]
and @scheme[high] bound on the uppermost level in the tree; in
@scheme[objects->density-tree], bounds are chosen to be the smallest
rectangular region enclosing all the objects.

}

@defproc[(fold (kons (-> object-density-tree? any/c any))
               (knil any/c)
               (tree density-tree?))
         any]{

Iterates over the @scheme[object-density-tree]s in @scheme[tree],
calling @scheme[kons] on each in turn. For the first application of
@scheme[kons], @scheme[knil] is used for the second argument; for
subsequent applications the result of the previous application is used
for the second argument.

}

@defproc[(find-object (->coords (-> any/c flvector))
                      (obj any/c)
                      (tree density-tree?))
         object-density-tree?]{

Returns the @scheme[object-density-tree] that would contain
@scheme[obj] in its bounds.  The @scheme[->coords] procedure is used
to find the coordinates of @scheme[obj].  This procedure can be used
for drawing samples from an approximate distribution described by a
sampled set of points: choose a point at random, and then draw
uniformly from within the region destribed by the
@scheme[object-density-tree] containing that point.

}

@section[#:tag "flvector-iter"]{Sequence Generators for Floating-Point Vectors}

@defmodule[(planet wmfarr/mcmc-denest:1:0/flvector-iter)]

A useful utility module for sequences of values from @scheme[flvector]s.  

@defproc[(in-flvector (flvec flvector?)) sequence?]{

Produces a sequence representing the values contained in
@scheme[flvec].  When @scheme[in-flvector] is used in the context of a
@scheme[for]-like loop, more efficient code can be generated than for
a generic sequence.

}

@section[#:tag "License"]{License}

This code is released under the
@link["http://www.gnu.org/licenses/gpl.html"]{GPL v3}.  A copy of the
GPL can be found in the @filepath{COPYING} file in the main directory
of this package.
