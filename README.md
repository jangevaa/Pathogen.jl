# Pathogen

[![Build Status](https://travis-ci.org/jangevaa/Pathogen.jl.svg?branch=master)](https://travis-ci.org/jangevaa/Pathogen.jl)

## Introduction

Pathogen.jl is a package that provides utilities for the simulation and
inference of pathogen phylodynamics, built in the [Julia
language](http://julialang.org). Specifically, Pathogen.jl presents an extension
to the individual level infectious diesase transmission models (ILMs) of Deardon
et al. (2010), to simultaneously model infectious disease transmission and
evolution. Pathogen genomic sequences are used in conjunction with the covariate
and disease state information of individuals to infer disease transmission
pathways, external disease pressure, and infectivity parameters.

At it's core, an ILM describes the probability an individual will become exposed
to a disease in a specific time interval based on their covariate information
and the disease state of the rest of the population. In a discrete time ILM, the
probability an individual $i$ will become exposed at time $t$ is as (Deardon et
al. 2010):

$P(i,t) = \begin{cases} 1 - \exp \left \{ \left[-\Omega_{S}(i)\sum_{j \in
I_{(t)}} \Omega_{T}(j)\kappa(i,j)\right] + \epsilon(i,t) \right \} \text{ if } i
\in S_{(t)}\\
0 \text{ otherwise}
\end{cases}$,

where:

* $I_{(t)}$ is the set of infectious individuals at time $t$
* $S_{(t)}$ is the set of susceptible individuals at time $t$
* $\Omega_{S}(i)$ represents risk factors associated with susceptible individual
$i$
* $\Omega_{T}(i)$ represents risk factors associated with infectious individual
$j$
* $\kappa(i,j)$ is an infection kernel representing risk factors involving both
a susceptible individual $i$, and an infectious individual $j$
* $\epsilon(i,t)$ represent external risk factors associated with susceptible
individual $i$

In an ILM, a susceptible individual's exposure risk is the sum of various risk
factors. The specific source of an infection in this framework, however, is
irrelevant. In the phylodynamic ILM extension, the source of an infection is of
relevance to the both the pathways of pathogen evolution and of pathogen
transmission. As such, it is required that seperate exposure probabilities are
described for each susceptible-infectious individual combination. Furthermore,
as an individual may only be infected from a single source, we must determine
which of the competiting sources (if any) are responsible for their exposure to
the pathogen. To do this we are required to model in continuous time, resulting
in the following exposure probability:

$P(i, j, t) = \begin{cases} \Omega_{S}(i)\Omega_{T}(j)\kappa(i,j) \exp \left\{
-\Omega_{S}(i)\Omega_{T}(j)\kappa(i,j)\delta(t)\right\} \text{ if } i \in
S_{(t)}, j \in  I_{(t)} \\
0 \text{ otherwise}
\end{cases}
$,

where $\delta(t)$ is the time from the previous event.

The risk of exposure from external factors to an susceptible individual can also
be defined seperately as:

$P(i, t) = \begin{cases} \epsilon(i,t) \exp \left \{ \epsilon(i,t) \delta(t)
\right \} \text{ if } i \in S_{(t)} \\
0 \text{ otherwise}
\end{cases}
$.

An individual's transition from an exposed to an infectious state is assumed to
be exponentially distributed with probability:

$P(i, t) = \begin{cases} \rho(i) \exp \left \{ \rho(i) \delta(t) \right \}
\text{ if } i \in E_{(t)} \\
0 \text{ otherwise}
\end{cases}
$.

where,

* $E_{(t)}$ is the set of exposed individuals at time $t$,
* $\rho(i)$ is the infectivity rate (inverse of the mean latent period)
associated with an exposed individual $i$.

Similarly, an individual's transition from an infectious state to a removed
state is assumed to be exponentially distributed with probability:

$P(i, t) = \begin{cases} \gamma(i) \exp \left \{ \gamma(i) \delta(t) \right \}
\text{ if } i \in I_{(t)} \\
0 \text{ otherwise}
\end{cases}
$.

where,

* $\gamma(i)$ is the removal rate (inverse of the mean infectious period)
associated with an infectious individual $i$.


## Simulation framework


## Inference framework


## Future work


## References



