# Pattern formation and nonlinear waves close to a 1:1 resonant Turing and Turing--Hopf instability

This repository contains the code for the numerically generated plots and the videos of pattern-forming fronts in the article "Pattern formation and nonlinear waves close to a 1:1 resonant Turing and Turing--Hopf instability" by Bastian Hilder and Christian Kuehn. The preprint can be found at [arxiv.org/abs/2508.21183](https://arxiv.org/abs/2508.21183).

### Abstract

In this paper, we analyse the dynamics of a pattern-forming system close to simultaneous Turing and Turing–Hopf instabilities, which have a 1:1 spatial resonance, that is, they have the same critical wave number. For this, we consider a system of coupled Swift–Hohenberg equations with dispersive terms and general, smooth nonlinearities. Close to the onset of instability, we derive a system of two coupled complex Ginzburg–Landau equations with a singular advection term as amplitude equations and justify the approximation by providing error estimates. We then construct space-time periodic solutions to the amplitude equations, as well as fast-travelling front solutions, which connect different space-time periodic states. This yields the existence of solutions to the pattern-forming system on a finite, but long time interval, which model the spatial transition between different patterns. The construction is based on geometric singular perturbation theory exploiting the fast travelling speed of the fronts. Finally, we construct global, spatially periodic solutions to the pattern-forming system by using centre manifold reduction, normal form theory and a variant of singular perturbation theory to handle fast oscillatory higher-order terms.

### Organisation

The repository is organised as follows:

1) The Mathematica file for the derivation of the formal amplitude equations can be found in "Derivation of amplitude equations"

2) The Matlab file to simulate the dynamics on the invariant manifold to describe the dynamics of spatially periodic solutions can be found in "Dynamics on invariant set"

3) The Matlab file to generate plots and videos of the resulting pattern-forming fronts can be found in "Pattern-forming fronts". This folder also contains the generated plots and videos.
