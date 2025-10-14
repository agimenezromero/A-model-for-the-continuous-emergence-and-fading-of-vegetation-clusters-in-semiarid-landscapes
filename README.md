# A model for the continuous emergence and fading of vegetation clusters in semiarid landscapes
**A stochastic individual-based model of vegetation dynamics in semiarid ecosystems**

This repository contains the simulation code and example notebooks developed for the study:

> *A model for the continuous emergence and fading of vegetation clusters in semiarid landscapes* (Submitted)

The model explores how the interplay between **local facilitation** (Allee effect), **competition**, and **environmental stochasticity** influences the spatial organization and persistence of vegetation in semiarid ecosystems. By explicitly representing individual plants and their spatial interactions, the simulations provide mechanistic insight into the processes driving the irregular emergence and disappearance of vegetation clusters.

---

## Overview

The code implements a **spatially explicit, stochastic individual-based model (IBM)** in which:

- Birth and death rates depend on the **local density** of individuals, mediated by **Gaussian competition and facilitation kernels**.  
- **Environmental fluctuations** are represented by stochastic alternation between low and high facilitation regimes, simulating irregular precipitation events.  
- Vegetation clustering is quantified using a **pair-correlation function**, allowing comparison with spatially uniform (mean-field) dynamics.  

---

## Repository Structure
