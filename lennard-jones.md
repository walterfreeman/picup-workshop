---
layout: page
permalink: lj.html
category: top
title: Lennard-Jones gas 
use_math: true
---

## The Lennard-Jones gas

This is a "graduate-level" exercise set: mapping the phase diagram of the truncated Lennard-Jones gas. It will likely require substantial LLM assistance to code 
in an hour and a half, and it is probably a good demo of the physics you can learn without actually writing the code yourself. Choose this option if you are interested
in exploring what students can do with substantial LLM assistance.

## Introduction

The Lennard-Jones force can be defined as the superposition of a short-ranged attractive force and a very-short-ranged repulsive force:

$F_{LJ} = {\left(\frac{r}{r_0}\right)^{\alpha}} - {\left(\frac{r}{r_0}\right)^{\beta}}$

where $\alpha$ is -13 (very short ranged repulsion) and $\beta$ is -7 (longer ranged attraction); $r_0$ is the distance at which the force is zero (the equilibrium radius). This is used as a model of atoms that interact only via the van der Waals force.

Switching to nondimensional units where $r_0 = 1$, this gives a potential

$V = -\int \vec F \cdot d\vec r = \frac{1}{12 r^{12}} - \frac{1}{6r^6}$

### The truncated L-J force

This force can be truncated to be exactly zero beyond some cutoff distance $r_{\rm max}$, where $r_{\rm max}$ is several times $r_0$ so that the force is already small by that point. This allows for significant computational efficiency gains since long-range forces go from "small" to "exactly zero".

To keep the potential self-consistent, we must shift it upward at all radii by an amount equal to

$\frac{1}{6r_{\rm max}^6} - \frac{1}{12 r_{\rm max}^{12}}$

so that it remains zero at infinity.

## The Exercises

*(These are intentionally left abbreviated since you are likely doing this with LLM assistance! You can do this in any way that you and your robot helper think is appropriate to your computing environment.)*

1. Write some code that simulates the truncated Lennard-Jones force using the leapfrog or Euler-Cromer integrator for a set of $N$ particles in two or three dimensions confined to a box. You may use any programming environment you have on your computer or trinket.gopicup.org, a high-performance web-based Python interpreter.

   This code should also visualize what your particles are doing in realtime.

2. How many particles can your code simulate with reasonable performance?

3. You will see that the equilibrium behavior of your simulation depends only on the initial total energy, since it is a microcanonical ensemble coming to equilibrium. Can you make your model reproduce gas-like behavior? Liquid-like behavior?

4. Confirm that your simulation reproduces the ideal gas law in the limit where it is hot and diffuse.

5. Map out the phase transition as a function of temperature: plot the equilibrium value of $\frac{NkT}{PV}$ and see if the transition is apparent. Does the sharpness of the transition depend on the number of particles in your simulation?

6. Now see if you can generate a *solid* phase for your simulation in 3D. (You will either need to be very clever about your initial conditions or add a dissipative force to further bleed off energy.) Determine the freezing temperature kT and the Bravais lattice type of the solid phase.
