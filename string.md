---
layout: page
permalink: string.html 
title: Nonlinear vibrating string 
category: top
use_math: true
---

## Overview

Here's an exercise that might be sophomore or junior-level depending on how much help students get from their LLMs. 
This is a good stress test of your system prompts -- have you calibrated them so that the LLM will help the students with the things you want,
but avoid doing things for them that you want them to learn on their own?

Note that you will probably not have time to write the code for this project yourself in an hour and a half unless you are a very experienced programmer.
The robots can do it quickly, though!

### Part 1: Modeling a vibrating string

We can model a stretched, vibrating guitar string as a set of $N$ Hooke's law strings, each with an equilibrium length $r_0$ and a spring constant $k$, which connect $N+1$
point masses each of mass $m$ which are stretched over a total distance $L$, with the two endmost point masses fixed in place. (Note that $L$ must be somewhat bigger than $Nr_0$.)

Write computer code to simulate and animate the motion of such a system.

### Part 2: Analyzing the string's behavior at low amplitude

At low amplitude, an ideal vibrating string is linear with independent normal modes which do not couple. The $n$'th transverse normal mode has period

$\tau_n = \frac{2L}{nv} = \frac{2L}{n} \sqrt{\frac{\mu}{T}}$

where $\mu$ is the linear mass density of the string and $T$ is the tension it is subjected to.

Verify for several normal modes that your string's oscillations have the expected period.

### Part 3: Analyzing the shift in period as the amplitude increases

Now increase the amplitude of the vibrations of your string. Measure the shift in period for different normal modes for a range of amplitudes. Make a log-log plot of the 
shift in period vs. amplitude. Comment on the features of your graph. Part of your graph will be linear, showing a power-law behavior; what is the power?

### Part 4: Coupling between modes

Write code that monitors the motion of one of the points in your string and computes its Fourier transform. Excite your string in the fundamental $n=1$ mode 
at low amplitude and verify that the Fourier peak is where you expect it to be.

Then increase the amplitude so that you are no longer in the low-amplitude regime. Verify that this peak has shifted consistent with the results you got in Part 3. But you will also notice *new* peaks showing up. Where 
did they come from?

Some of these peaks will be near-integer multiples of the frequency of the fundamental. Where did they come from, and why are they not *exactly* integer multiples of the 
fundamental?

However, you will notice other peaks that are *not* integer multiples of the frequency of the fundamental. What's going on there? (Hint: Are there modes of vibration
we haven't discussed yet?)


