---
layout: page
permalink: pendulum.html
title: Pendulum exercises  
category: top
use_math: true
---

## Overview

As an example of an exercise that begins with simple physics but that can quickly get into advanced topics, we've chosen the simple pendulum. This exercise is one of my 
favorites; we're extending it a little bit given the fact that you all will be doing it with LLM help this time.

Learning Goals:

* To build a computational model of a simple pendulum using the Euler-Cromer algorithm. 
* To produce graphs of the angular displacement and angular velocity of the pendulum as a function of time from the results of the computational model. 
* To assess the accuracy of the computational results of the model.
* To identify the limitations of the small angle approximation for the pendulum. 


## The simple pendulum

### Part 1: Creating the model

Consider a simple pendulum with a point mass $m$ at the end of a string of length $L$. Performing a torque analysis to find the angular acceleration of the system, we find that $$\frac{d^2\theta}{dt^2} = -\frac{g}{L}\sin \theta$$.

This equation can be solved analytically only using the small angle approximation $\sin \theta \approx \theta$, in which case $\theta(t) = A \cos (\omega t + \phi)$, with $\omega = \sqrt{g/L}$.

First, write some code that simulates and animates a swinging pendulum. You can do this using any tool you want; the PICUP Trinket deployment is a good fallback
if you don't have access to a local programming environment.

Verify that your simulation reproduces the expected period ($\tau = 2 \pi \sqrt{\frac{L}{g}}$) at small angles.

### Part 2: Analyzing the departure from small angles qualitatively

Increase the initial amplitude from small values (e.g. 0.1 radian) to large ones (3 radian) gradually. How does the qualitative behavior of the pendulum change?

Now make a plot of the angle vs. time for different initial amplitudes. You will notice *two* separate effects as the amplitude increases -- the character of the oscillations
and their period both change. 

### Part 3: Analyzing the shift in the period quantitatively

Now make a log-log plot of the shift in period vs. initial amplitude over a wide range of initial amplitudes (starting at $10^{-5}$ radian and going up to 3 radians). 
Explain its features (some of them, but not all of them, will be physics!), and in particular, determine the slope of the long linear region. What does it tell you?
