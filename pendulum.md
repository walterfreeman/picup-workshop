---
layout: page
permalink: pendulum.html
category: top
title: Pendulum exercises  
use_math: true
---

## Overview

As an example of an exercise that begins with simple physics (and that you all can code from scratch within this workshop) but that can quickly get into advanced topics, we've chosen the simple pendulum.

In this workshop, we imagine that you all will create a pendulum simulation using any tool you like, and then explore one of two paths: either examining the effects of damping and driving, or studying perturbation theory via the breakdown of the small-angle approximation.

Learning Goals:

* To build a computational model of a simple pendulum using the Euler-Cromer algorithm. 
* To produce graphs of the angular displacement and angular velocity of the pendulum as a function of time from the results of the computational model. 
* To assess the accuracy of the computational results of the model.
* To identify the limitations of the small angle approximation for the pendulum. 


## The simple pendulum

### Part 1: Creating the model

Consider a simple pendulum with a point mass $m$ at the end of a string of length $L$. Performing a torque analysis to find the angular acceleration of the system, we find that $$\frac{d^2\theta}{dt^2} = -\frac{g}{L}\sin \theta$$.

This equation can be solved analytically only using the small angle approximation $\sin \theta \approx \theta$, in which case $\theta(t) = A \cos (\omega t + \phi)$, with $\omega = \sqrt{g/L}$.

First, solve this numerically using the Euler-Cromer method, and generate a plot of $\theta(t)$ vs. $t$

* $\omega_{n} = \omega_{n-1} + \alpha_{n-1} \Delta t$
* $\theta_{n} = \theta_{n-1} + \omega_n \Delta t$.

Use the following parameters:

* Mass: $m = 1 {\rm kg}$
* Length: $L = 1 \rm m$
* Gravitational acceleration: $g = 9.81 \rm m/\rm s^2$

and the following initial conditions:

* Initial angular displacement: $A=\theta_0=0.1$ radian
* Initial angular velocity: $\omega_0 = 0$
* Time step: $\Delta t = 1 \rm s$

Change $\Delta t$ to smaller values until the solution converges (i.e. it does not change significantly as you make $\Delta t$ smaller). 

You can do this in any tool you like. For those who'd like to do this in a spreadsheet, there is a <a href="DampedDrivenPendulumTemplate.xlsx">spreadsheet template</a>. If you'd like to do this in Glowscript (a web-based version of Visual Python), there is
a Trinket below.

<iframe src="https://trinket.io/embed/glowscript/afb30d4fc5" width="100%" height="600" frameborder="0" marginwidth="0" marginheight="0" allowfullscreen></iframe>

From here, you can either continue with Parts 2 and 3 below to investigate the effects of the damped driven pendulum, or investigate the details of how the small-angle approximation breaks down and how this relates to <a href="perturb.html">perturbation theory</a>.

### Part 2: Adding damping

Friction and air resistance can be included as a damping term, which is proportional to the angular velocity: $\alpha_{\rm damp} = - \frac{b}{mL^2}\omega$, where $b$ is a constant that characterizes the damping effects.

1. Modify your angular acceleration to include this term. Use $b = 0.2\, {\rm kg}\, {\rm m^2/\rm s}$, then plot $\theta(t)$ vs. $t$.
2. How many oscillations will it take for the amplitude to decrease to $1/e$ of its initial value?

### Part 3: Adding driving

The pendulum can be driven by a sinusoidal torque, which will add a third term (this one positive) to the angular acceleration:

$\alpha_{\rm drive} = \frac{\tau_d}{mL^2} \cos (\omega_d t)$

1. Modify your angular acceleration to include this term. Use $\omega_d = 1$ rad/s and $\tau_d$ = 0.3 Nm, then plot $\theta(t)$ vs $t$ for at least five oscillations. 
2. Using the small angle approximation, the natural angular frequency of this system is $\omega_0=\sqrt{g/L}$. What happens to the amplitude of the oscillations when the driving frequency is...
  * $0.5\omega_0$?
  * $0.9\omega_0$?
  * $1.0\omega_0$? 
  * $1.1\omega_0$?

If you're done with this, you can investigate the <a href="perturb.html">perturbation-theory exercises</a>.
