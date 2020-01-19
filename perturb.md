---
layout: page
permalink: perturb.html
category: top
title: Perturbation Theory and Pendula
use_math: true
---

## Overview 

Using computational tools doesn't just let us do a better job teaching students the same stuff we were already going to teach them. It gives us access to new ways to teach them *new things* that would otherwise be impossible.

(Snowshoes vs. boots)

In this workshop segment, I want to explore an example of how you can add to your learning goals to look for opportunities to encourage expert-level thinking among students that you couldn't otherwise do.

## An example: the swinging pendulum

The standard presentation of the pendulum goes something like this: "Here's the ODE for a swinging pendulum. This has no analytic solution if $\theta$ is large so let's assume it's small.
Then it reduces to a simple harmonic oscillator."

The previous presentation discussed how we can teach students about simple harmonic oscillators -- a crucial thing for students to understand in order to develop expert-level understanding of physics -- more effectively by using 
numerical methods. 

But ... if we're using numerical methods, we don't have to make the small-$\theta$ assumption! And, if we don't, what *else* could we teach them that could help them think about physics 
in expert ways sooner?

### Perturbation theory and perturbative reasoning

How do we think about the pendulum? It's a perturbative problem. We (physicists) understand that the SHO approximation is just the solution to the unperturbed problem. Perturbation theory is essential for students to understand. So why not teach them the beginnings of perturbation theory?

The usual answer: "we don't because the math is too hard". Our answer: "we have a computer to do the math for us!"

The basics don't involve any complicated analytical math:

* Simple problems often have simple solutions (SHO)
* Problems that are close to simple problems have solutions that are close to those simple solutions
* The mathematics of "close to" is a power series: if a problem is "simple plus small adjustment", the solution is "simple plus power series in size of adjustment"

The following is an exercise that students in an introductory class can do, written in a style appropriate to present to students. I've made the following choices here:

* Audience: first-year physics majors or talented physics students who would benefit from understanding perturbative reasoning
* Tools: Glowscript Python via Trinket
* Format: instructor provides minimally-working code skeleton; students fill in missing physics and computational-science elements without getting into the minutiae of Glowscript

Many other formats are feasible, too. In fact, in my class last semester, I made the following choices to fit my audience:

* Audience: students in a computational physics course 
* Tools: C on Linux, using anim/gnuplot for visualization
* Format: students write all code from scratch

In this workshop, feel free to write your own code with whatever tools you like, or fiddle with my Trinkets below. If you've done the previous in a spreadsheet, some of the later steps may require more manual work than you'd like, so
you may want to switch to the Trinkets below.

## The exercise

### Part 1: qualitative study of the swinging pendulum at large amplitude

If you've not gotten the pendulum simulation implemented in your chosen tool already, do that first. If you're using my Trinkets, you'll need to:

* Think about what values for g and L you want to use, and fill them in
* Think about a sensible value for the timestep dt, and add it. (Remember computers are fast!)
* Set some initial conditions for your pendulum (and the stopwatch you are using to time it)
* Fill in the physics update: how do theta, omega, and t change during each timestep?

Then, whatever tool you're using, investigate what happens as the amplitude becomes large:

* Play with $\theta_0$ and see how the behavior changes. What happens for an amplitude of 3 radians? Is this physical?
* Look at the shapes of the $\theta$ vs. $t$ graphs for amplitudes of 3 radians, 1.5 radians, and in the small-amplitude limit. How do they differ? Is this consistent with the SHO approximation for the pendulum?
* What happens in the limit where the amplitude approaches $\pi$? Is this consistent with the SHO approximation?
      
<iframe src="https://trinket.io/embed/glowscript/afb30d4fc5" width="100%" height="600" frameborder="0" marginwidth="0" marginheight="0" allowfullscreen></iframe>

### Part 2: quantitatively measuring the period

 Modify the condition on the while loop: how can you make the pendulum stop after one-half swing?
 The code in the above Trinket will then print out the period and end. How does the period depend on the amplitude?


### Part 3a (if you have plenty of time and want to hack yourself): period vs. amplitude

* Fill in the missing pieces of the below code to run many simulations, one after the other, and compute the "period anomaly" for each
* Investigate the resulting plot of period anomaly vs. amplitude. How does this relate to broader physics ideas? In particular:
  * Why does the plot look "stair-stepped" at very low amplitude?
  * Why is the anomaly quadratic in the angle for intermediate amplitude?
  * Why does the anomaly increase faster than angle squared for very large angle?

  <iframe src="https://trinket.io/embed/glowscript/6e12850f82" width="100%" height="600" frameborder="0" marginwidth="0" marginheight="0" allowfullscreen></iframe>

### Part 3b (if you are short on time and want to skip to the analysis, for the workshop): period vs. amplitude

* Run the below code. What is it doing?
* The "period anomaly" is the difference between the simulated period and the one predicted by the small-angle approximation formula in the textbooks. The code will plot it in real time. You may wish to make a log-log plot to better analyze your data. Explain the following features:
  * Why does the plot look "stair-stepped" at very low amplitude?
  * Why is the anomaly quadratic in the angle for intermediate amplitude?
  * Why does the anomaly increase faster than angle squared for very large angle?

<iframe src="https://trinket.io/embed/glowscript/a00ebb6c42" width="100%" height="600" frameborder="0" marginwidth="0" marginheight="0" allowfullscreen></iframe>


