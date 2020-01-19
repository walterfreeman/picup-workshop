---
layout: page
permalink: showcase.html
category: top
title: Showcase 
use_math: true
---

Here are links to some of the upper-division projects we've showed. Walter's sample codes use the `anim` animation tool, available for Linux or Mac <a href="https://walterfreeman.github.io/anim.tar">here</a>.

### Quantum wave packet

When I (Walter) was an undergraduate I never quite understood how the Schroedinger equation replicated classical behavior. Sure, I wrote the equations down for expectation values and so forth, but I wish I'd gotten to do this assignment (written by Larry):

The assignment <a href="wavepacket-assignment.pdf">[PDF]</a>

<a href="wavepacket-template.html">A viewable HTML template</a> (for those without Jupyter)

<a href="wavepacket.ipynb">Jupyter notebook</a>

### Nonlinear vibrating string (uses `anim`)

In the computational physics course at Syracuse, the nonlinear swinging pendulum is just a warmup for the nonlinear vibrating string. It turns out there's a lot of physics here that most of us don't get a chance to study -- and the musicians don't know either!

This project is actually online as a <a href="https://www.compadre.org/PICUP/exercises/exercise.cfm?I=151&A=vibrating-string">PICUP Exercise Set</a>. Or, you can see the version I've assigned my students. It's online in two parts -- see <a href="https://walterfreeman.github.io/phys307/projects/string-1.pdf">part 1</a> and <a href="https://walterfreeman.github.io/phys307/projects/string-2.pdf">part 2</a>.

Here is the C++ code that shows <a href="string-parallel.cpp">many strings oscillating at once</a>, or the simpler (and faster) C code that animates only <a href="string-plot.c">one at a time</a>; these generate the plots of period shift vs. amplitude.

Here's some Python code that simulates <a href="stringparallel.py">multiple strings at once.</a> It doesn't measure the period, but is an illustration of how to do this with NumPy.

### Heat conduction

Here's Kelly Roos' simulation of <a href="heat_equation_2D_source.m">heat conduction in 2D</a> in MATLAB.

### Thermodynamics of the Lennard-Jones gas(uses `anim`)

To compile the Lennard-Jones gas simulation, you'll need two files:

* the <a href="gas3d.cpp">main code</a>
* a small <a href="vector.h">vector library header file</a>

Download these and put them in the same directory, then compile with `g++` as usual.
