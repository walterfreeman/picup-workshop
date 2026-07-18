---
layout: page
permalink: orbits.html
category: top
title: Orbital dynamics 
use_math: true
---

## Overview

As an example of a computational exercise that you might do in an introductory course and that is calibrated to be completed without that much LLM help,
let's look at planetary orbits. This will be particularly interesting as a playground for how well system prompt guardrails work -- if you want AI
to help students with programming but explicitly not computational science or physics.

### Part 1: Building the model

Write some code that simulates and animates planetary orbits around a star. Instead of using SI units here, use sensible units for simulating things in our 
Solar System: AU for distance, years for time, and solar masses for mass.

In these units, what is the value of G?

Test your model by simulating an Earthlike orbit. (You will need to think about the initial conditions to use for Earth.) Does it behave like we know Earth's orbit does?

### Part 2: Testing Kepler's laws

Now verify Kepler's laws for your planet. Each law has a qualitative and a quantitative thing you can test. For the quantitative tests, make numerical measurements
on your simulation and determine to what precision your simulation follows Kepler's laws. (Why are they not followed perfectly?)

* Kepler's first law:
  * Qualitative: Are the orbits elliptical?
  * Quantitative: Do the orbits close precisely?
* Kepler's second law:
  * Qualitative: Does the planet go faster when it is closer to the Sun?
  * Quantitative: "Equal areas in equal times" is just a geometric way of phrasing the conservation of angular momentum. To what precision does your simulation conserve angular momentum?
* Kepler's third law:
  * Qualitative only: How closely do your orbits follow the square/cube law?


### Part 3: An extension

Replace your star with a binary star. Verify that a planet orbiting a binary star far away follows Kepler's laws, and explore their breakdown as the planet gets closer to the
stars. Can you quantify the extent to which Kepler's laws break down as the size of the planetary orbit decreases? How does the eccentricity of the stars' orbit and the planet's orbit affect this?
