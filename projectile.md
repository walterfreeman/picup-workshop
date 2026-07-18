---
layout: page
permalink: projectile.html
category: top
title: Projectile motion 
use_math: true
---

## Projectile motion with and without drag 

This is an introductory exercise that beginning students might do with little LLM assistance.

## Introduction

Consider a projectile launched at an angle $\theta$ above the horizontal at an initial velocity $v_0$ on flat ground. 

1. Write code that simulates and visualizes its trajectory without drag, and compare with the analytical result for how far it travels before 
landing back on the ground. Does your simulation reproduce the textbook result $d = \frac{v_0^2 \sin \theta \cos \theta}{2g}$?

2. Now add a quadratic drag force $F_{\rm drag} = (-\hat v) \beta v^2$, where $\beta = \frac{1}{2} \rho C_D A$; $\rho$ is the density of air, $C_D$ is a dimensionless
drag coefficient that depends only on the shape of the object, and $A$ is the object's cross-sectional area.

Consider a soccer ball kicked at the World Cup. Make a plot of range vs. initial velocity with and without drag for a variety of initial angles.

3. Now consider what happens if there is a headwind of $v_w = 8 \rm m/\rm s$. Show with your simulation data that there is some truth behind the 
traditional wisdom that it is better to keep the ball "low" when playing into the wind.

4. Now consider a more complicated problem: trajectories that extend high enough into the atmosphere that $\rho$ is not constant. Consider a projectile with
   $v_0 = 800 \rm m/\rm s$, $m = 50$ kg, $A = 0.075 \rm m^2$, and $C_D = 0.2$. (This is roughly characteristic of the artillery used in Ukraine.) Determine the 
maximum range of these projectiles with an atmosphere model in which $\rho = \rho_0 e^{-h/h_0}$ and $h_0$ is 10 km. Also determine the error made in the distance
by not taking into account the changing density of the atmosphere at high altitude.
