# Packet Speed and Cost in Mobile Wireless Delay-Tolerant Networks
#### by Riccardo Cavallari, Stavros Toumpis, Roberto Verdone, Ioannis Kontoyiannis

This repository contains the C source code of the simulator and the MATLAB scripts used in the paper "Packet Speed and Cost in Mobile Wireless Delay-Tolerant Networks". The paper can be found on arXiv at the following link https://arxiv.org/abs/1801.02229.

See the paper for details.

```
sim.exe takes 13 parameters:
    lambda:         density of the PPP.
    v0:             speed of the nodes.
    r0:             rate at which the nodes change direction.
    tw:             Theta_w (parameter of f_D(x)). (Anything between 1.57 and 1.58 translates to pi/2.)
    simDur:         duration of the simulation, can be in simulation time or in stages (see stageOrTime)
    dT:             simulation time step. E.g., if simDur is time, simulation will last simDur/dT steps.
    potential:      potential function. "a" for the only potential function implemented so far.
    boundary:       shape of the forwarding region. "d" for disk, "e" for ellipse
    a:              1st parameter of the forwarding region. It is the radius for the disk; the half-axis length for the ellipse.
    b               2nd parameter of the forwarding region. Not used for the disk; it is the eccentricity for the ellipse.
    stageOrTime:    if 0 simuDur represents time, if 1 represents number of stages.
    log_stages:     1 to produce log_stages.txt, 0 to not produce it. This file logs all the quantities of each stage (can get huge if there are a lot of stages)
    seedRand: -1 to use time(NULL), a non-negative number to use it as seed (for debugging).
```   
