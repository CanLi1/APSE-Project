Bubbling Fluidized Bed - DAE - Model
====================================
*by David M Thierry*
*July 2016*

This folder contains models for the BFB written for Pyomo.
This model was based on the work by Minghzao Yu on:

> Yu, Mingzhao, David C. Miller, and Lorenz T. Biegler. "Dynamic reduced order models for simulating bubbling fluidized bed adsorbers." Industrial & Engineering Chemistry Research 54.27 (2015): 6959-6974.

The model is intended to be used in optimization problems. Discretization on space done by Orthogonal collocation. Non-linear system of equations solves for a dummy objective function. It uses a previously known solution as initial guess.

Main files:

	* `ref_ss_2.py` : the main model
	* `ss_sript_v0_test.py` : the script that runs the model
	* `dat_ss_1.dat` : the data file that contains the initial guess
	* `cpoinsc.py` : generates the roots of orthogonal polynomials
