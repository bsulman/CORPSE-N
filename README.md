# CORPSE-mycorrhizae

This repository contains code for running the carbon-nitrogen version of the "Carbon Organisms Rhizosphere and Protection in the Soil Environment" (CORPSE) model of plant-microbe-soil interactions. The goal is to implement mycorrhizal processes including plant allocation of carbon to symbiotic microbes as well as mycorrhizal transfers of nutrients to plants. This work is supported by a grant from the U.S. Department of Energy Terrestrial Ecosystem Science Program (PI: Caitlin Hicks Pries; Modeling lead: Benjamin Sulman).

## Files:
`CORPSE_deriv.py`: Code defining the CORPSE model. The function `CORPSE_deriv` calculates the rate of change of all model pools given the current state of SOM pools along with parameter values, temperature, and moisture. The file also includes a list and description of required parameters and a utility function for adding up carbon/nitrogen over model pools.

`CORPSE_integrate.py`: Code for integrating CORPSE simulations over time. Currently includes two alternative methods for integrating the model over time. `run_CORPSE_ODE` uses the python ODE solver and `run_CORPSE_iterator` uses explicit time stepping. Both these approaches are doing a simple one-compartment model that does not separate the rhizosphere and does not include plant interactions. This file includes some example runs using both methods at the end which demonstrate that the two approaches produce the same results.

`gradient_sim.py`: Simulations used to generate graphs of modeled soil C and N pools across gradients of mycorrhizal associations, temperature, and soil clay content used in the original proposal.

## Requirements:
The code has been tested with python 3.9.5, numpy 1.20.3, matplotlib 3.4.2, pandas 1.2.4, and scipy 1.6.3