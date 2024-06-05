# Low dimensional model of plane Couette flow

This code performs modal reduction and Galerkin projection of plane couette flow (PCF). This code is an extension for formulation in Zhao. et al (2021) to more than 1 POD mode for each basis. 

## Problem now
 - In modal decomposition, one must mannually set the number of wavenumber in span and streamwise direction.
 - In moal decomposition, the code is not aware how many dns file (i.e xx.nc) are there to load. That need to be mannually set.

## File prefixes
- t_: test file/scripts
- verify_: verification file/script to test code
- vis_: data/scripts for visualisation only
- dns_: dns related data/scripts
- gp_: Galerkin projection related data/scripts
