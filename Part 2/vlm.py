#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main driver showing how to model Unified Engineering Flight Competition wings
using a vortex lattice method

Author: David Darmofal, MIT
Date: April 2, 2021
"""

import UEFC_wing
import matplotlib.pyplot as plt

# Set-up a wing object. Inputs: wingspan, root chord, tip chord, root incidence
# angle, tip incidence angle, dihedral angle.
# Note: wing twist is defined as (agroot - agtip). Therefore, this wing has a
# twist of +5 degrees, and the tip is at a lower incidence angle than the root.
root_angle = 3.13
washout_diff = -5.0
tip_angle = root_angle + washout_diff
PV = UEFC_wing.UEFC_wing(b=1.153, croot=0.186, ctip=0.093, agroot=12., agtip=4.5, dihedral=10.)

# Plot the wing geometry
PV.plotgeom()

# Solve the flow around the wing at a desired CL
G, alpha = PV.solve(CL=0.65) # The solution is returned in G and the required angle of attack in alpha
# G, alpha = PV.solve(alpha=3.0) # If you wanted to solve at a desired angle of attack (in degrees)

# Get the aspect ratio and surface area
AR = PV.get_AR()
S  = PV.get_S()

print('AR = {:.2f}'.format(AR))
print('S  = {:.3f}'.format(S))
print()

# Get the cl vs. y (2D lift coefficient vs. wing span location) distribution
cl, y = PV.calccldist(G)

# Plot the cl distribution (needed to determine stall behavior)
PV.plotcl(G, plotclccbar=False)

# # Calculate key aerodynamic performance parameters using solution
CL, CDi, e0, clmax= PV.calc_aeroperf(G)

print('alpha  = {:.2f} deg'.format(alpha))   # Angle of attack
print('CL     = {:.2f}'.format(CL))          # 3D lift coefficient
print('clmax  = {:.2f}'.format(clmax))       # Maximum 2D lift coefficient
print('CDi    = {:.4f}'.format(CDi))         # Induced drag coefficient
print('e0     = {:.2f}'.format(e0))          # Span efficiency in level flight

plt.show()
