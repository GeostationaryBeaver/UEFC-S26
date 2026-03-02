import numpy as np

def GetWingWeight(UEFC, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
    
    # Calculate the wing weight from UEFC parameters and S and AR
    tau      = UEFC.tau
    l        = UEFC.taper
    g        = UEFC.g
    rhofoam  = UEFC.rhofoam
    dihedral = (np.pi/180) * UEFC.dihedral  # Convert to radians
    
    Afac = 0.66 # Area of airfoil assumed to be Afac*tau*c^2
    
    Wwing = (4/3)*Afac*rhofoam*g*tau/np.cos(dihedral)*S**1.5*AR**(-0.5)*(l**2+l+1)/(l+1)**2
    
    return Wwing
