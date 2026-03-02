import numpy as np

def Getspaneff(UEFC, opt_vars, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
    
    # Get span efficiency including effects of sideslip due to flying
    # in a circular path
    
    N = opt_vars[0]
    R = UEFC.R
    b = UEFC.wing_dimensions(AR, S)["Span"]
    
    rbar = 0.5*b/R/N
    
    dihedral = (np.pi/180) * UEFC.dihedral  # Convert to radians
    
    CL = UEFC.lift_coefficient(opt_vars, AR, S)
    
    beta = CL/dihedral*(1+4/AR)/(2*np.pi)*rbar
    
    spaneff = UEFC.e0 * (1-0.5*rbar**2)*(np.cos(beta))**2
    
    return spaneff
