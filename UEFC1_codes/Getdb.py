def Getdb(UEFC, opt_vars, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

    # Calculate the tip d/b from UEFC parameters and opt_vars, AR, S
    l       = UEFC.taper
    tau     = UEFC.tau
    Efoam   = UEFC.Efoam
    epsilon = UEFC.max_camber()
    Wfuse   = UEFC.fuselage_weight(AR, S)
    Wpay    = UEFC.payload_weight(opt_vars, AR, S)
    N       = opt_vars[0]

    db = (0.018*N*(Wfuse+Wpay)/(Efoam*tau*(tau**2+0.7*epsilon**2))*
         (1+l)**3*(1+2*l)*AR**3/S)

    return db
