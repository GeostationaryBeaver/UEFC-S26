def Getepsilon(UEFC):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
    
    # Determine max camber/chord as a function of max thickness/chord
    # for the assumed RC plane airfoil shapes
    tau     = UEFC.tau
    epsilon = 0.1 - 0.5*tau
    
    return epsilon
