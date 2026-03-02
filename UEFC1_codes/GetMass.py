def GetMass(UEFC, opt_vars, AR, S):
    
    # Returns the total mass, as well as a breakdown, in g. 

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM 
    
    g           = UEFC.g
    weight_data = UEFC.weight(opt_vars, AR, S)
    
    mass_data = {
            "Total": weight_data["Total"] * 1000 / g,
            "Breakdown": {}
            }
    
    for element, weight in weight_data["Breakdown"].items():
        mass_data["Breakdown"][element] = weight * 1000 / g

    return mass_data
    