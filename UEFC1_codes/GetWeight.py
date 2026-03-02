def GetWeight(UEFC, opt_vars, AR, S):
    
    # Returns the total weight, as well as a breakdown, in N. 

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
    
    Wfuse = UEFC.fuselage_weight(AR, S)
    Wwing = UEFC.wing_weight(AR, S)
    Wpay  = UEFC.payload_weight(opt_vars, AR, S)
    
    W = Wfuse + Wwing + Wpay
    
    return {
            "Total":     W,
            "Breakdown": {
                    "Fuselage": Wfuse,
                    "Wing":     Wwing,
                    "Payload":  Wpay,
                    }
            }