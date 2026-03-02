from GetV import GetV

def GetObjective(UEFC, opt_vars, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FUNCTION FOR THIS PROBLEM

    # Calculate the objective function (payload mass x turn rate^3) in g/s^3

    # mpay  = opt_vars[2]
    # Omega = UEFC.turn_rate(opt_vars, AR, S)
    # b = UEFC.wing_dimensions(AR, S)['Span']
    #mwing = UEFC.wing_weight(AR,S) * 1000 / UEFC.g

    #obj = mpay * Omega / b.

    # Calculate the objective function (velocity) in m/s
    # print("objective called")
    obj = GetV(UEFC, opt_vars, AR, S)

    return obj



