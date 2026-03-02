def GetWpay(UEFC, opt_vars, AR, S):

    # Return payload weight in N

    mpay_g = UEFC.mpay_g  # Mass of payload (alone) in grams

    Wpay = UEFC.g * mpay_g / 1000

    return Wpay
