import numpy as np

def GetCDfuse(UEFC, opt_vars, AR, S):

    # You need to finish this file

    # Fuselage drag model including S-dependence (calibrated to Plane Vanilla)
    SPV       = 0.225
    CDA_fuse0 = 0.002
    CDA_fuseS = 0.002

    # calculate CDfuse from given variables
    # CDA0   = np.nan
    CDfuse = 1/S * (CDA_fuse0 + CDA_fuseS *(S/SPV))

    return CDfuse

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    aircraft = UEFC()
    aircraft.mpay_g = 300. # set payload mass (g)
    aircraft.tau = 0.12
    aircraft.taper = 0.7
    aircraft.dbmax = 0.10
    AR = 10
    S = 0.225
    opt_vars = np.array([1.1]) # load factor
    CDfuse = GetCDfuse(aircraft, opt_vars, AR, S)
    CLOSE_TOL = 1E-10
    assert check_close(CDfuse, 0.017777777777777778, CLOSE_TOL)

    aircraft = UEFC()
    aircraft.mpay_g = 280. # set payload mass (g)
    aircraft.tau = 0.10
    aircraft.taper = 0.4
    aircraft.dbmax = 0.05
    AR = 11
    S = 0.7
    opt_vars = np.array([1.06]) # load factor
    CDfuse = GetCDfuse(aircraft, opt_vars, AR, S)
    assert check_close(CDfuse, 0.011746031746031746, CLOSE_TOL)
    print(f"==> All GetCDfuse tests have passed!")

if __name__ == "__main__":
    tests()