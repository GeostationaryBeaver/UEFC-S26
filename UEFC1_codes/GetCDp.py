import numpy as np

def GetCDp(UEFC, opt_vars, AR, S):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
    tau      = UEFC.tau
    dihedral = (np.pi/180) * UEFC.dihedral  # Convert to radians

    # Determine profile drag (function of CL, tau, Re)
    cd0    = 0.020*(1+tau**2)
    cd1    = -0.005/(1+6*tau)
    cd2    = 0.160/(1+60*tau)
    cd8    = 1.0
    cl0    = 1.25 - 3*tau
    Re_ref = 1E5
    Re_a   = -0.75

    CL   = UEFC.lift_coefficient(opt_vars, AR, S)
    V    = UEFC.flight_velocity(opt_vars, AR, S)
    cbar = UEFC.wing_dimensions(AR, S)["Mean chord"]
    Re   = (UEFC.rho * V * cbar) / UEFC.mu  # Reynolds number

    cl2d   = CL/np.cos(dihedral)
    cdpfac = cd0 + cd1*(cl2d-cl0) + cd2*(cl2d-cl0)**2 + cd8*(cl2d-cl0)**8
    CDp    = (cdpfac*(Re/Re_ref)**Re_a)/np.cos(dihedral)

    return CDp

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
    CDp = GetCDp(aircraft, opt_vars, AR, S)
    CLOSE_TOL = 1E-10
    assert check_close(CDp, 0.02500368734511909, CLOSE_TOL)

    aircraft = UEFC()
    aircraft.mpay_g = 280. # set payload mass (g)
    aircraft.tau = 0.10
    aircraft.taper = 0.4
    aircraft.dbmax = 0.05
    AR = 11
    S = 0.7
    opt_vars = np.array([1.06]) # load factor
    CDp = GetCDp(aircraft, opt_vars, AR, S)
    assert check_close(CDp, 0.021979298341213087, CLOSE_TOL)

    print(f"==> All GetCDp tests have passed!")

if __name__ == "__main__":
    tests()