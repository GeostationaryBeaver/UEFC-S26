import numpy as np
def GetCL(UEFC, opt_vars, AR, S):

    # You need to finish this file

    # Calculate the lift coefficient from UEFC parameters and opt_vars, AR, S
    rho = UEFC.rho
    N   = opt_vars[0]

    V  = UEFC.flight_velocity(opt_vars, AR, S)
    W  = UEFC.weight(opt_vars, AR, S)["Total"]

    # calculate CL from given variables
    CL = N*W / (0.5 * rho * V**2 * S)

    return CL

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    CLOSE_TOL = 1E-10
    aircraft = UEFC()
    aircraft.mpay_g = 250. # set payload mass (g)
    aircraft.tau = 0.12
    aircraft.taper = 0.7
    aircraft.dbmax = 0.10
    AR = 9
    S = 0.354
    N = 1.1
    opt_vars = np.array([N])
    CL = GetCL(aircraft, opt_vars, AR, S)
    assert check_close(CL, 0.8023760145690171, CLOSE_TOL)

    N = 1.05
    AR = 12
    S = 0.9
    opt_vars = np.array([N])
    CL = GetCL(aircraft, opt_vars, AR, S)
    assert check_close(CL, 0.7843885762197533, CLOSE_TOL)

    print(f"==> All GetCL tests have passed!")


if __name__ == "__main__":
    tests()