import numpy as np

def GetRequiredThrust(UEFC, opt_vars, AR, S):

    # Calculate the required thrust from UEFC parameters and opt_vars, AR, S
    rho = UEFC.rho

    # Calculate CD
    CD = UEFC.drag_coefficient(opt_vars, AR, S)["Total"]

    # Calculate q (dynamic pressure)
    V = UEFC.flight_velocity(opt_vars, AR, S)
    q = 0.5 * rho * V * V

    # Calculate required thrust from CD, q, S
    T = CD * q * S

    return T

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    CLOSE_TOL = 1E-10
    aircraft = UEFC()
    aircraft.rho = 1.225     # air density kg/m^3


    S = 0.354
    AR = 9
    N = 1.1
    opt_vars = np.array([N])
    Treq = GetRequiredThrust(aircraft, opt_vars, AR, S)
    assert check_close(Treq, 0.6051271411702426, CLOSE_TOL)

    S = 0.9
    AR = 12
    N = 1.05
    opt_vars = np.array([N])
    Treq = GetRequiredThrust(aircraft, opt_vars, AR, S)
    assert check_close(Treq, 0.9285875411863036, CLOSE_TOL)

    print(f"==> All GetRequiredThrust tests have passed!")


if __name__ == "__main__":
    tests()
