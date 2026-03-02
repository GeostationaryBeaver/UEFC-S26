import numpy as np

def GetV(UEFC, opt_vars, AR, S):

    # Calculate speed from N, g, R during a level turn.
    # N = load factor = Lift / Weight.
    g = UEFC.g
    N = opt_vars[0]
    R = UEFC.R

    V = np.nan  # Holds in a level turn

    V = np.sqrt(R*g*np.sqrt(N**2-1))

    return V

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    CLOSE_TOL = 1E-10
    aircraft = UEFC()
    aircraft.g = 9.81

    aircraft.R = 12.5

    S = 0.225
    AR = 10
    N = 1.5
    opt_vars = np.array([N])
    V = GetV(aircraft, opt_vars, AR, S)
    assert check_close(V, 11.708924710256525, CLOSE_TOL)

    S = 0.9
    AR = 12
    N = 1.01
    opt_vars = np.array([N])
    V = GetV(aircraft, opt_vars, AR, S)
    assert check_close(V, 4.16954364829975, CLOSE_TOL)

    print(f"==> All GetV tests have passed!")


if __name__ == "__main__":
    tests()