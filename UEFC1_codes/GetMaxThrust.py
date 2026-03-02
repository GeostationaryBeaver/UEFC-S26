import numpy as np

def GetMaxThrust(UEFC, V):

    # YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
    rho = UEFC.rho

    ct0 =  0.0282
    ct1 =  -0.0573
    ct2 =  -0.2022

    Tmax_static = 2               # Maximum thrust desired at static conditions
    Rprop       = 0.1016          # Propeller radius (m)
    Aprop       = np.pi * (Rprop**2)          # Propeller disk area

    Omega  = np.sqrt(Tmax_static/(0.5*rho*Rprop**2*Aprop*ct0))

    # calculate Lambda, CT, and then Tmax from given variables
    Lambda = V/(Rprop * Omega)  # Advance ratio

    CT = ct0 + ct1*Lambda + ct2*(Lambda**2) # Thrust coefficient

    Tmax = CT * 0.5 * rho * (Omega * Rprop)**2 *Aprop

    return Tmax

# DO NOT MODIFY THIS
def check_close(truth_val, test_val, close_tol):
    return np.abs(truth_val - test_val) < close_tol

def tests() -> None:
    # DO NOT CHANGE THE VALUES HERE
    from GetUEFC        import UEFC
    CLOSE_TOL = 1E-10
    aircraft = UEFC()
    aircraft.mpay_g = 300. # set payload mass (g)
    aircraft.tau = 0.12
    aircraft.taper = 0.7
    aircraft.dbmax = 0.10

    V = 0.0
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 2, CLOSE_TOL)

    V = 7.3
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 1.2895060476076137, CLOSE_TOL)

    V = 10.0
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 0.9182809157548276, CLOSE_TOL)

    V = 12.0
    Tmax = GetMaxThrust(aircraft, V)
    assert check_close(Tmax, 0.6055462567668711, CLOSE_TOL)

    print(f"==> All GetMaxThrust tests have passed!")

if __name__ == "__main__":
    tests()