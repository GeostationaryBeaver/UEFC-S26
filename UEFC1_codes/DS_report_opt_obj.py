#  YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM
#  BESIDES SCRIPT PORTION (BELOW LINE 93)
import numpy as np
from GetUEFC import UEFC
from opt_obj import opt_obj

def report_opt_obj(aircraft, AR, S):

    # This function is a wrapper for opt_obj. Calling it will print out the
    # optimized performance, operating conditions, etc found after running
    # opt_obj.  It calls opt_obj and then prints out useful information.

    opt_vars, obj, success = opt_obj(aircraft, AR, S)

    if success:

        wing_dimensions = aircraft.wing_dimensions(AR, S)

        N     = opt_vars[0]
        V     = aircraft.flight_velocity(opt_vars, AR, S)
        mpay  = aircraft.mpay_g

        print()
        print("Results Summary from opt_obj\n")
        print("obj = V      = %0.2f m/s"   % obj)
        print("mpay         = %0.0f g"     % mpay)
        # print("R            = %0.2f m"     % R)
        # print("V            = %0.2f m/s"   % V)
        # print("Omega        = %0.2f rad/s" % Omega)
        print()

        print("Geometry")
        print("----------------------------------------------\n")
        print("AR      = %5.3f"       % AR)
        print("S       = %5.3f sq. m" % S)
        print("b       = %5.3f m"     % wing_dimensions["Span"])
        print("cbar    = %5.3f m"     % wing_dimensions["Mean chord"])
        print("cr      = %5.3f m"     % wing_dimensions["Root chord"])
        print("ct      = %5.3f m"     % wing_dimensions["Tip chord"])
        print("lambda  = %5.3f"       % aircraft.taper)
        print("tau     = %5.3f"       % aircraft.tau)
        print("eps     = %5.3f\n"     % aircraft.max_camber())

        mass_data = aircraft.mass(opt_vars,AR, S)

        print("Masses")
        print("----------------------------------------------\n")
        print("W/g     = %4.0f g"   % mass_data["Total"])
        print("Wfuse/g = %4.0f g"   % mass_data["Breakdown"]["Fuselage"])
        print("Wwing/g = %4.0f g"   % mass_data["Breakdown"]["Wing"])
        print("Wpay/g  = %4.0f g\n" % mass_data["Breakdown"]["Payload"])

        CL      = aircraft.lift_coefficient(opt_vars, AR, S)
        CD_data = aircraft.drag_coefficient(opt_vars, AR, S)
        e       = aircraft.span_efficiency( opt_vars, AR, S)

        print("Aerodynamic performance")
        print("----------------------------------------------")
        print("N       = %5.3f"   % N)
        print("CL      = %5.3f"   % CL)
        print("CLdes   = %5.3f"   % aircraft.CLdes)
        print("CD      = %5.3f"   % CD_data["Total"])
        print("CDfuse  = %5.3f"   % CD_data["Breakdown"]["Fuselage"])
        print("CDp     = %5.3f"   % CD_data["Breakdown"]["Wing"])
        print("CDi     = %5.3f"   % CD_data["Breakdown"]["Induced"])
        print("CDpay   = %5.3f"   % CD_data["Breakdown"]["Payload"])
        print("e0      = %5.3f"   % aircraft.e0)
        print("e       = %5.3f\n" % e)

        T_req = aircraft.required_thrust(opt_vars, AR, S)
        T_max = aircraft.maximum_thrust(V)

        print("Thrust")
        print("----------------------------------------------\n")
        print("T       = %5.3f N"   % T_req)
        print("Tmax    = %5.3f N\n" % T_max)

        db = aircraft.wing_tip_deflection(opt_vars, AR, S)

        print("Bending")
        print("----------------------------------------------\n")
        print("d/b     = %5.3f" % db)
        print("d/bmax  = %5.3f" % aircraft.dbmax)

    else:

        print("\nError in opt_obj: success = " + str(success))
        print("  Usually this is because the airplane could not fly while " +
              "meeting all constraints.\n")

    return


if __name__ == "__main__":

    # Simple test case. Feel free to modify this part of the file.
    aircraft = UEFC()

    # Payload weight
    aircraft.mpay_g   = np.nan     # payload weight in grams

    # Geometry parameters
    S  = np.nan                  # Wing area (m^2)
    AR = np.nan                  # Wing aspect ratio
    aircraft.taper    = np.nan   # taper ratio
    aircraft.dihedral = np.nan   # Wing dihedral (degrees)
    aircraft.tau      = np.nan   # thickness-to-chord ratio

    # Tail parameters
    aircraft.Sh = np.nan # Wing area of horizontal tail (m^2)
    aircraft.Sv = np.nan # Wing area of vertical tail (m^2)

    # Fuselage parameters
    aircraft.l_AR = np.nan  # Fuselage wingspan to length ratio (-)

    # Aerodynamic parameters
    aircraft.CLdes    = np.nan    # maximum CL wing will be designed to fly at (in cruise)
    aircraft.e0       = np.nan    # Span efficiency for straight level flight

    # Wing bending and material properties
    aircraft.dbmax    = np.nan  # tip displacement bending constraint
    aircraft.rhofoam  = 32.     # kg/m^3. high load foam
    aircraft.Efoam    = 19.3E6  # Pa.     high load foam

    report_opt_obj(aircraft, AR, S)
