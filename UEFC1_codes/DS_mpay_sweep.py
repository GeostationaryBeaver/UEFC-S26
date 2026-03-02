from GetUEFC import UEFC
import numpy as np
import os
from DS_report_opt_obj import report_opt_obj
from opt_obj import opt_obj
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.style.use(os.path.join(os.path.dirname(__file__), "uefc.mplstyle"))

# You should not need to change the values within this function
def mpay_sweep(aircraft: UEFC,
               AR: float,
               S: float,
               mpay_start: float = 0.0,
               mpay_end: float = 0.0,
               mpay_num: int = 0,
               show_plot: bool = True,
               ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    assert mpay_num >= 0, "mpay_num must be an integer greater than or equal to 0"
    if mpay_num == 0:
        mpay_array = np.array([mpay_start])
    else:
        mpay_array = np.linspace(mpay_start, mpay_end, mpay_num)

    obj_array = []
    mpayout_array = []
    CL_array = []
    CD_array = []
    T_req_array = []
    T_max_array = []
    db_array = []
    N_array = []
    for j, mpay in enumerate(mpay_array):
        aircraft.mpay_g = mpay
        opt_vars_maxObj, obj_max, success = opt_obj(aircraft, AR, S)
        if success:
            CL = aircraft.lift_coefficient(opt_vars_maxObj, AR, S)
            CD = aircraft.drag_coefficient(opt_vars_maxObj, AR, S)["Total"]
            obj = obj_max
            T_req = aircraft.required_thrust(opt_vars_maxObj, AR, S)
            T_max = aircraft.maximum_thrust(obj)
            db = aircraft.wing_tip_deflection(opt_vars_maxObj, AR, S)
            N_array.append(opt_vars_maxObj[0])
            obj_array.append(obj)
            mpayout_array.append(mpay)
            CL_array.append(CL)
            CD_array.append(CD)
            T_req_array.append(T_req)
            T_max_array.append(T_max)
            db_array.append(db)
        else:
            print(f"WARNING: For mpay = {mpay:.1f} g, the optimizer did not converge. This likely means that current plane cannot fly with the requested payload.")
    # Cast to a numpy array
    obj_array = np.array(obj_array)
    mpayout_array = np.array(mpayout_array)
    CL_array = np.array(CL_array)
    CD_array = np.array(CD_array)
    T_req_array = np.array(T_req_array)
    T_max_array = np.array(T_max_array)
    db_array = np.array(db_array)
    N_array = np.array(N_array)

    # plotting
    if show_plot:
        # plotting
        plt.close()
        markers = ['x', '*', '.', '+', "^"]
        tab10_colors = cm.tab10.colors
        #colors = ["blue", "red", "green", "msgenta", "cyan"]
        #label_base = '$(\\delta/b)_{{\\mathrm{max}}}'

        fig, ax = plt.subplots(3, 2, figsize=(9, 9))

        # Objective function: velocity
        ax[0,0].plot(mpayout_array, obj_array, marker=markers[0], color=tab10_colors[0])
        ax[0,0].grid(True)
        #ax[0,0].set_xlabel(f"Payload mass $m_{{\\mathrm{{pay}}}}$ [g]")
        ax[0,0].set_ylabel(f"Velocity $V$ [m/s]")

        # Wing tip deflection
        ax[0,1].plot(mpayout_array, db_array, marker=markers[4], color=tab10_colors[4] )
        ax[0,1].axhline(y=aircraft.dbmax, color='r', linestyle='--')
        ax[0,1].grid(True)
        #ax[0,1].set_xlabel(f"Payload mass $m_{{\\mathrm{{pay}}}}$ [g]")
        ax[0,1].set_ylabel(f"Tip bending $\\delta/b$")

        # Lift coefficient
        ax[1,0].plot(mpayout_array, CL_array, marker=markers[1], color=tab10_colors[1] )
        ax[1,0].axhline(y=aircraft.CLdes, color='r', linestyle='--')
        ax[1,0].grid(True)
        #ax[1,0].set_xlabel(f"Payload mass $m_{{\\mathrm{{pay}}}}$ [g]")
        ax[1,0].set_ylabel(f"Lift coefficient $C_L$")

        # Drag coefficient
        ax[1,1].plot(mpayout_array, CD_array, marker=markers[1], color=tab10_colors[1] )
        ax[1,1].grid(True)
        #ax[1,1].set_xlabel(f"Payload mass $m_{{\\mathrm{{pay}}}}$ [g]")
        ax[1,1].set_ylabel(f"Drag coefficient $C_D$")

        # Thrust
        ax[2,0].plot(mpayout_array, T_max_array, color=tab10_colors[2], marker=markers[2], label="Thrust Max", linestyle="solid")
        ax[2,0].plot(mpayout_array, T_req_array, color=tab10_colors[3], marker=markers[3], label="Thrust Required", linestyle="dotted")
        ax[2,0].legend()
        ax[2,0].grid(True)
        ax[2,0].set_xlabel(f"Payload mass $m_{{\\mathrm{{pay}}}}$ [g]")
        ax[2,0].set_ylabel(f"Thrust [N]")

        # Load Factor
        ax[2,1].plot(mpayout_array, N_array, marker=markers[0], color=tab10_colors[0] )
        ax[2,1].grid(True)
        ax[2,1].set_xlabel(f"Payload mass $m_{{\\mathrm{{pay}}}}$ [g]")
        ax[2,1].set_ylabel(f"Load factor [-]")

        suptitle = f"$AR = {AR:.1f}$, $S = {S:.3f}$ m$^2$, \n $C_{{L_{{\\mathrm{{des}}}}}} = {aircraft.CLdes:.2f}$, $\\lambda = {aircraft.taper:.2f}$, $\\tau = {aircraft.tau:.2f}$, $(\\delta/b)_{{\\mathrm{{max}}}} = {aircraft.dbmax:.2f}$"
        fig.suptitle(suptitle)
        plt.show()
    return mpayout_array, obj_array, CL_array, CD_array, T_req_array, T_max_array, db_array, N_array

if __name__ == "__main__":

    # Simple test case. Feel free to modify this part of the file.

    # Payload
    mpay_start = 0.0  # g
    mpay_end   = 500  # g
    mpay_num   = 25

    aircraft = UEFC()

    # Geometry parameters
    S  = 0.354                 # Wing area (m^2)
    AR = 9                 # Wing aspect ratio
    aircraft.taper    = 0.45  # taper ratio
    aircraft.dihedral = 10 # Wing dihedral (degrees)
    aircraft.tau      = .10  # thickness-to-chord ratio

    # Tail parameters
    aircraft.Sh = 0.04 # Wing area of horizontal tail (m^2)
    aircraft.Sv = 0.03 # Wing area of vertical tail (m^2)

    # Fuselage parameters
    aircraft.l_AR = 1.63 # Fuselage length to wingspan ratio (-)

    # Aerodynamic parameters
    aircraft.CLdes = 0.75  # maximum CL wing will be designed to fly at (in cruise)
    aircraft.e0    = 1.00  # Span efficiency for straight level flight

    # Wing bending and material properties
    aircraft.dbmax   = 0.06  # tip displacement bending constraint
    aircraft.rhofoam = 32.     # kg/m^3. high load foam
    aircraft.Efoam   = 19.3E6  # Pa.     high load foam


    mpay, obj, CL, CD, T_req, T_max, db, N = mpay_sweep(aircraft,
                                                        AR, S,
                                                        mpay_start=mpay_start,
                                                        mpay_end=mpay_end,
                                                        mpay_num=mpay_num,
                                                        show_plot=True)