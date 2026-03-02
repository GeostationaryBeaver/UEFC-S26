# YOU SHOULD NOT NEED TO CHANGE THIS FILE FOR THIS PROBLEM

# This script will search over (AR,S) to determine optimal designs for each
# AR,S considered. Then, it plots contours of obj, d/b, CL, T, and Tmax.
# An asterick (*) is placed at the location in (AR,S) which has the
# highest objective (obj).  Finally, scan_ARS also prints out what the
# performance, operating conditions, weight breakdowns, etc are for this
# optimized aircraft.

import numpy as np
import os
from matplotlib import pyplot as plt
plt.style.use(os.path.join(os.path.dirname(__file__), "uefc.mplstyle"))

# you may have to comment this line out if you are running from the command window
#from IPython    import get_ipython

from GetUEFC           import UEFC
from opt_obj           import opt_obj
from DS_report_opt_obj import report_opt_obj

def scan_ARS(aircraft: UEFC,
             AR_start: float = 5.0,
             AR_end: float = 15.0,
             S_start: float = 0.1,
             S_end: float = 0.7,
             num_division: int = 41,
             show_plots: bool = True,
             savefig: str = None
             ):
    nAR = nS = num_division
    ARarray = np.linspace(AR_start,   AR_end,  nAR)  # Aspect-ratio values
    Sarray  = np.linspace(S_start, S_end, nS)   # Surface-area values (m^2)

    ARvals, Svals = np.meshgrid(ARarray, Sarray, indexing="ij")  # 2D.

    # initialize 2D output arrays
    objvals   = np.zeros((nAR, nS))  # Objective function (g/s/m)
    mpayvals  = np.zeros((nAR, nS))  # Payload mass (g)
    pfvals    = np.zeros((nAR, nS))  # Payload fraction (-)
    Omegavals = np.zeros((nAR, nS))  # Turn rate (rad/s)
    bvals     = np.zeros((nAR, nS))  # Wingspan (m)
    Rvals     = np.zeros((nAR, nS))  # Turn radius (m)
    CLvals    = np.zeros((nAR, nS))  # Lift coefficient (-)
    CDvals    = np.zeros((nAR, nS))  # Drag coefficient (-)
    Tvals     = np.zeros((nAR, nS))  # Required thrust (N)
    Tmaxvals  = np.zeros((nAR, nS))  # Maximum thrust (N)
    dbvals    = np.zeros((nAR, nS))  # Wingtip deflection / wingspan
    Nvals     = np.zeros((nAR, nS))  # load factor

    # Sweep over (AR, S)
    for iAR,AR in enumerate(ARarray):
        for iS,S in enumerate(Sarray):

            # Determine max objective
            opt_vars, obj, success = opt_obj(aircraft, AR, S)

            if success:

                V = aircraft.flight_velocity(opt_vars, AR, S)

                objvals[iAR,iS]   = obj
                pfvals[iAR,iS]    = aircraft.payload_fraction(opt_vars, AR, S)
                Omegavals[iAR,iS] = aircraft.turn_rate(opt_vars, AR, S)
                bvals[iAR,iS]     = aircraft.wing_dimensions(AR, S)['Span']
                CLvals[iAR,iS]    = aircraft.lift_coefficient(opt_vars, AR, S)
                CDvals[iAR,iS]    = aircraft.drag_coefficient(opt_vars, AR, S)["Total"]
                Tvals[iAR,iS]     = aircraft.required_thrust(opt_vars, AR, S)
                Tmaxvals[iAR,iS]  = aircraft.maximum_thrust(V)
                dbvals[iAR,iS]    = aircraft.wing_tip_deflection(opt_vars, AR, S)
                Nvals[iAR,iS]     = opt_vars[0]

        if show_plots:
            print("Completed %3.1f%% of (AR,S) scan" % (100*(iAR+1)/nAR))

    # Find and print the optimal point (where objective is maximized)
    obj_opt           = np.max(objvals)
    (iAR_opt, iS_opt) = np.unravel_index(objvals.argmax(), objvals.shape)

    ARopt = ARvals[iAR_opt, iS_opt]
    Sopt  = Svals[iAR_opt, iS_opt]

    # # Plotting commands
    if show_plots:

        print("\n")
        print("Objective Maximizing aircraft characteristics:")
        print("----------------------------------------------")
        report_opt_obj(aircraft, ARopt, Sopt)

        # you may have to comment out these next 4 lines out if you are running from the command window
        # plt.ion()
        # plt.rc('axes', axisbelow=True)
        # plt.show()
        # get_ipython().run_line_magic('matplotlib', 'qt')
        marker=(8,2,0)  # 8-sided asterisk

        fig, ax = plt.subplots(3, 2, figsize=(9, 9)) #, dpi=150

        # Helper function for consistent contour plots
        def contourf(ax,vals,levels):
            cs = ax.contourf(ARvals, Svals, vals, levels=levels)
            cs.set_edgecolor("black")
            clabels = ax.clabel(cs, colors="black", fontsize=8.0, fmt="%.2f")
            [txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=0)) for txt in clabels]
            ax.plot(ARopt, Sopt, marker=marker, color="red", markersize=20)

        # Objective function: velocity
        levels = np.linspace(np.min(objvals[objvals > 0]), np.max(objvals), 9)
        contourf(ax[0,0],objvals,levels)
        #ax[0,0].set_xlabel("Aspect ratio (-)")
        ax[0,0].set_ylabel("Wing area ($m^2$)")
        ax[0,0].set_title("Objective = V (m/s)")

        # Objective function: payload fraction
        levels = np.linspace(np.min(pfvals[pfvals > 0]), np.max(pfvals), 9)
        contourf(ax[0,1],pfvals,levels)
        #ax[0,0].set_xlabel("Aspect ratio (-)")
        ax[0,1].set_title("Payload fraction (-)")

        # Lift coefficient
        levels = np.linspace(np.min(CLvals[CLvals > 0]), np.max(CLvals), 11)
        contourf(ax[1,0],CLvals,levels)
        #ax[1,0].set_xlabel("Aspect ratio (-)")
        ax[1,0].set_ylabel("Wing area ($m^2$)")
        ax[1,0].set_title("Lift coefficient (-)")

        # Drag coefficient
        levels = np.linspace(np.min(CDvals[CDvals > 0]), np.max(CDvals), 11)
        contourf(ax[1,1],CDvals,levels)
        #ax[1,1].set_xlabel("Aspect ratio (-)")
        #ax[1,1].set_ylabel("Wing area ($m^2$)")
        ax[1,1].set_title("Drag coefficient (-)")

        # Maximum thrust
        # Levels here are the same as for Tvals, so thast thrust contours can be
        # compared more easily
        levels = np.linspace(np.min(Tvals[Tvals > 0]), np.max(Tvals), 11)
        contourf(ax[2,0],Tmaxvals,levels)
        ax[2,0].set_xlabel("Aspect ratio (-)")
        ax[2,0].set_ylabel("Wing area ($m^2$)")
        ax[2,0].set_title("Maximum Thrust (N)")

        # Wing tip deflection
        levels = np.linspace(np.min(dbvals[dbvals > 0]), np.max(dbvals), 9)
        contourf(ax[2,1],dbvals,levels)
        ax[2,1].set_xlabel("Aspect ratio (-)")
        #ax[2,1].set_ylabel("Wing area ($m^2$)")
        ax[2,1].set_title("Wing Tip Deflection (-)")

        suptitle = f"$V_{{opt}} = {obj_opt:.2f}$ m/s, $AR_{{opt}} = {ARopt:.1f}$, $S_{{opt}} = {Sopt:.3f}$ m$^2$, \n $C_{{L_{{\\mathrm{{des}}}}}} = {aircraft.CLdes:.2f}$, $\\lambda = {aircraft.taper:.2f}$, $\\tau = {aircraft.tau:.2f}$, $(\\delta/b)_{{\\mathrm{{max}}}} = {aircraft.dbmax:.2f}$"
        fig.suptitle(suptitle)

        if savefig is not None:
            fig.savefig(savefig)

        plt.show()

    return obj_opt, ARopt, Sopt

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

    num_division = 41
    scan_ARS(aircraft, np.nan, np.nan, np.nan, np.nan, num_division, show_plots=True)




