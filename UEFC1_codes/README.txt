This file describes how to use the Python aircraft sizing for the UEFC Project. While you can certainly look at any of the provided Python functions, the functions you are most likely to use are briefly described below.

GetUEFC: Contains the class definition for the UEFC aircraft. Most of its attached methods take in a 1-element array called opt_vars, which represents the optimization variables.
	opt_vars[0]: Load factor (-)

GetCDpay: Sets the payload-dependent drag coefficient increment.  Currently, this is set to zero (i.e. there is no drag caused by the payload).  Clearly, this is almost certainly incorrect and you can include a payload drag increment here.

GetWfuse: Calculates fuselage weight. The constants in here may need to be adjusted to fit your estimated airplane.

GetCDfuse: Calculates fuselage drag increment. You will need to implement this.

GetCDi: Calculates induced drag increment. You will need to implement this.

GetCL: Calculates lift coefficient for steady circular flight. You will need to implement this.

GetMaxThrust: Calculates speed-dependent max available thrust. You will need to implement this.

scan_ARS(UEFC, AR_start, AR_end, S_start, S_end, num_division, show_plots): This function will search over (AR,S) to determine optimal design velocity for each AR,S considered.  Then, it plots contours of velocity,  R, CL, T, Tmax, and d/b. An asterisk (*) is placed at the location in (AR,S) which has the highest obj.  Finally, scan_ARS also prints out what the performance, operating conditions, weight breakdowns, etc are for this highest obj aircraft. It can also be called as a script.

opt_obj(UEFC,AR,S): this function determines the maximum objective (velocity) achievable for an airplane with the inputted values of AR, S.  This function is called repeatedly by scan_ARS as it scans over (AR,S).

report_opt_obj(UEFC,AR,S): this function is a wrapper for opt_obj.  Calling it will printout the optimized performance, operating conditions, etc found after running opt_obj.  It calls opt_obj for you and then prints out useful information.

mpay_sweep(UEFC, AR, S, mpay_start, mpay_end, mpay_num, show_plot): This is a function that will search for the optimum plane at a given mpay. It can also be called as a script.