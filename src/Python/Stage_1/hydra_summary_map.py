#=======================================================================
# This function creates a map between column index and the appropriate
# quantity of interest in summary.dat
#=======================================================================

def hydra_summary_map():

    data        = {}

    idval               = 0
    data['id']          = idval; idval = idval+1;       # 1      
    data['DL']          = idval; idval = idval+1;       # 2
    data['sigma']       = idval; idval = idval+1;       # 3
    data['Nb']          = idval; idval = idval+1;       # 4
    data['nub']         = idval; idval = idval+1;       # 5
    data['omgrat']      = idval; idval = idval+1;       # 6
    data['vtip']        = idval; idval = idval+1;       # 7
    data['wing_AR']     = idval; idval = idval+1;       # 8
    data['wing_LF']     = idval; idval = idval+1;       # 9
    data['wing_CL']     = idval; idval = idval+1;       # 10
    data['Weight']      = idval; idval = idval+1;       # 11
    data['Power']       = idval; idval = idval+1;       # 12
    data['Radius']      = idval; idval = idval+1;       # 13
    data['Lift Offset'] = idval; idval = idval+1;       # 14
    data['wing_b']      = idval; idval = idval+1;       # 15
    data['RTFrac']      = idval; idval = idval+1;       # 16
    data['Fuel']        = idval; idval = idval+1;       # 17
    data['EmptyWt']     = idval; idval = idval+1;       # 18
    data['CTsigma']     = idval; idval = idval+1;       # 19
    data['Payload']     = idval; idval = idval+1;       # 20
    data['Cost']        = idval; idval = idval+1;       # 21
    data['valid']       = idval; idval = idval+1;       # 22

    return data 