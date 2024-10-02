### This module contains date handling routines used in the other ecco_access modules.


import numpy as np


def date_adjustment(ShortName,StartDate,EndDate,CMR_query=True):
    """
    Adjusts StartDate and EndDate, augmenting where day or month may be missing.
    Returns either strings ready for NASA Earthdata CMR query (CMR_query = True), 
    or numpy.datetime64 values.
    """
    
    pass
    
    
    # # Adjust StartDate and EndDate

    if StartDate=='yesterday':
        StartDate = yesterday()
    if EndDate==-1:
        EndDate = StartDate
    elif StartDate=='yesterday':
        StartDate = yesterday()
    elif EndDate=='today':
        EndDate = today()

    if len(StartDate) == 4:
        StartDate += '-01-01'
    elif len(StartDate) == 7:
        StartDate += '-01'
    elif len(StartDate) != 10:
        sys.exit('\nStart date should be in format ''YYYY'', ''YYYY-MM'', or ''YYYY-MM-DD''!\n'\
                 +'Program will exit now !\n')

    if len(EndDate) == 4:
        EndDate += '-12-31'
    elif len(EndDate) == 7:
        EndDate = str(np.datetime64(str(np.datetime64(EndDate,'M')+np.timedelta64(1,'M'))+'-01','D')\
                      -np.timedelta64(1,'D'))
    elif len(EndDate) != 10:
        sys.exit('\nEnd date should be in format ''YYYY'', ''YYYY-MM'', or ''YYYY-MM-DD''!\n'\
                 +'Program will exit now !\n')
    
    # for snapshot datasets, move EndDate one day later
    if 'SNAPSHOT' in ShortName:
        EndDate = str(np.datetime64(EndDate,'D') + np.timedelta64(1,'D'))
    
    # CMR request adjustments
    if CMR_query:
        SingleDay_flag = False
        if (('MONTHLY' in ShortName) or ('DAILY' in ShortName)):
            if np.datetime64(EndDate,'D') - np.datetime64(StartDate,'D') \
              > np.timedelta64(1,'D'):
                # for monthly and daily datasets, do not include the month or day before
                StartDate = str(np.datetime64(StartDate,'D') + np.timedelta64(1,'D'))
            else:
                # for single day ranges we need to make the adjustment
                # after the CMR request
                SingleDay_flag = True
        
        return StartDate,EndDate,SingleDay_flag
    
    else:
        StartDate = np.datetime64(StartDate,'D')
        EndDate = np.datetime64(EndDate,'D')
        
        return StartDate,EndDate
