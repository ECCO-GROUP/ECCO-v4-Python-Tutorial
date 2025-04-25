### This module contains routines to search the text of ECCO variable lists for matching datasets


import requests
import sys
from os.path import join


def varlist_file_parse(file):
    """Parse content of varlist file into dictionary"""

    response = requests.get(file)
    varlist_lines = response.text.split('\n')
    content_dict = {}
    curr_shortname = ''
    multiline_var = False
    for line in varlist_lines:
        if len(line) == 0:
            continue
        if (('ShortName' in line) and ('Variable Name' in line) and ('Description' in line)):
            varname_pos = line.index('Variable Name')
            descrip_pos = line.index('Description')
            varname_colwidth = descrip_pos - varname_pos
        elif 'ECCO_' in line:
            if '*' in line:
                first_ast_pos = line.index('*')
                curr_shortname = line[:first_ast_pos].replace(' ','')
                content_dict[curr_shortname] = {'Note':line[first_ast_pos:]}
            else:
                curr_shortname = line
                content_dict[curr_shortname] = {}
        elif 'Note:' in line:
            content_dict[curr_shortname]['Note'] = line[6:]
        elif len(curr_shortname) > 0:
            if len(line[:varname_pos].replace(' ','')) > 0:
                continue
            len_varname = line[varname_pos:].find(' ')
            if len_varname == -1:
                len_varname = len(line[varname_pos:])
            if len_varname >= varname_colwidth:
                varname = line[varname_pos:]
                multiline_var = True
                continue
            elif not multiline_var:
                varname = line[varname_pos:descrip_pos].replace(' ','')
            descrip = line[descrip_pos:]
            content_dict[curr_shortname][varname] = descrip
            if multiline_var:
                multiline_var = False

    return content_dict,varname_pos,descrip_pos



###================================================================================================================


def query_shortnames(query,content_dict,\
                     grid,test_grids,\
                     time_res,test_timeres):
    """Search varlist file for case-insensitive query results"""
    
    # make query case-insensitive
    query_ci = query.casefold()
    
    shortnames_match = []
    for shortname,curr_content in content_dict.items():

        # handle query results from the "mixed" grid and time res files
        # if single grid or time res was specified
        if ((grid is not None) and ('mixed' in test_grids)):
            if 'Note' in curr_content.keys():
                if '*' in curr_content['Note']:
                    if ((grid == 'native')\
                      and ('native' not in curr_content['Note'])):
                        continue
                    elif ((grid == 'latlon')\
                      and ('lat-lon' not in curr_content['Note'])):
                        continue
        if ((time_res != 'all') and ('all' in test_timeres)):
            if 'Note' in curr_content.keys():
                if '*' in curr_content['Note']:
                    if time_res not in curr_content['Note']:
                        continue

        # add query result to list
        if query_ci in shortname.casefold():
            shortnames_match.append(shortname)
        elif query_ci in str(curr_content).casefold():
            shortnames_match.append(shortname)
    
    return shortnames_match



###================================================================================================================


def print_varlist_query_results(query,shortnames_match,\
                                content_dict,varname_pos,descrip_pos):
    """Print text output of query"""

    varname_colwidth = descrip_pos - varname_pos
    query_resp_str = 'ShortName Options for query "'+query+'":\n'\
                        +(' '*varname_pos)\
                        +'Variable Name'.ljust(descrip_pos-varname_pos)\
                        +'Description (units)\n\n'
    for count,shortname in enumerate(shortnames_match):
        shortname_line = 'Option '+str(count+1)+': '+shortname
        if (('Note' in content_dict[shortname].keys())\
          and (' *' in content_dict[shortname]['Note'])):
            shortname_line += ('    '+content_dict[shortname]['Note'])
        else:
            shortname_line += '    *'
            if 'LLC' in shortname:
                shortname_line += 'native grid,'
            elif 'DEG' in shortname:
                shortname_line += 'lat-lon grid,'
            if 'MONTHLY' in shortname:
                shortname_line += 'monthly means*'
            elif 'DAILY' in shortname:
                shortname_line += 'daily means*'
            elif 'SNAPSHOT' in shortname:
                shortname_line += 'snapshots at daily intervals*'
        shortname_line += '\n'
        query_resp_str += shortname_line
        for varname,descrip in content_dict[shortname].items():
            max_descrip_len_perline = 50
            if varname == 'Note':
                max_descrip_len_perline += descrip_pos - 6
            if len(descrip) > max_descrip_len_perline:
                # put line breaks in long descriptions
                descrip_withbreaks = ''
                descrip_remaining = descrip + ''
                while len(descrip_remaining) > 0:
                    descrip_words = descrip_remaining[:(max_descrip_len_perline+1)].split(' ')
                    if len(descrip_remaining) <= max_descrip_len_perline:
                        curr_descrip = descrip_remaining
                    else:
                        curr_descrip = " ".join(descrip_words[:-1])
                    descrip_withbreaks += curr_descrip
                    descrip_remaining = descrip_remaining[(len(curr_descrip)+1):]
                    if len(descrip_remaining) > 0:
                        if varname == 'Note':
                            descrip_withbreaks += '\n      '
                        else:
                            descrip_withbreaks += ('\n'+(' '*descrip_pos))
                descrip = descrip_withbreaks
            if varname == 'Note':
                if ' *' in descrip:
                    continue
                else:
                    var_line = 'Note: '+descrip+'\n'
            elif len(varname) >= varname_colwidth:
                var_line = (' '*varname_pos)+varname+'\n'
                var_line += (' '*descrip_pos)+descrip+'\n'
            else:
                var_line = (' '*varname_pos)+(varname.ljust(varname_colwidth,' '))\
                            +descrip+'\n'
            query_resp_str += var_line
        query_resp_str += '\n'

    print(query_resp_str)
    
    return None



###================================================================================================================


def ecco_podaac_varlist_query(query,version,grid=None,time_res='all'):
    """
    
    This function queries the ECCO variable lists to find matches to "query" that match the grid and time_res types.
    The function then takes input from the user to determine which matching result to use.
    If no matches are found, then a ValueError raised.
    
    Parameters
    ----------    
    query: str, a text string being used to query ShortNames, variable names,
           and descriptions.
    
    version: ('v4r4','v4r5'), ECCO version to search variable lists for
    
    grid: ('native','latlon',None), specifies whether to query datasets with output
          on the native grid or the interpolated lat/lon grid.
          The default None will query both types of grids (and datasets with no spatial
          dimension), unless specified otherwise in a query dict (e.g., the example above).
    
    time_res: ('monthly','daily','snapshot','all'), specifies which time resolution 
              to include in query and downloads. 'all' includes all time resolutions, 
              and datasets that have no time dimension, such as the grid parameter 
              and mixing coefficient datasets.
    
    Returns
    -------
    shortname_match: str, ShortName matching the ECCO variable list query.
    
    """
    
    
    pass
    
    
    if version not in ['v4r4','v4r5']:
        raise ValueError('ECCO '+version+' is not currently available from PO.DAAC')
    
    # paths to variable list files
    varlist_url_root = 'https://raw.githubusercontent.com/ECCO-GROUP/ECCO-v4-Python-Tutorial/master/ecco_access/varlist/'
    if version == 'v4r4':
        varlist_url_ids = {'native,monthly':'v4r4_nctiles_monthly_varlist.txt',\
                           'native,daily':'v4r4_nctiles_daily_varlist.txt',\
                           'native,snapshot':'v4r4_nctiles_snapshots_varlist.txt',\
                           'latlon,monthly':'v4r4_latlon_monthly_varlist.txt',\
                           'latlon,daily':'v4r4_latlon_daily_varlist.txt',\
                           'mixed,all':'v4r4_tseries_grid_varlist.txt'}
    elif version == 'v4r5':
        varlist_url_ids = {'native,monthly':'v4r5_nctiles_monthly_varlist.txt'}
    
    
    # set keys of grid types and time resolutions to search
    
    if grid is None:
        test_grids = ['native','latlon','mixed']
        if time_res != 'all':
            test_timeres = [time_res,'all']
    else:
        test_grids = [grid]
    if time_res == 'all':
        if grid is not None:
            test_grids += ['mixed']
        test_timeres = ['monthly','daily','snapshot','all']
    elif grid is not None:
        test_timeres = [time_res]
    
    grid_timeres_keys = []
    for gridtype in test_grids:
        for timeres in test_timeres:
            if (((gridtype == 'mixed') and (timeres != 'all'))\
              or ((gridtype != 'mixed') and (timeres == 'all'))):
                # include time series/grid varlist
                grid_timeres_keys.append('mixed,all')
            elif ((gridtype == 'latlon') and (timeres == 'snapshot')):
                continue        
            else:
                grid_timeres_keys.append(gridtype+','+timeres)    
    
    # remove any duplicate keys
    keys_no_duplicates = []
    keys_no_duplicates = [key for key in grid_timeres_keys\
                              if key not in keys_no_duplicates]
    grid_timeres_keys = keys_no_duplicates
    
    # build content dictionary containing all varlists being queried
    content_dict = {}
    for curr_key in grid_timeres_keys:
        if curr_key not in varlist_url_ids.keys():
            continue
        curr_content_dict,varname_pos,descrip_pos = varlist_file_parse(\
                                    varlist_url_root+varlist_url_ids[curr_key])
        content_dict = {**content_dict,**curr_content_dict}
    
    # find matches to query in the varlists
    shortnames_match = query_shortnames(query,content_dict,\
                                        grid,test_grids,\
                                        time_res,test_timeres)
    
    # print to screen the query match results
    print_varlist_query_results(query,shortnames_match,\
                                content_dict,varname_pos,descrip_pos)
    
    if len(shortnames_match) == 0:
        raise ValueError('No valid matches to query found; please try again.')
    if len(shortnames_match) == 1:
        option_proceed = input('Proceed with option 1? [y/n]: ')
        if option_proceed.casefold() == 'y':
            option_num = 1
        else:
            raise ValueError('No valid matches to query found; please try again.')
    else:
        option_num = input('Please select option [1-'+str(len(shortnames_match))+']: ')
    shortname_match = shortnames_match[int(option_num)-1]
    print('Using dataset with ShortName: '+shortname_match)
    
    return shortname_match
