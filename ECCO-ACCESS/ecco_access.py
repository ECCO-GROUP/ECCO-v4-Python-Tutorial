### This function allows users to query ECCO variables and datasets, and then gain access via direct download, or opening files remotely on S3


def ecco_podaac_access(query,version='v4r4',grid=None,time_res='all',\
                StartDate=None,EndDate=None,\
                mode='download_ifspace',download_root_dir=None,**kwargs):
    """
    
    This function queries and accesses ECCO datasets from PO.DAAC. The core query and download functions are adapted from Jupyter notebooks 
    created by Jack McNelis and Ian Fenty 
    (https://github.com/ECCO-GROUP/ECCO-ACCESS/blob/master/PODAAC/Downloading_ECCO_datasets_from_PODAAC/README.md)
    and modified by Andrew Delman (https://ecco-v4-python-tutorial.readthedocs.io).
    
    Parameters
    ----------    
    query: str, list, or dict, defines datasets or variables to access.
           If query is str, it specifies either a dataset ShortName (which is 
           assumed if the string begins with 'ECCO_'), or a text string that 
           can be used to search the ShortNames, variable names, and descriptions.
           A query may also be a list of multiple ShortNames and/or text searches, 
           or a dict that contains grid,time_res specifiers as keys and ShortNames 
           or text searches as values, e.g.,
           {'native,monthly':['ECCO_L4_SSH_LLC0090GRID_MONTHLY_V4R4',
                              'THETA']}
           will query the native grid monthly SSH datasets, and all native grid 
           monthly datasets with variables or descriptions matching 'THETA'.
    
    version: ('v4r4'), specifies ECCO version to query
    
    grid: ('native','latlon',None), specifies whether to query datasets with output
          on the native grid or the interpolated lat/lon grid.
          The default None will query both types of grids, unless specified 
          otherwise in a query dict (e.g., the example above).
    
    time_res: ('monthly','daily','snapshot','all'), specifies which time resolution 
              to include in query and downloads. 'all' includes all time resolutions, 
              and datasets that have no time dimension, such as the grid parameter 
              and mixing coefficient datasets.

    
    StartDate,EndDate: str, in 'YYYY', 'YYYY-MM', or 'YYYY-MM-DD' format, 
                       define date range [StartDate,EndDate] for download.
                       EndDate is included in the time range (unlike typical Python ranges).
                       ECCOv4r4 date range is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable closed budgets
                       within the specified date range.
    
    mode: str, one of the following:
          'ls' or 'query': Query dataset ShortNames and variable names/
                           descriptions only; no downloads.
          's3_ls' or 's3_query': Query dataset ShortNames and variable names/
                                 descriptions only; return paths on S3.
          'download': Download datasets using NASA Earthdata URLs
          'download_ifspace': Check storage availability before downloading.
                              Download only if storage footprint of downloads 
                              <= max_avail_frac*(available storage)
          'download_subset': Download spatial and temporal subsets of datasets 
                             via Opendap; query help(ecco_podaac_download_subset)
                             to see keyword arguments that can be used in this mode.
          The following modes work within the AWS cloud only:
          's3_open': Access datasets on S3 without downloading.
          's3_get': Download from S3 (to AWS EC2 instance).
          's3_get_ifspace': Check storage availability before downloading; 
                            download if storage footprint 
                            <= max_avail_frac*(available storage).
                            Otherwise data are opened "remotely" from S3 bucket.
          's3_fsspec': Use `fsspec` json files (generated with `kerchunk`) 
                       for expedited loading of datasets.

    download_root_dir: str, defines parent directory to download files to.
                       Files will be downloaded to directory download_root_dir/ShortName/.
                       If not specified, parent directory defaults to '~/Downloads/ECCO_V4r4_PODAAC/'.
    
    Additional keyword arguments*:
    *This is not an exhaustive list, especially for 
    'download_subset' mode; use help(ecco_podaac_download_subset) to display 
    options specific to that mode
    
    max_avail_frac: float, maximum fraction of remaining available disk space to 
                    use in storing ECCO datasets.
                    If storing the datasets exceeds this fraction, an error is returned.
                    Valid range is [0,0.9]. If number provided is outside this range, it is replaced by the closer 
                    endpoint of the range.
    
    n_workers: int, number of workers to use in concurrent downloads. Benefits typically taper off above 5-6.
    
    force_redownload: bool, if True, existing files will be redownloaded and replaced;
                            if False (default), existing files will not be replaced.

    return_granules: bool, if True (default), str or list of queried or 
                           downloaded granules/files (including ones that 
                           were already on disk and not replaced) is returned.
                           if False, the function returns nothing.

    Returns
    -------
    download_files: str, list, or dict, queried or downloaded file(s) 
                    with either URLs (if in 'query' mode), or paths that can be 
                    passed directly to xarray (open_dataset or open_mfdataset).
                    A str is returned if query finds only one granule/file.
                    A list is returned if query finds multiple granules in the 
                    same dataset.
                    A dict (with ShortNames as keys) is returned if the query 
                    finds granules in multiple datasets.
                    Only returned if return_granules=True (default).
    
    """
    
    pass
    
    
    ## query varlists as needed to obtain shortnames
    
    def shortnames_find(query_list,grid,time_res):
        shortnames_list = []
        for query_item in query_list:
            if 'ECCO_' in query_item:
                shortnames_list.append(query_item)
            else:
        
        return shortnames_list
    
    
    if isinstance(query,str):
        query = [query]
    if isinstance(query,dict):
        shortnames = []
        for gridtime_spec,curr_query in query.items():
            if isinstance(curr_query,str):
                curr_query = [curr_query]
            shortnames += shortnames_find(curr_query,\
                                          grid=curr_grid,\
                                          time_res=curr_time_res)
    else:
        shortnames = shortnames_find(query,grid=grid,time_res=time_res)
    
    
    ## query NASA Earthdata CMR and download granules
          'ls' or 'query': Query dataset ShortNames and variable names/
                           descriptions only; no downloads.
          's3_ls' or 's3_query': Query dataset ShortNames and variable names/
                                 descriptions only; return paths on S3.
          'download': Download datasets using NASA Earthdata URLs
          'download_ifspace': Check storage availability before downloading.
                              Download only if storage footprint of downloads 
                              <= max_avail_frac*(available storage)
          'download_subset': Download spatial and temporal subsets of datasets 
                             via Opendap; query help(ecco_podaac_download_subset)
                             to see keyword arguments that can be used in this mode.
          The following modes work within the AWS cloud only:
          's3_open': Access datasets on S3 without downloading.
          's3_get': Download from S3 (to AWS EC2 instance).
          's3_get_ifspace': Check storage availability before downloading; 
                            download if storage footprint 
                            <= max_avail_frac*(available storage).
                            Otherwise data are opened "remotely" from S3 bucket.
          's3_fsspec': Use `fsspec` json files (generated with `kerchunk`) 
                       for expedited loading of datasets.
    
    possible_mode_list = "['ls','query','s3_ls','s3_query','download',\n"\
                         +"'download_ifspace','download_subset',\n"\
                         +"'s3_open','s3_get','s3_get_ifspace','s3_fsspec']"
    # set some default keyword arguments
    kwargs_dict = {}
    if (('n_workers' not in locals()) and (mode != 'download_subset')):
        kwargs_dict['n_workers'] = 6
    if 'force_redownload' not in locals():
        kwargs_dict['force_redownload'] = False
    
    
    # download or otherwise access granules, depending on mode
    
    if mode in ['download_ifspace','s3_get_ifspace']:
        if 'max_avail_frac' not in locals():
            kwargs_dict['max_avail_frac'] = 0.5
        if mode == 'download_ifspace':
            granule_files = ecco_podaac_download_diskaware(\
                               shortnames,StartDate,EndDate,**kwargs_dict)
        elif mode == 's3_get_ifspace':
            granule_files = ecco_podaac_s3_get_diskaware(\
                               shortnames,StartDate,EndDate,**kwargs_dict)
        else:
            raise ValueError('Invalid mode specified; please specify one of the following:'\
              +'\n'+possible_mode_list)
    else:
        granule_files = {}
        for shortname in shortnames:
            if mode in ['ls','query']:
                urls = ecco_podaac_query(shortname,StartDate,EndDate)
                granule_files[shortname] = urls
            elif mode in ['s3_ls','s3_query']:
                granule_files[shortname] = ecco_podaac_s3_query(\
                                              shortname,StartDate,EndDate)
            elif mode == 'download':
                kwargs_dict['return_downloaded_files'] = True
                granule_files[shortname] = ecco_podaac_download(\
                                              shortname,StartDate,EndDate,\
                                              download_root_dir=download_root_dir,\
                                              **kwargs_dict)
            elif mode == 'download_subset':
                if 'n_workers' not in locals():
                    kwargs_dict['n_workers'] = 4
                kwargs_dict['return_downloaded_files'] = True
                granule_files[shortname] = ecco_podaac_download_subset(\
                                              shortname,StartDate,EndDate,\
                                              **kwargs_dict)
            elif mode == 's3_open':
                granule_files[shortname] = ecco_podaac_s3_open(\
                                              shortname,StartDate,EndDate)
            elif mode == 's3_get':
                kwargs_dict['return_downloaded_files'] = True
                granule_files[shortname] = ecco_podaac_s3_get(\
                                              shortname,StartDate,EndDate,\
                                              download_root_dir=download_root_dir,\
                                              **kwargs_dict)
            elif mode == 's3_fsspec':
                
            else:
                raise ValueError('Invalid mode specified; please specify one of the following:'\
                  +'\n'+possible_mode_list)
    
    
    # return granule/file list
    
    if 'return_granules' not in locals():
        return_granules = True
    if return_granules:
        for shortname in granule_files.keys():
            if len(granule_files[shortname]) == 1:
                # if only 1 file is downloaded, return a string of filename instead of a list
                granule_files = granule_files[0]
        
        return granule_files