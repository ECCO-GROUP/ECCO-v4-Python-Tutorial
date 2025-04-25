### This function allows users to query ECCO variables and datasets, and then gain access via direct download, or opening files remotely on S3

from .ecco_download import ecco_podaac_query
from .ecco_download import ecco_podaac_download
from .ecco_download import ecco_podaac_download_diskaware
from .ecco_download import ecco_podaac_download
from .ecco_download import ecco_podaac_download_subset

from .ecco_s3_retrieve import ecco_podaac_s3_query
from .ecco_s3_retrieve import ecco_podaac_s3_open
from .ecco_s3_retrieve import ecco_podaac_s3_open_fsspec
from .ecco_s3_retrieve import ecco_podaac_s3_get
from .ecco_s3_retrieve import ecco_podaac_s3_get_diskaware

from .ecco_acc_dates import date_adjustment

from .ecco_varlist import ecco_podaac_varlist_query


import requests


def ecco_podaac_access(query,version='v4r4',grid=None,time_res='all',\
                StartDate=None,EndDate=None,snapshot_interval=None,\
                mode='download_ifspace',download_root_dir=None,**kwargs):
    """
    
    This function queries and accesses ECCO datasets from PO.DAAC. The core query and download functions 
    are adapted from Jupyter notebooks created by Jack McNelis and Ian Fenty 
    (https://github.com/ECCO-GROUP/ECCO-ACCESS/blob/master/PODAAC/Downloading_ECCO_datasets_from_PODAAC/README.md)
    and modified by Andrew Delman (https://ecco-v4-python-tutorial.readthedocs.io).
    
    Parameters
    ----------    
    query: str, list, or dict, defines datasets or variables to access.
           If query is str, it specifies either a dataset ShortName (if query 
           matches a NASA Earthdata ShortName), or a text string that can be 
           used to search the ECCO ShortNames, variable names, and descriptions.
           A query may also be a list of multiple ShortNames and/or text searches, 
           or a dict that contains grid,time_res specifiers as keys and ShortNames 
           or text searches as values, e.g.,
           {'native,monthly':['ECCO_L4_SSH_LLC0090GRID_MONTHLY_V4R4',
                              'THETA']}
           will query the native grid monthly SSH datasets, and all native grid 
           monthly datasets with variables or descriptions matching 'THETA'.
    
    version: ('v4r4','v4r5'), specifies ECCO version to query.
             Currently 'v4r5' only works with ['s3_open','s3_get','s3_get_ifspace'] modes, 
             or if the files are already stored in download_root_dir/ShortName/.
             Otherwise an error is returned.
             'v4r5' only has grid='native' and time_res='monthly' data files available.
    
    grid: ('native','latlon',None), specifies whether to query datasets with output
          on the native grid or the interpolated lat/lon grid.
          The default None will query both types of grids (and datasets with no spatial
          dimension), unless specified otherwise in a query dict (e.g., the example above).
    
    time_res: ('monthly','daily','snapshot','all'), specifies which time resolution 
              to include in query and downloads. 'all' includes all time resolutions, 
              and datasets that have no time dimension, such as the grid parameter 
              and mixing coefficient datasets.

    
    StartDate,EndDate: str, in 'YYYY', 'YYYY-MM', or 'YYYY-MM-DD' format, 
                       define date range [StartDate,EndDate] for download.
                       EndDate is included in the time range (unlike typical Python ranges).
                       Full ECCOv4r4 date range (default) is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable 
                       closed budgets within the specified date range.
    
    snapshot_interval: ('monthly', 'daily', or None), if snapshot datasets are included in ShortNames, 
                       this determines whether snapshots are included for only the beginning/end 
                       of each month ('monthly'), or for every day ('daily').
                       If None or not specified, defaults to 'daily' if any daily mean ShortNames 
                       are included and 'monthly' otherwise.
    
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
                             via Opendap; query help(ecco_access.ecco_podaac_download_subset)
                             to see keyword arguments that can be used in this mode.
          The following modes work within the AWS cloud only:
          's3_open': Access datasets on S3 without downloading.
          's3_open_fsspec': Use json files (generated with `fsspec` and `kerchunk`) 
                            for expedited opening of datasets.
          's3_get': Download from S3 (to AWS EC2 instance).
          's3_get_ifspace': Check storage availability before downloading; 
                            download if storage footprint 
                            <= max_avail_frac*(available storage).
                            Otherwise data are opened "remotely" from S3 bucket.

    download_root_dir: str, defines parent directory to download files to.
                       Files will be downloaded to directory download_root_dir/ShortName/.
                       If not specified, parent directory defaults to '~/Downloads/ECCO_V4r4_PODAAC/',
                       or '~/Downloads/ECCO_V4r5_PODAAC/' if version == 'v4r5'.
    
    Additional keyword arguments*:
    *This is not an exhaustive list, especially for 
    'download_subset' mode; use help(ecco_access.ecco_podaac_download_subset) to display 
    options specific to that mode
    
    max_avail_frac: float, maximum fraction of remaining available disk space to 
                    use in storing ECCO datasets.
                    If storing the datasets exceeds this fraction, an error is returned.
                    Valid range is [0,0.9]. If number provided is outside this range, it is replaced 
                    by the closer endpoint of the range.
    
    jsons_root_dir: str, for s3_open_fsspec mode only, the root/parent directory where the 
                    fsspec/kerchunk-generated jsons are found.
                    jsons are generated using the steps described here:
                    https://medium.com/pangeo/fake-it-until-you-make-it-reading-goes-netcdf4-data-on-aws-s3
                    as-zarr-for-rapid-data-access-61e33f8fe685
                    and stored as {jsons_root_dir}/MZZ_{GRIDTYPE}_{TIME_RES}/{SHORTNAME}.json.
                    For v4r4, GRIDTYPE is '05DEG' or 'LLC0090GRID'.
                    TIME_RES is one of: ('MONTHLY','DAILY','SNAPSHOT','GEOMETRY','MIXING_COEFFS').
    
    n_workers: int, number of workers to use in concurrent downloads. Benefits typically taper off above 5-6.
    
    force_redownload: bool, if True, existing files will be redownloaded and replaced;
                            if False (default), existing files will not be replaced.
    
    show_noredownload_msg: bool, if True (default), and force_redownload=False, 
                               display message for each file that is already 
                               downloaded (and therefore not re-downloaded); 
                               if False, these messages are not shown.
    
    prompt_request_payer: bool, if True (default), user is prompted to approve 
                                (by entering "y" or "Y") any access to a 
                                requester pays bucket, otherwise request is canceled; 
                                if False, data access proceeds without prompting.
    
    return_granules: bool, if True (default), str or list of queried or 
                           downloaded granules/files (including ones that 
                           were already on disk and not replaced) is returned.
                           if False, the function returns nothing.

    Returns
    -------
    granule_files: dict with ShortNames as keys; values are URLs or S3 paths
                   (if in 'query' mode), or paths of files that can be 
                   passed directly to xarray (open_dataset or open_mfdataset).
                   Values are of type str if query finds only one granule/file
                   for that ShortName; of type list if query finds 
                   multiple granules in the same dataset; 
                   or of type fsspec.mapping.FSMap if mode = 's3_open_fsspec'.
                   Only returned if return_granules=True (default).
    
    """
    
    
    pass
    
    
    ## query varlists as needed to obtain shortnames
    
    def shortnames_find(query_list,version,grid,time_res):
        shortnames_list = []
        if version == 'v4r5':
            if ((grid == 'native') or (grid is None)):
                if ((time_res == 'monthly') or (time_res == 'all')):
                    import s3fs
                    from os.path import split
                    s3 = s3fs.S3FileSystem(anon=False,\
                                           requester_pays=True)
                    s3_datasets_list = [split(dataset_path)[-1]\
                                          for dataset_path in \
                                          s3.ls("s3://ecco-model-granules/netcdf/V4r5/native/mon_mean/")]
                else:
                    raise ValueError("'"+time_res+"' time res can not currently be accessed for v4r5.\n"\
                                     +"ecco_access can currently access only 'monthly' time_res v4r5 netCDF files.")
            else:
                raise ValueError("'"+grid+"' grid can not currently be accessed for v4r5.\n"\
                                 +"ecco_access can currently access only 'native' grid v4r5 netCDF files.")
                
        for query_item in query_list:
            if version == 'v4r5':
                # see if the query is an existing dataset ID
                # if not, then do a text search of the ECCO variable lists
                if query_item in s3_datasets_list:
                    shortnames_list.append(query_item)
                else:
                    shortname_match = ecco_podaac_varlist_query(query_item,version,grid,time_res)
                    shortnames_list.append(shortname_match)
            else:    
                # see if the query is an existing NASA Earthdata ShortName
                # if not, then do a text search of the ECCO variable lists
                response = requests.get(url="https://cmr.earthdata.nasa.gov/search/collections.json", 
                                        params={'ShortName':query_item})
                if len(response.json()['feed']['entry']) > 0:
                    shortnames_list.append(query_item)
                else:
                    shortname_match = ecco_podaac_varlist_query(query_item,version,grid,time_res)
                    shortnames_list.append(shortname_match)
        
        return shortnames_list
    
    
    if isinstance(query,str):
        query = [query]
    if isinstance(query,dict):
        shortnames = []
        for gridtime_spec,curr_query in query.items():
            try:
                curr_grid,curr_time_res = gridtime_spec.split(',')
            except:
                raise ValueError("Keys of dict 'query' must be of the form grid,time_res\n'\
                                 +'with 1 comma in the middle")
            if isinstance(curr_query,str):
                curr_query = [curr_query]
            shortnames += shortnames_find(curr_query,\
                                          version,\
                                          grid=curr_grid,\
                                          time_res=curr_time_res)
    else:
        shortnames = shortnames_find(query,version,grid=grid,time_res=time_res)
    
    
    ## query NASA Earthdata CMR and download granules
    
    possible_mode_list = "['ls','query','s3_ls','s3_query','download',\n"\
                         +"'download_ifspace','download_subset',\n"\
                         +"'s3_open','s3_get','s3_get_ifspace','s3_open_fsspec']"
    
    # set some default keyword arguments
    if (('n_workers' not in kwargs.keys()) and (mode != 'download_subset')):
        kwargs['n_workers'] = 6
    if 'force_redownload' not in kwargs.keys():
        kwargs['force_redownload'] = False

    # remove unneeded keyword arguments
    if mode == 's3_open_fsspec':
        for kwarg in list(kwargs.keys()):
            if kwarg != 'jsons_root_dir':
                del kwargs[kwarg]
    elif mode == 's3_open':
        for kwarg in list(kwargs.keys()):
            if kwarg in ['n_workers','force_redownload','show_noredownload_msg']:
                del kwargs[kwarg]
    else:
        if 'jsons_root_dir' in kwargs.keys():
            del kwargs['jsons_root_dir']
    
    
    # download or otherwise access granules, depending on mode
    
    if mode in ['download_ifspace','s3_get_ifspace']:
        if 'max_avail_frac' not in kwargs.keys():
            kwargs['max_avail_frac'] = 0.5
        if mode == 'download_ifspace':
            granule_files = ecco_podaac_download_diskaware(\
                               shortnames,StartDate,EndDate,version,snapshot_interval,\
                               download_root_dir=download_root_dir,**kwargs)
        elif mode == 's3_get_ifspace':
            granule_files = ecco_podaac_s3_get_diskaware(\
                               shortnames,StartDate,EndDate,version,snapshot_interval,\
                               download_root_dir=download_root_dir,**kwargs)
        else:
            raise ValueError('Invalid mode specified; please specify one of the following:'\
              +'\n'+possible_mode_list)
    else:
        if 'max_avail_frac' in kwargs.keys():
            del kwargs['max_avail_frac']
        granule_files = {}
        
        # determine value of snapshot_interval if None or not specified
        if snapshot_interval == None:
            snapshot_interval = 'monthly'
            for curr_shortname in shortnames:
                if 'DAILY' in curr_shortname:
                    snapshot_interval = 'daily'
                    break
        
        for shortname in shortnames:
            
            if mode in ['ls','query']:
                urls,sizes = ecco_podaac_query(shortname,StartDate,EndDate,version,snapshot_interval)
                granule_files[shortname] = urls
            elif mode in ['s3_ls','s3_query']:
                granule_files[shortname] = ecco_podaac_s3_query(\
                                              shortname,StartDate,EndDate,version,snapshot_interval)
            elif mode == 'download':
                kwargs['return_downloaded_files'] = True
                granule_files[shortname] = ecco_podaac_download(\
                                              shortname,StartDate,EndDate,version,snapshot_interval,\
                                              download_root_dir=download_root_dir,\
                                              **kwargs)
            elif mode == 'download_subset':
                if 'n_workers' not in kwargs.keys():
                    kwargs['n_workers'] = 4
                kwargs['return_downloaded_files'] = True
                granule_files[shortname] = ecco_podaac_download_subset(\
                                              shortname,StartDate,EndDate,snapshot_interval,\
                                              download_root_dir=download_root_dir,\
                                              **kwargs)
            elif mode == 's3_open':
                granule_files[shortname] = ecco_podaac_s3_open(\
                                              shortname,StartDate,EndDate,version,snapshot_interval,\
                                              **kwargs)
            elif mode == 's3_open_fsspec':
                # granule_files will consist of mapper objects rather than URL/path or file lists
                granule_files[shortname] = ecco_podaac_s3_open_fsspec(\
                                              shortname,version,**kwargs)
            elif mode == 's3_get':
                kwargs['return_downloaded_files'] = True
                granule_files[shortname] = ecco_podaac_s3_get(\
                                              shortname,StartDate,EndDate,version,snapshot_interval,\
                                              download_root_dir=download_root_dir,\
                                              **kwargs)
            else:
                raise ValueError('Invalid mode specified; please specify one of the following:'\
                  +'\n'+possible_mode_list)
    
    
    # return granule/file list
    
    if 'return_granules' not in kwargs.keys():
        return_granules = True
    if return_granules:
        for shortname in granule_files.keys():
            if isinstance(granule_files[shortname],list):
                if ((len(granule_files[shortname]) == 1) and (mode != 's3_open_fsspec')):
                    # if only 1 file is downloaded, return a string of filename instead of a list
                    granule_files[shortname] = granule_files[shortname][0]
        
        return granule_files



###================================================================================================================



def ecco_podaac_to_xrdataset(query,version='v4r4',grid=None,time_res='all',\
                             StartDate=None,EndDate=None,snapshot_interval=None,\
                             mode='download_ifspace',download_root_dir=None,**kwargs):
    """
    
    This function queries and accesses ECCO datasets from PO.DAAC. The core query and download functions 
    are adapted from Jupyter notebooks created by Jack McNelis and Ian Fenty 
    (https://github.com/ECCO-GROUP/ECCO-ACCESS/blob/master/PODAAC/Downloading_ECCO_datasets_from_PODAAC/README.md)
    and modified by Andrew Delman (https://ecco-v4-python-tutorial.readthedocs.io).
    It is similar to ecco_podaac_access, except instead of a list of URLs or files, 
    an xarray Dataset with all of the queried ECCO datasets is returned.

    Parameters
    ----------    
    query: str, list, or dict, defines datasets or variables to access.
           If query is str, it specifies either a dataset ShortName (if query 
           matches a NASA Earthdata ShortName), or a text string that can be 
           used to search the ECCO ShortNames, variable names, and descriptions.
           A query may also be a list of multiple ShortNames and/or text searches, 
           or a dict that contains grid,time_res specifiers as keys and ShortNames 
           or text searches as values, e.g.,
           {'native,monthly':['ECCO_L4_SSH_LLC0090GRID_MONTHLY_V4R4',
                              'THETA']}
           will query the native grid monthly SSH datasets, and all native grid 
           monthly datasets with variables or descriptions matching 'THETA'.
    
    version: ('v4r4','v4r5'), specifies ECCO version to query.
             Currently 'v4r5' only works with ['s3_open','s3_get','s3_get_ifspace'] modes, 
             or if the files are already stored in download_root_dir/ShortName/.
             Otherwise an error is returned.
             'v4r5' only has grid='native' and time_res='monthly' data files available.
    
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
                       Full ECCOv4r4 date range (default) is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable 
                       closed budgets within the specified date range.
    
    snapshot_interval: ('monthly', 'daily', or None), if snapshot datasets are included in ShortNames, 
                       this determines whether snapshots are included for only the beginning/end 
                       of each month ('monthly'), or for every day ('daily').
                       If None or not specified, defaults to 'daily' if any daily mean ShortNames 
                       are included and 'monthly' otherwise.
    
    mode: str, one of the following:
          'download': Download datasets using NASA Earthdata URLs
          'download_ifspace': Check storage availability before downloading.
                              Download only if storage footprint of downloads 
                              <= max_avail_frac*(available storage)
          'download_subset': Download spatial and temporal subsets of datasets 
                             via Opendap; query help(ecco_access.ecco_podaac_download_subset)
                             to see keyword arguments that can be used in this mode.
          The following modes work within the AWS cloud only:
          's3_open': Access datasets on S3 without downloading.
          's3_open_fsspec': Use json files (generated with `fsspec` and `kerchunk`) 
                            for expedited opening of datasets.
          's3_get': Download from S3 (to AWS EC2 instance).
          's3_get_ifspace': Check storage availability before downloading; 
                            download if storage footprint 
                            <= max_avail_frac*(available storage).
                            Otherwise data are opened "remotely" from S3 bucket.

    download_root_dir: str, defines parent directory to download files to.
                       Files will be downloaded to directory download_root_dir/ShortName/.
                       If not specified, parent directory defaults to '~/Downloads/ECCO_V4r4_PODAAC/',
                       or '~/Downloads/ECCO_V4r5_PODAAC/' if version == 'v4r5'.
    
    Additional keyword arguments*:
    *This is not an exhaustive list, especially for 
    'download_subset' mode; use help(ecco_access.ecco_podaac_download_subset) to display 
    options specific to that mode
    
    max_avail_frac: float, maximum fraction of remaining available disk space to 
                    use in storing ECCO datasets.
                    If storing the datasets exceeds this fraction, an error is returned.
                    Valid range is [0,0.9]. If number provided is outside this range, it is replaced 
                    by the closer endpoint of the range.
    
    jsons_root_dir: str, for s3_open_fsspec mode only, the root/parent directory where the 
                    fsspec/kerchunk-generated jsons are found.
                    jsons are generated using the steps described here:
                    https://medium.com/pangeo/fake-it-until-you-make-it-reading-goes-netcdf4-data-on-aws-s3
                    as-zarr-for-rapid-data-access-61e33f8fe685
                    and stored as {jsons_root_dir}/MZZ_{GRIDTYPE}_{TIME_RES}/{SHORTNAME}.json.
                    For v4r4, GRIDTYPE is '05DEG' or 'LLC0090GRID'.
                    TIME_RES is one of: ('MONTHLY','DAILY','SNAPSHOT','GEOMETRY','MIXING_COEFFS').
    
    n_workers: int, number of workers to use in concurrent downloads. Benefits typically taper off above 5-6.
    
    force_redownload: bool, if True, existing files will be redownloaded and replaced;
                            if False (default), existing files will not be replaced.
    
    show_noredownload_msg: bool, if True (default), and force_redownload=False, 
                               display message for each file that is already 
                               downloaded (and therefore not re-downloaded); 
                               if False, these messages are not shown.
    
    prompt_request_payer: bool, if True (default), user is prompted to approve 
                                (by entering "y" or "Y") any access to a 
                                requester pays bucket, otherwise request is canceled; 
                                if False, data access proceeds without prompting.

    Returns
    -------
    ds_out: xarray Dataset or dict of xarray Datasets (with ShortNames as keys), 
            containing all of the accessed datasets.
            This function does not work with the query modes: 'ls','query','s3_ls','s3_query'.
    """
    
    pass
    
    
    import numpy as np
    import xarray as xr
    

    # raise error if mode is ls/query only
    if mode in ['ls','query','s3_ls','s3_query']:
        raise ValueError("ecco_podaac_to_xrdataset does not work with 'ls'/'query' modes. \n"\
                         +"Please use ecco_podaac_access with these modes.")
        
        return -1
    
    # submit access query (and download if needed)
    access_output = ecco_podaac_access(query,version,grid,time_res,\
                                       StartDate,EndDate,snapshot_interval,\
                                       mode,download_root_dir,**kwargs)
    
    # determine value of snapshot_interval if None or not specified
    if snapshot_interval == None:
        snapshot_interval = 'monthly'
        for curr_shortname in access_output.keys():
            if 'DAILY' in curr_shortname:
                snapshot_interval = 'daily'
                break
    
    # open xarray datasets
    ds_out = {}
    for shortname,access_out in access_output.items():
        if mode == 's3_open_fsspec':
            curr_ds = xr.open_dataset(access_out,engine='zarr',chunks='auto',consolidated=False)
            if 'time' in curr_ds.dims:
                # isolate time range specified
                startdate,enddate = date_adjustment(shortname,\
                                                    StartDate,EndDate,CMR_query=False)
                time_values = curr_ds.time.values.astype('datetime64[D]')
                in_time_range = np.logical_and(time_values >= startdate,\
                                               time_values <= enddate).nonzero()[0]
                curr_ds = curr_ds.isel(time=in_time_range)
                if (('SNAPSHOT' in shortname) and (snapshot_interval == 'monthly')):
                    month_bounds_list = np.arange(np.datetime64('1992-01','M'),\
                                                  np.datetime64('2040-01','M'),\
                                                  np.timedelta64(1,'M'))\
                                                  .astype('datetime64[D]')
                    time_values = curr_ds.time.values.astype('datetime64[D]')
                    time_subind = list(np.arange(0,len(time_values)).astype('int64'))
                    for count,time_val in enumerate(time_values):
                        if time_val not in month_bounds_list:
                            time_subind.remove(count)
                    curr_ds = curr_ds.isel(time=time_subind)
        else:
            if not isinstance(access_out,list):
                access_out = [access_out]
            curr_ds = xr.open_mfdataset(access_out,\
                                        compat='override',data_vars='minimal',coords='minimal',\
                                        parallel=True)
        ds_out[shortname] = curr_ds
    
    # if only one ShortName is involved, then extract dataset from dictionary
    if len(ds_out) == 1:
        ds_out = list(ds_out.values())[0]
    
    return ds_out
