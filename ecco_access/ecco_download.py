### This module contains routines to download ECCO datasets using Python requests


from .ecco_acc_dates import date_adjustment

## Initalize Python libraries
import numpy as np
import pandas as pd
import requests
import shutil
import time as time

# for concurrent simulatenous downloads
from concurrent.futures import ThreadPoolExecutor
from getpass import getpass
from http.cookiejar import CookieJar
from io import StringIO
from itertools import repeat
from pathlib import Path
from platform import system
from netrc import netrc
from os.path import basename, isfile, isdir, join, expanduser
import sys
# progress bar
from tqdm import tqdm
# library to download files
from urllib import request    


def ecco_podaac_query(ShortName,StartDate,EndDate,version,snapshot_interval='monthly'):
    
    """
    
    This routine queries NASA Earthdata to find ECCO datasets hosted by PO.DAAC. It is adapted from the Jupyter notebooks 
    created by Jack McNelis and Ian Fenty 
    (https://github.com/ECCO-GROUP/ECCO-ACCESS/blob/master/PODAAC/Downloading_ECCO_datasets_from_PODAAC/README.md)
    and modified by Andrew Delman (https://ecco-v4-python-tutorial.readthedocs.io).
    
    Parameters
    ----------    
    ShortName: str, the ShortName that identifies the dataset on PO.DAAC.
    
    StartDate,EndDate: str, in 'YYYY', 'YYYY-MM', or 'YYYY-MM-DD' format, 
                       define date range [StartDate,EndDate] for download.
                       EndDate is included in the time range (unlike typical Python ranges).
                       ECCOv4r4 date range is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable closed budgets
                       within the specified date range.

    version: ('v4r4', 'v4r5'), which version of ECCO to query.
    
    snapshot_interval: ('monthly', 'daily'), if the dataset corresponding to ShortName is a snapshot, 
                       determines whether snapshots are included for only the beginning/end of each month 
                       ('monthly'), or for every day ('daily'). Defaults to 'monthly'.


    Returns
    -------
    urls,sizes: list, file URLs that match the query and sizes; this list can be used to download the files
    
    """
    
    pass
    
    
    #=====================================================
    
    ### Define Helper Subroutines
    
    ### Helper subroutine to log into NASA EarthData
    
    # not pretty but it works
    def setup_earthdata_login_auth(url: str='urs.earthdata.nasa.gov'):
        # look for the netrc file and use the login/password
        try:
            username, _, password = netrc(file=_netrc).authenticators(url)
    
        # if the file is not found, prompt the user for the login/password
        except (FileNotFoundError, TypeError):
            print('Please provide Earthdata Login credentials for access.')
            username, password = input('Username: '), getpass('Password: ')
        
        manager = request.HTTPPasswordMgrWithDefaultRealm()
        manager.add_password(None, url, username, password)
        auth = request.HTTPBasicAuthHandler(manager)
        jar = CookieJar()
        processor = request.HTTPCookieProcessor(jar)
        opener = request.build_opener(auth, processor)
        request.install_opener(opener)
    
    ### Helper subroutines to make the API calls to search CMR and parse response
    def set_params(params: dict):
#         params.update({'scroll': "true", 'page_size': 2000})
        params.update({'page_size': 2000})
        return {par: val for par, val in params.items() if val is not None}
    
    def get_results(params: dict, headers: dict=None):
        response = requests.get(url="https://cmr.earthdata.nasa.gov/search/granules.json", 
                                params=set_params(params),
                                headers=headers).json()
        return response
    
    
#     def get_granules(params: dict):
#         time_start = np.array([]).astype('datetime64[ns]')
#         response, headers = get_results(params=params)
# #         scroll = headers['CMR-Scroll-Id']
#         hits = int(headers['CMR-Hits'])
#         if hits==0:
#             raise Exception("No granules matched your input parameters.")
#         df = pd.read_csv(StringIO(response.text)) 
#         while hits > df.index.size:
# #             response, _ = get_results(params=params, headers={'CMR-Scroll-Id': scroll})
#             response, _ = get_results(params=params)
#             data = pd.read_csv(StringIO(response.text))
#             df = pd.concat([df, data])
#         return df
    
    def get_granules(params: dict, ShortName: str, SingleDay_flag: bool):
        time_start = np.array([]).astype('datetime64[ns]')
        urls = []
        sizes = []
        completed_query = False
        while completed_query == False:
            response = get_results(params=params)
            if 'feed' in response.keys():
                for curr_entry in response['feed']['entry']:
                    time_start = np.append(time_start,np.datetime64(curr_entry['time_start'],'ns'))
                    sizes.append(curr_entry['granule_size'])
                    for curr_link in curr_entry['links']:
                        if ".nc" in curr_link['title'][-3:]:
                            urls.append(curr_link['href'])
                            break
            elif 'errors' in response.keys():
                raise Exception(response['errors'][0])
            
            if len(response['feed']['entry']) < 2000:
                completed_query = True
            else:
                # do another CMR search since previous search hit the allowed maximum
                # number of entries (2000)
                params['temporal'] = str(np.datetime64(response['feed']['entry'][-1]['time_end'],'D')\
                                         + np.timedelta64(1,'D'))+params['temporal'][10:]

        # reduce granule list to single day if only one day in requested range
        if (('MONTHLY' in ShortName) or ('DAILY' in ShortName)):
            if ((SingleDay_flag == True) and (len(urls) > 1)):
                day_index = np.argmin(np.abs(time_start - np.datetime64(StartDate,'D')))
                urls = urls[day_index:(day_index+1)]

        return urls,sizes
    
    
    
    #=====================================================
    
    
    
    # # set default StartDate or EndDate if not previously provided
    if StartDate == None:
        StartDate = '1992-01-01'
    if EndDate == None:
        EndDate = '2099-12-31'
    
    # # Adjust StartDate and EndDate to CMR query values
    StartDate,EndDate,SingleDay_flag = date_adjustment(ShortName,\
                                         StartDate,EndDate,CMR_query=True)    
    
    
    ## Log into Earthdata using your username and password
    
    # Predict the path of the netrc file depending on os/platform type.
    _netrc = join(expanduser('~'), "_netrc" if system()=="Windows" else ".netrc")
    
    # actually log in with this command:
    setup_earthdata_login_auth()
    
    
    # Query the NASA Common Metadata Repository to find the URL of every granule associated with the desired ECCO Dataset and date range of interest.
    
    # create a Python dictionary with our search criteria:  `ShortName` and `temporal`
    input_search_params = {'ShortName': ShortName,
                           'temporal': ",".join([StartDate, EndDate])}
    
    
    if version == 'v4r5':
        urls = glob.glob(join(download_root_dir,ShortName,'*.nc'))
        sizes = [0]*len(urls)
    else:
        ### Query CMR for the desired ECCO Dataset
        
        # grans means 'granules', PO.DAAC's term for individual files in a dataset
        urls,gran_sizes = get_granules(input_search_params,ShortName,SingleDay_flag)
        
    #     ## Prepare results of query
    #     
    #     # reduce granule list to single day if only one day in requested range
    #     if (('MONTHLY' in ShortName) or ('DAILY' in ShortName)):
    #         if ((SingleDay_flag == True) and (len(grans['Granule UR']) > 1)):
    #             day_index = np.argmin(np.abs(np.asarray(grans['Start Time'])\
    #                        .astype('datetime64[ns]') - np.datetime64(StartDate,'D')))
    #             grans = grans[day_index:(day_index+1)]
    #     
    #     # convert the rows of the 'Online Access URLS' column to a Python list
    #     urls = grans['Online Access URLs'].tolist()
        
        # estimate granule sizes where this info is missing from CMR
        sizes = (2**20)*np.asarray(gran_sizes).astype('float64')
        sizes = np.where(sizes > (2**10),sizes,np.nan) 
        if np.sum(~np.isnan(sizes)) >= 1:
            sizes = np.where(~np.isnan(sizes),sizes,np.nanmean(sizes))
        else:
            input_search_params['temporal'] = ['1992-01-01','2017-12-31']
            _,gran_sizes_all = get_granules(input_search_params)
            sizes_all = (2**20)*np.asarray(grans_all['Size']).astype('float64')
            sizes_all = np.where(sizes_all > (2**10),sizes_all,np.nan) 
            sizes = np.where(~np.isnan(sizes),sizes,np.nanmean(sizes_all))
        sizes = list(sizes)
    #     urls = grans['Online Access URLs'].tolist()
        
        # for snapshot datasets with monthly snapshot_interval, only include snapshots at beginning/end of months
        if 'SNAPSHOT' in ShortName:
            if snapshot_interval == 'monthly':
                import re
                url_sizes_dict = {url:size for url,size in zip(urls,sizes)}
                for url,size in zip(urls,sizes):
                    snapshot_date = re.findall("_[0-9]{4}-[0-9]{2}-[0-9]{2}",url)[0][1:]
                    if snapshot_date[8:] != '01':
                        del url_sizes_dict[url]
                urls = list(url_sizes_dict.keys())
                sizes = list(url_sizes_dict.values())
    
    return urls,sizes


###================================================================================================================


def download_files_wrapper(urls, download_dir, version, n_workers, force_redownload, show_noredownload_msg):
    """Wrapper for downloading functions"""
    
    pass
    
    #=====================================================
    
    ### Define Helper Subroutines
    
    ### Helper subroutine to gracefully download single files and avoids re-downloading if file already exists.
    # To force redownload of the file, pass **True** to the boolean argument *force* (default **False**)\n,
    def download_file(url: str, output_dir: str, version: str, force: bool=False, show_noredownload_msg: bool=True):
        """url (str): the HTTPS url from which the file will download
        output_dir (str): the local path into which the file will download
        version (str): the version of ECCO to be downloaded
        force (bool): download even if the file exists locally already
        show_noredownload_msg (bool): show no re-download messages (vs. not showing messages)
        """
        if not isdir(output_dir):
            raise Exception(f"Output directory doesnt exist! ({output_dir})")
        
        target_file = join(output_dir, basename(url))
        
        # if the file has already been downloaded, skip
        if isfile(target_file) and force is False:
            if show_noredownload_msg:
                print(f'\n{basename(url)} already exists, and force=False, not re-downloading')
            return target_file,0
        elif version == 'v4r5':
            raise Exception('Capability to download v4r5 files via HTTPS has not been added yet.\n'\
                            +'Please download files from ECCO Drive or use v4r4 files.')
        
        with requests.get(url) as r:
            if not r.status_code // 100 == 2: 
                raise Exception(r.text)
                return target_file,0
            else:
                with open(target_file, 'wb') as f:
                    total_size_in_bytes = int(r.headers.get('content-length', 0))
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk:
                            f.write(chunk)
    
        return target_file,total_size_in_bytes
    
    
    ### Helper subroutine to download all urls in the list `dls`
    def download_files_concurrently(dls, download_dir, version, n_workers, force=False, show_noredownload_msg=True):
        start_time = time.time()
    
        # use thread pool for concurrent downloads
        with ThreadPoolExecutor(max_workers=n_workers) as executor:
    
            # tqdm makes a cool progress bar
            results = list(tqdm(executor.map(download_file, dls, repeat(download_dir),\
                                repeat(version), repeat(force), repeat(show_noredownload_msg)),\
                                total=len(dls), desc='DL Progress',\
                                ascii=True, ncols=75, file=sys.stdout))
            
            # add up the total downloaded file sizes
            total_download_size_in_bytes = 0
            for result in results:
                total_download_size_in_bytes += result[-1]
            # calculate total time spent in the download
            total_time_download = time.time() - start_time
            
            print('\n=====================================')
            print(f'total downloaded: {np.round(total_download_size_in_bytes/1e6,2)} Mb')
            print(f'avg download speed: {np.round(total_download_size_in_bytes/1e6/total_time_download,2)} Mb/s')
            print('Time spent = ' + str(total_time_download) + ' seconds')
            print('\n')
            
            # return list of downloaded files
            downloaded_files = []
            for result in results:
                downloaded_files.append(result[0])

            return downloaded_files
    
    
    #=====================================================
    
    try:
        # Attempt concurrent downloads, but if error arises switch to sequential downloads
        ### Method 1: Concurrent downloads        
       
        # Force redownload (or not) depending on value of force_redownload
        downloaded_files = download_files_concurrently(urls, download_dir, version, n_workers, force_redownload, show_noredownload_msg)
        
    except:
        ### Method 2: Sequential Downloads
        
        start_time = time.time()
        
        # Download each URL sequentially in a for loop.

        downloaded_files = []
        total_download_size_in_bytes = 0
        
        # loop through all urls
        for u in urls:
            u_name = u.split('/')[-1]
            print(f'downloading {u_name}')
            result = download_file(url=u, output_dir=download_dir, force=force_redownload, show_noredownload_msg=show_noredownload_msg)
            downloaded_files.append(result[0])
            total_download_size_in_bytes += result[-1]
        
        # calculate total time spent in the download
        total_time_download = time.time() - start_time
        
        print('\n=====================================')
        print(f'total downloaded: {np.round(total_download_size_in_bytes/1e6,2)} Mb')
        print(f'avg download speed: {np.round(total_download_size_in_bytes/1e6/total_time_download,2)} Mb/s')
        print('Time spent = ' + str(total_time_download) + ' seconds')
        print('\n')
    
    return downloaded_files



###================================================================================================================


def ecco_podaac_download(ShortName,StartDate,EndDate,version,snapshot_interval='monthly',download_root_dir=None,\
                         n_workers=6,force_redownload=False,show_noredownload_msg=True,\
                         return_downloaded_files=False):
    """
    
    This routine downloads ECCO datasets from PO.DAAC. It is adapted from the Jupyter notebooks 
    created by Jack McNelis and Ian Fenty 
    (https://github.com/ECCO-GROUP/ECCO-ACCESS/blob/master/PODAAC/Downloading_ECCO_datasets_from_PODAAC/README.md)
    and modified by Andrew Delman (https://ecco-v4-python-tutorial.readthedocs.io).
    
    Parameters
    ----------    
    ShortName: str, the ShortName that identifies the dataset on PO.DAAC.
    
    StartDate,EndDate: str, in 'YYYY', 'YYYY-MM', or 'YYYY-MM-DD' format, 
                       define date range [StartDate,EndDate] for download.
                       EndDate is included in the time range (unlike typical Python ranges).
                       ECCOv4r4 date range is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable closed budgets
                       within the specified date range.
    
    version: ('v4r4', 'v4r5'), the version of ECCO to download.

    snapshot_interval: ('monthly', 'daily'), if the dataset corresponding to ShortName is a snapshot, 
                       determines whether snapshots are included for only the beginning/end of each month 
                       ('monthly'), or for every day ('daily'). Defaults to 'monthly'.

    download_root_dir: str, defines parent directory to download files to.
                       Files will be downloaded to directory download_root_dir/ShortName/.
                       If not specified, parent directory defaults to '~/Downloads/ECCO_V4r4_PODAAC/',
                       or '~/Downloads/ECCO_V4r5_PODAAC/' if version == 'v4r5'.
    
    n_workers: int, number of workers to use in concurrent downloads. Benefits typically taper off above 5-6.
    
    force_redownload: bool, if True, existing files will be redownloaded and replaced;
                            if False (default), existing files will not be replaced.
    
    show_noredownload_msg: bool, if True (default), and force_redownload=False, 
                               display message for each file that is already 
                               downloaded (and therefore not re-downloaded); 
                               if False, these messages are not shown.

    return_downloaded_files: bool, if True, string or list of downloaded file(s) (including files that were already on disk
                                   and not replaced) is returned.
                                   if False (default), the function returns nothing.

    Returns
    -------
    downloaded_files: str or list, downloaded file(s) with local path that can be passed 
                      directly to xarray (open_dataset or open_mfdataset).
                      Only returned if return_downloaded_files=True.
    
    """
    
    pass
    
    
    # set default download parent directory
    if download_root_dir==None:
        if version == 'v4r4':
            download_root_dir = join(expanduser('~'),'Downloads','ECCO_V4r4_PODAAC')
        elif version == 'v4r5':
            download_root_dir = join(expanduser('~'),'Downloads','ECCO_V4r5_PODAAC')

    # define the directory where the downloaded files will be saved
    download_dir = Path(download_root_dir) / ShortName
    
    # create the download directory if it does not already exist
    if isdir(download_dir) == True:
        print(f'Download to directory {download_dir}')
    else:
        print(f'Creating download directory {download_dir}')
    download_dir.mkdir(exist_ok = True, parents=True)
    
    # query CMR for granules matching the request
    urls,sizes = ecco_podaac_query(ShortName,StartDate,EndDate,version,snapshot_interval)
    
    # Download the granules
    
    downloaded_files = download_files_wrapper(urls, download_dir, n_workers, force_redownload, show_noredownload_msg)
    
    if return_downloaded_files:
        if len(downloaded_files) == 1:
            # if only 1 file is downloaded, return a string of filename instead of a list
            downloaded_files = downloaded_files[0]
        return downloaded_files
    


###================================================================================================================


def ecco_podaac_download_diskaware(ShortNames,StartDate,EndDate,version,snapshot_interval=None,\
                                   download_root_dir=None,max_avail_frac=0.5,n_workers=6,force_redownload=False,show_noredownload_msg=True):
    
    """
    
    This function estimates the storage footprint of ECCO datasets, given ShortName(s), a date range, and which 
    files (if any) are already present.
    If the footprint of the files to be downloaded (not including files already on the instance or re-downloads) 
    is <= the max_avail_frac specified of the environment's available storage, they are downloaded to the user's machine.
    Otherwise, an error is returned.
    
    Parameters
    ----------
    
    ShortNames: str or list, the ShortName(s) that identify the dataset on PO.DAAC.
    
    StartDate,EndDate: str, in 'YYYY', 'YYYY-MM', or 'YYYY-MM-DD' format, 
                       define date range [StartDate,EndDate] for download.
                       EndDate is included in the time range (unlike typical Python ranges).
                       ECCOv4r4 date range is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable closed budgets
                       within the specified date range.
    
    version: ('v4r4', 'v4r5'), the version of ECCO to download.
    
    snapshot_interval: ('monthly', 'daily', or None), if snapshot datasets are included in ShortNames, 
                       this determines whether snapshots are included for only the beginning/end of each month 
                       ('monthly'), or for every day ('daily').
                       If None or not specified, defaults to 'daily' if any daily mean ShortNames are included 
                       and 'monthly' otherwise.
    
    download_root_dir: str, defines parent directory to download files to.
                       Files will be downloaded to directory download_root_dir/ShortName/.
                       If not specified, parent directory defaults to '~/Downloads/ECCO_V4r4_PODAAC/',
                       or '~/Downloads/ECCO_V4r5_PODAAC/' if version == 'v4r5'.

    max_avail_frac: float, maximum fraction of remaining available disk space to use in storing current ECCO datasets.
                    If storing the datasets exceeds this fraction, an error is returned.
                    Valid range is [0,0.9]. If number provided is outside this range, it is replaced by the closer 
                    endpoint of the range.
    
    n_workers: int, number of workers to use in concurrent downloads. Benefits typically taper off above 5-6.
               Applies only if files are downloaded.
    
    force_redownload: bool, if True, existing files will be redownloaded and replaced;
                            if False, existing files will not be replaced.
                            Applies only if files are downloaded.
    
    show_noredownload_msg: bool, if True (default), and force_redownload=False, 
                               display message for each file that is already 
                               downloaded (and therefore not re-downloaded); 
                               if False, these messages are not shown.
    
    
    Returns
    -------
    downloaded_files: dict, with keys: ShortNames and values: downloaded or opened file(s) with path on local machine,
                     that can be passed directly to xarray (open_dataset or open_mfdataset).
    
    """

    pass

    import shutil
    
    # force max_avail_frac to be within limits [0,0.9]
    max_avail_frac = np.fmin(np.fmax(max_avail_frac,0),0.9)
    
    # determine value of snapshot_interval if None or not specified
    if snapshot_interval == None:
        snapshot_interval = 'monthly'
        for curr_shortname in ShortNames:
            if 'DAILY' in curr_shortname:
                snapshot_interval = 'daily'
                break

    # set default download parent directory
    if download_root_dir==None:
        if version == 'v4r4':
            download_root_dir = join(expanduser('~'),'Downloads','ECCO_V4r4_PODAAC')
        elif version == 'v4r5':
            download_root_dir = join(expanduser('~'),'Downloads','ECCO_V4r5_PODAAC')

    # add up total size of files that would be downloaded
    dataset_sizes = np.array([])
    urls_list_all = []
    for curr_shortname in ShortNames:
        
        # get list of files
        urls,sizes = ecco_podaac_query(curr_shortname,StartDate,EndDate,version,snapshot_interval)
        
        # define the directory where the downloaded files will be saved
        download_dir = Path(download_root_dir) / curr_shortname
        
        # create the download directory if it does not already exist
        if isdir(download_dir) == True:
            print(f'Download to directory {download_dir}')
        else:
            print(f'Creating download directory {download_dir}')
        download_dir.mkdir(exist_ok = True, parents=True)
        
        # compute size of current dataset
        curr_dataset_size = 0
        for url,size in zip(urls,sizes):
            if not isfile(join(download_dir,basename(url))):
                curr_dataset_size += size

        dataset_sizes = np.append(dataset_sizes,curr_dataset_size)
        urls_list_all.append(urls)
            

    # query available disk space at download location
    query_disk_completed = False
    query_dir = [download_root_dir][0]
    while query_disk_completed == False:
        try:
            avail_storage = shutil.disk_usage(query_dir).free
            query_disk_completed = True
        except:
            try:
                query_dir = join(*os.path.split(query_dir)[:-1])
            except:                
                print('Error: can not detect available disk space for download_root_dir: '+download_root_dir)
                return -1

    # fraction of available storage that would be occupied by downloads
    sizes_sum = np.sum(dataset_sizes)
    avail_frac = sizes_sum/avail_storage

    print(f'Size of files to be downloaded to instance is {(1.e-3)*np.round((1.e3)*sizes_sum/(2**30))} GB,\n'\
                +f'which is {.01*np.round((1.e4)*avail_frac)}% of the {(1.e-3)*np.round((1.e3)*avail_storage/(2**30))} GB available storage.')

    downloaded_files = {}
    if avail_frac <= max_avail_frac:
        # proceed with file downloads
        print('Proceeding with file downloads via NASA Earthdata URLs')
        for curr_shortname,urls_list in zip(ShortNames,urls_list_all):
        
            # define the directory where the downloaded files will be saved
            download_dir = Path(download_root_dir) / curr_shortname
            
            # download files
            curr_downloaded_files = download_files_wrapper(urls_list, download_dir, version, n_workers, force_redownload, show_noredownload_msg)

            if len(curr_downloaded_files) == 1:
                # if only 1 file is downloaded, return a string of filename instead of a list
                curr_downloaded_files = curr_downloaded_files[0]

            downloaded_files[curr_shortname] = curr_downloaded_files

    else:
        # raise exception and prevent download
        raise Exception('Download size is larger than specified fraction of available storage.\n'\
              +'Please modify your request to reduce its storage footprint.')

    return downloaded_files


    
    
###================================================================================================================


def ecco_podaac_download_subset(ShortName,StartDate=None,EndDate=None,snapshot_interval=None,\
                                n_workers=4,force_redownload=False,show_noredownload_msg=True,\
                                vars_to_include='all',vars_to_omit=None,\
                                times_to_include='all',\
                                k_isel=[0,50,1],tile_isel=[0,13,1],j_isel=[0,90,1],i_isel=[0,90,1],\
                                Z_isel=[0,50,1],latitude_isel=[0,360,1],longitude_isel=[0,720,1],\
                                netcdf4=True,include_latlon_coords=True,\
                                download_or_list='download',\
                                list_filename='files_to_download.txt',\
                                download_root_dir=None,subset_file_id='',\
                                return_downloaded_files=False):
  
    """
    
    Downloads subsets of ECCOv4r4 datasets from PO.DAAC using OPeNDAP.
    This routine downloads ECCO datasets from PO.DAAC. It is adapted by Andrew Delman from the 
    ecco_podaac_download routine derived from the Jupyter notebooks created by Jack McNelis and Ian Fenty,
    with some code from the OPeNDAP subsetting download script by Toshio Mike Chin and Y. Jiang 
    (https://github.com/nasa/podaac_tools_and_services/blob/master/subset_opendap/subset_dataset.py).
    
    Parameters
    ----------
    
    ShortName: str, the ShortName that identifies the dataset on PO.DAAC.
    
    StartDate,EndDate: str, in 'YYYY', 'YYYY-MM', or 'YYYY-MM-DD' format, 
                       define date range [StartDate,EndDate] for download.
                       EndDate is included in the time range (unlike typical Python ranges).
                       ECCOv4r4 date range is '1992-01-01' to '2017-12-31'.
                       For 'SNAPSHOT' datasets, an additional day is added to EndDate to enable closed budgets
                       within the specified date range.
                       If StartDate or EndDate are not specified, they are inferred from times_to_include;
                       if times_to_include is also not specified an error is returned.
    
    snapshot_interval: ('monthly', 'daily', or None), if ShortName refers to a snapshot dataset, 
                       this determines whether snapshots are included for only the beginning/end of each month 
                       ('monthly'), or for every day ('daily').
                       If None or not specified, defaults to the type of the date(s) in times_to_include,
                       or StartDate if times_to_include not specified.
    
    n_workers: int, number of workers to use in concurrent downloads.
    
    force_redownload: bool, if True, existing files will be redownloaded and replaced;
                            if False, existing files will not be replaced.
    
    show_noredownload_msg: bool, if True (default), and force_redownload=False, 
                               display message for each file that is already 
                               downloaded (and therefore not re-downloaded); 
                               if False, these messages are not shown.
    
    vars_to_include: list or tuple, names of data variables to include in the downloaded files.
                                    Dimension and coordinate variables are automatically included,
                                    except for the lat/lon coordinate variables when include_latlon_coords=False.
                                    Default is 'all', i.e., to include all data variables in the dataset.
    
    vars_to_omit: list or tuple, names of variables to exclude from the downloaded files.
                                 Default is None, i.e., to include all variables in the dataset.
                                 If both vars_to_include and vars_to_omit are specified,
                                 vars_to_include takes precedence, unless 
                                 vars_to_include='all' in which case vars_to_omit takes precedence.
    
    times_to_include: 'all' or list, tuple, or NumPy array.
                      Indicates the specific days or months to be downloaded, within the StartDate,EndDate 
                      time range specified previously.
                      If a list/tuple/NumPy array is given, it must consist either of strings of the format 
                      'YYYY', 'YYYY-MM', or 'YYYY-MM-DD', or of NumPy datetime64 objects, 
                      e.g., np.datetime64('YYYY-MM-DD').
                      This may be useful for downloading specific years, 
                      specific times of the year from multiple years, or specific days of the month.
                      If a 'YYYY' string or np.datetime64[Y] object is given, all months or days in the given year
                      will be included.
                      If a 'YYYY-MM' string or np.datetime64[M] object is given but the ShortName indicates 
                      daily temporal resolution, all of the days in that month will be included.
                      If a 'YYYY-MM-DD' string or np.datetime64[D] object is given but the ShortName indicates 
                      monthly temporal resolution, the given string/object will be truncated to 'YYYY-MM'.
                      For 'SNAPSHOT' datasets where a year/month string or np.datetime64 object type is included, 
                      the first of the following month will also be included 
                      (to enable budget closure for the last month).
                      Default is 'all', which downloads all files within the StartDate,EndDate time range.
    
    k_isel,tile_isel,j_isel,i_isel,
    Z_isel,latitude_isel,longitude_isel: 3-element lists, tuples, or NumPy arrays.
                                         Enables spatial subsetting, either in the native grid or lat/lon domain, 
                                         by defining the indices to download for each dimension
                                         in the format [start,end,stride] (using Python indexing conventions
                                         where 0 is the first index and end is not included).
                                         Note: only index ranges with consistent spacing can be downloaded 
                                         (e.g., downloading tiles 0,1,3,4 would need to be done either with
                                         tile_isel=[0,5,1] or as two separate downloads [0,2,1] and [3,5,1]).
                                         Defaults to the full range of each dimension.
                                         If indices are specified but the dimension does not exist in the files 
                                         (e.g., tile_isel is specified but the ShortName is for a lat/lon regridded
                                         dataset), the index specification is ignored.
    
    netcdf4: bool, indicates whether to download files as NetCDF4 or (classic) NetCDF3 files.
    
    include_latlon_coords: bool, indicates whether to include lat/lon coordinate variables in the 
                           native grid downloaded files.
                           Default is True. For the download of a large number of files (especially daily files),
                           False is recommended to reduce the size of the download.
                           Use the grid file instead to obtain the lat/lon coordinate variables.
                           If downloading the grid file, or if downloading a lat/lon re-mapped data file, 
                           this option is ignored and the coordinates are included regardless.
    
    download_or_list: ('download' or 'list'), indicates whether to download the files,
                      or output download URLs to a text file to be downloaded later (e.g., using wget or curl).
                      Default is 'download'.
                      The options after this apply to either 'list' or 'download',
                      if not relevant they are ignored.
    
    if download_or_list == 'list':
        
        list_filename: str, path and filename of text file to write download URLs to.
                       Default is 'urls_to_download.txt' in the current working directory.
                       If list_filename already exists, output will be appended to existing file.
    
    if download_or_list == 'download':
        
        download_root_dir: str, defines parent directory to download files to.
                           Files will be downloaded to directory download_root_dir/ShortName/.
                           If not specified, parent directory defaults to '~/Downloads/ECCO_V4r4_PODAAC/',
                           or '~/Downloads/ECCO_V4r5_PODAAC/' if version == 'v4r5'.
        subset_file_id: str, identifier appended to each downloaded file to identify it as a subset.
                        Default is to not append an identifier.

        return_downloaded_files: bool, if True, string or list of downloaded file(s) (including files that were already on disk
                                 and not replaced) is returned.
                                 if False (default), the function returns nothing.

    Returns
    -------
    downloaded_files: str or list, downloaded file(s) with local path that can be passed 
                      directly to xarray (open_dataset or open_mfdataset).
                      Only returned if download_or_list='download' and return_downloaded_files=True.
    
    """
    
    pass
    
    import sys,os
    from datetime import date, timedelta
    from math import floor,ceil
    import requests
    import re

    import ftplib
    import numpy as np
    import pandas as pd

    if sys.version_info >= (3,0):
        import subprocess
        from urllib import request
    else:
        import commands
        import urllib
    from concurrent.futures import ThreadPoolExecutor
    from itertools import repeat
    from getpass import getpass
    from http.cookiejar import CookieJar
    from io import StringIO

    from pathlib import Path
    from platform import system
    from netrc import netrc
    from os.path import basename, isfile, isdir, join, expanduser
    from tqdm import tqdm
    import time    
    
    #=====================================================

    ### Subroutines

    ### Helper subroutine to log into NASA EarthData

    # not pretty but it works
    def setup_earthdata_login_auth(url: str='urs.earthdata.nasa.gov'):
    # look for the netrc file and use the login/password
        try:
            username, _, password = netrc(file=_netrc).authenticators(url)
            
        # if the file is not found, prompt the user for the login/password
        except (FileNotFoundError, TypeError):
            print('Please provide Earthdata Login credentials for access.')
            username, password = input('Username: '), getpass('Password: ')
            
        manager = request.HTTPPasswordMgrWithDefaultRealm()
        manager.add_password(None, url, username, password)
        auth = request.HTTPBasicAuthHandler(manager)
        jar = CookieJar()
        processor = request.HTTPCookieProcessor(jar)
        opener = request.build_opener(auth, processor)
        request.install_opener(opener)
        
    
    ### Helper subroutines to make the API calls to search CMR and parse response
    def set_params(params: dict):
        params.update({'scroll': "true", 'page_size': 2000})
        return {par: val for par, val in params.items() if val is not None}
    
    
    def get_results(params: dict, headers: dict=None):
        response = requests.get(url="https://cmr.earthdata.nasa.gov/search/granules.json", 
                                  params=set_params(params),
                                  headers=headers)
        return response, response.headers
    
    
    def get_granules(params: dict):
        import json

        response, headers = get_results(params=params)
        json_dict = json.loads(response.text)
        gran_listings = json_dict['feed']['entry']
        opendap_urls = list()
        for curr_gran in gran_listings:
            if "POCLOUD" in curr_gran['data_center']:
                for curr_links in curr_gran['links']:
                    if ('title' in curr_links.keys()) and ("OPeNDAP" in curr_links['title']):
                        opendap_urls.append(curr_links['href'])
        
        if len(opendap_urls)==0:
            raise Exception("No granules matched your input parameters.")
        
        return opendap_urls
    
    
    ### Create datetime arrays from times_to_include and granule URLs to parse which granules should be included in download
    
    def datetimes_to_include(times_to_include,ShortName,StartDate,EndDate,snapshot_interval):
        # create array of datetimes indicated by times_to_include
        if 'MONTHLY' in ShortName:
            include_datetimes = np.array([]).astype('datetime64[M]')
        else:
            include_datetimes = np.array([]).astype('datetime64[D]')
        for curr_time in times_to_include:
            if isinstance(curr_time,np.datetime64) == False:
                if len(curr_time) == 4:
                    curr_time = np.datetime64(curr_time,'Y')
                elif len(curr_time) == 7:
                    curr_time = np.datetime64(curr_time,'M')
                elif len(curr_time) == 10:
                    curr_time = np.datetime64(curr_time,'D')
                else:
                    curr_time = np.datetime64('NaT')
            if 'MONTHLY' in ShortName:
                if curr_time.dtype == 'datetime64[Y]':
                    include_datetimes = np.append(include_datetimes,\
                                                  np.arange(curr_time.astype('datetime64[M]'),\
                                                  (curr_time+np.timedelta64(1,'Y')).astype('datetime64[M]'),\
                                                  np.timedelta64(1,'M')))
                else:
                    include_datetimes = np.append(include_datetimes,\
                                                  curr_time.astype('datetime64[M]'))
            elif 'DAILY' in ShortName:
                if curr_time.dtype == 'datetime64[Y]':
                    include_datetimes = np.append(include_datetimes,\
                                                  np.arange(curr_time.astype('datetime64[D]'),\
                                                  (curr_time+np.timedelta64(1,'Y')).astype('datetime64[D]'),\
                                                  np.timedelta64(1,'D')))
                elif curr_time.dtype == 'datetime64[M]':
                    curr_month = int(str(curr_time)[5:7])
                    include_datetimes = np.append(include_datetimes,\
                                                  np.arange(curr_time.astype('datetime64[D]'),\
                                                  (curr_time+np.timedelta64(1,'M')).astype('datetime64[D]'),\
                                                  np.timedelta64(1,'D')))
                else:
                    include_datetimes = np.append(include_datetimes,\
                                                  curr_time.astype('datetime64[D]'))
            elif 'SNAPSHOT' in ShortName:
                # include first day of succeeding month for snapshot datasets
                if curr_time.dtype == 'datetime64[Y]':
                    include_datetimes = np.append(include_datetimes,\
                                                  np.arange(curr_time.astype('datetime64[M]'),\
                                                  (curr_time+np.timedelta64(1,'Y')).astype('datetime64[M]') \
                                                  + np.timedelta64(1,'M'),\
                                                  np.timedelta64(1,'M')))
                elif ((snapshot_interval == 'monthly') \
                  or (curr_time.dtype == 'datetime64[M]') or (len(StartDate) <= 7)):

                    curr_month = int(str(curr_time)[5:7])
                    include_datetimes = np.append(include_datetimes,\
                                                  np.arange(curr_time.astype('datetime64[M]'),\
                                                  curr_time.astype('datetime64[M]') \
                                                  + np.timedelta64(2,'M'),\
                                                  np.timedelta64(1,'M')).astype('datetime64[D]'))
                else:
                    include_datetimes = np.append(include_datetimes,\
                                                  curr_time.astype('datetime64[D]'))
        
        return include_datetimes
    
    def datetimes_grans(grans_urls):
        # create array of datetimes (if present) associated with granules
        for url_count,url in enumerate(grans_urls):
            if 'MONTHLY' in ShortName:
                if url_count == 0:
                    grans_datetimes = np.array([]).astype('datetime64[M]')
                    include_datetimes = np.array([]).astype('datetime64[M]')
                curr_datetime_month = re.findall("_[0-9]{4}-[0-9]{2}",url)[0][1:]
                grans_datetimes = np.append(grans_datetimes,\
                                            np.datetime64(curr_datetime_month,'M'))
            elif (('DAILY' in ShortName) or ('SNAPSHOT' in ShortName)):
                if url_count == 0:
                    grans_datetimes = np.array([]).astype('datetime64[D]')
                    include_datetimes = np.array([]).astype('datetime64[D]')
                curr_datetime_day = re.findall("_[0-9]{4}-[0-9]{2}-[0-9]{2}",url)[0][1:]
                grans_datetimes = np.append(grans_datetimes,\
                                            np.datetime64(curr_datetime_day,'D'))
            else:
                if url_count == 0:
                    grans_datetimes = np.array([])
                    include_datetimes = np.array([])
                grans_datetimes = np.append(grans_datetimes,np.nan)
        
        return grans_datetimes
    
    
    def get_variable_info(url: str):
        url_appended=url+'.info'
        response = requests.get(url=url_appended)
        response_linesplit = response.text.splitlines()
        varinfo_dict = {}
        varinfo_section = False
        for line in response_linesplit:
            if 'Variables in this Dataset' in line:
                varinfo_section = True
            if varinfo_section == True:
                if "<b>" in line:
                    varname_start = line.find("<b>") + 3
                    varname_end = line.find("</b>")
                    varname = line[varname_start:varname_end]
                if ("[" in line) and ("]" in line) and (".." in line):
                    if "<br>" in line:
                        line = line.split("<br>")[0]
                        dims = line.split("[")[1:]
                        dims_dict = {}
                        for dim in dims:
                            dim_split = dim.split("=")
                            dimname = dim_split[0][:-1]
                            dim_bounds = dim_split[1].split("]")[0]
                            dim_bounds = dim_bounds.split("..")
                            dim_bounds[0] = int(dim_bounds[0][1:])
                            dim_bounds[1] = int(dim_bounds[1])
                            dims_dict = {**dims_dict,**{dimname:dim_bounds}}
                        varinfo_dict = {**varinfo_dict,**{varname:dims_dict}}
        
        return varinfo_dict
    
    
    def indices_dim_append(indices_str,input_ind,bounds):
        index_first = np.fmax(input_ind[0],bounds[0])
        index_last = np.fmin(input_ind[1]-1,bounds[-1])
        index_stride = input_ind[2]
        index_range = index_last - index_first
        if index_range % index_stride > 0:
            index_last = index_first + (index_stride*(index_range//index_stride))
        indices_append = "["+":".join([str(index_first),str(index_stride),str(index_last)])+"]"
        indices_str += indices_append

        return indices_str
    
    
    # Encode special characters in URLs to be used for downloads
    def encode_url(url: str):
        url = url.replace("%20","%2520")
        url = url.replace("[","%5B")
        url = url.replace("]","%5D")

        return url
    
    
    ### Helper subroutine to gracefully download single files and avoids re-downloading if file already exists.
    # To force redownload of the file, pass **True** to the boolean argument *force* (default **False**)\n,
    def download_file(url: str, output_file: str, force: bool=False):
        """url (str): the HTTPS url from which the file will download
        output_file (str): the filename (with path) of the output file
        force (bool): download even if the file exists locally already
        """
        
        output_dir,output_filename = os.path.split(output_file)
        if not isdir(output_dir):
            # if download directory does not already exist, create it
            os.mkdir(output_dir)
        
        # if the file has already been downloaded, skip    
        if isfile(output_file) and force is False:
            print(output_filename + ' already exists, and force=False, not re-downloading')
            return output_file,0
        
        with requests.get(url) as r:
            if not r.status_code // 100 == 2: 
                raise Exception(r.text)
                return 0
            else:
                with open(output_file, 'wb') as f:
                    total_size_in_bytes = len(r.content)
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk:
                            f.write(chunk)
        
        return output_file,total_size_in_bytes
    
    
    def download_wrapper(url: str, url_append: str, download_dir: str, subset_file_id: str,\
                         force: bool=False, show_noredownload_msg: bool=True):
        import os.path
        
        head, tail = os.path.split(url)
        ncout = join(download_dir,tail)
        
        if ncout.endswith('.bz2') or ncout.endswith('.gz'):
            ncout = ncout.rsplit( ".", 1 )[ 0 ]
        
        if len(subset_file_id) == 0:
            ncout += '.nc'
        else:
            ncout += '_'+subset_file_id+'.nc'
        
        url=url+url_append

        downloaded_files = []
        total_download_size_in_bytes = 0
        try:
            result = download_file(url=url, output_file=ncout, force=force,\
                                   show_noredownload_msg=show_noredownload_msg)
            downloaded_files.append(result[0])
            total_download_size_in_bytes += result[1]
            status_code = 0
        except:
            # retry download with progressively longer pauses between download attempts
            max_retries = 3
            n_retry = 1
            while n_retry <= max_retries:
                time.sleep(5*(n_retry**2))
                try:
                    result = download_file(url=url, output_file=ncout, force=force,\
                                           show_noredownload_msg=show_noredownload_msg)
                    downloaded_files.append(result[0])
                    total_download_size_in_bytes += result[1]
                    status_code = 0
                    break
                except:
                    n_retry += 1
            if n_retry > max_retries:
                status_code = -1
        
        return downloaded_files,total_download_size_in_bytes,status_code
    
    
    #=====================================================
    
    
    # retrieve NumPy datetimes associated with times_to_include
    subset_datetimes = False
    if isinstance(times_to_include,str) == False:
        subset_datetimes = True
    else:
        if times_to_include != 'all':
            subset_datetimes = True
    if subset_datetimes == True:
        include_datetimes = datetimes_to_include(times_to_include,\
                                                 ShortName,StartDate,EndDate,\
                                                 snapshot_interval)
    
    
    # # If no StartDate or EndDate provided, obtain them from times_to_include
    try:
        if StartDate == None:
            StartDate = str(np.nanmin(include_datetimes))
        if EndDate == None:
            EndDate = str(np.nanmax(include_datetimes))
    except:
        sys.exit('Error: must specify either StartDate and EndDate, or times_to_include')
    
    
    # # Adjust StartDate and EndDate to CMR query values
    StartDate,EndDate,SingleDay_flag = date_adjustment(ShortName,\
                                         StartDate,EndDate,CMR_query=True)
    
    
    # set default download parent directory
    if download_root_dir==None:
        if version == 'v4r4':
            download_root_dir = join(expanduser('~'),'Downloads','ECCO_V4r4_PODAAC')
        elif version == 'v4r5':
            download_root_dir = join(expanduser('~'),'Downloads','ECCO_V4r5_PODAAC')
    
    # define the directory where the downloaded files will be saved
    download_dir = Path(download_root_dir) / ShortName
    
    # create the download directory if it does not already exist
    if isdir(download_dir) == True:
        print(f'Download to directory {download_dir}')
    else:
        print(f'Creating download directory {download_dir}')
    download_dir.mkdir(exist_ok = True, parents=True)
        
    
    ### Log into Earthdata using your username and password
    
    # Predict the path of the netrc file depending on os/platform type.
    _netrc = join(expanduser('~'), "_netrc" if system()=="Windows" else ".netrc")

    # actually log in with this command:
    setup_earthdata_login_auth()
    

    ### Query the NASA Common Metadata Repository (CMR) to find the URL of every granule
    ### associated with the desired ECCO Dataset and date range of interest.

    # create a Python dictionary with our search criteria:  `ShortName` and `temporal`
    input_search_params = {'ShortName': ShortName,
                           'temporal': ",".join([StartDate, EndDate])}
    
    print ('\nPlease wait while program searches for the granules ...\n')

    # grans means 'granules', PO.DAAC's term for individual files in a dataset
    grans_urls = get_granules(input_search_params)
    
    # reduce granule URL list to single day if only one day in requested range
    if (('MONTHLY' in ShortName) or ('DAILY' in ShortName)):
        if ((SingleDay_flag == True) and (len(grans['Granule UR']) > 1)):
            grans_urls = grans_urls.sort()[-1]
    
    if len(grans_urls) == 0:
        sys.exit('No granules with OPeNDAP access found for dataset: '+ShortName+'\nProgram will exit now !\n')
    
    
    ### Further time granule selection
    if subset_datetimes == True:
        grans_datetimes = datetimes_grans(grans_urls)
        
        # cycle through granule URLs and remove if not in include_datetimes array
        grans_datetimes_copy = list(tuple(grans_datetimes))
        grans_urls_copy = list(tuple(grans_urls))
        gran_count = 0
        for gran_datetime,gran_url in zip(grans_datetimes_copy,grans_urls_copy):
            if gran_datetime not in include_datetimes:
                grans_datetimes = np.delete(grans_datetimes,gran_count)
                grans_urls.remove(gran_url)
            else:
                gran_count += 1
  
    num_grans = len(grans_urls)
    print (f'\nTotal number of matching granules: {num_grans}')  
    
    
    # get variable names and dimension info
    varinfo_dict = get_variable_info(grans_urls[0])
    
    
    ### remove variables from varinfo_dict that will not be downloaded
    dim_names = ['i','i_g','j','j_g','tile','k','k_l','k_u','k_p1',\
                 'latitude','longitude','Z','time','nb','nv']
    latlon_coord_names_in_llc90 = ['XC','YC','XG','YG','XC_bnds','YC_bnds']  
    coord_names = latlon_coord_names_in_llc90 + ['time_bnds','Zl','Zu','Zp1','Z_bnds',\
                                                 'latitude_bnds','longitude_bnds']
    dim_coord_names = dim_names + coord_names
    if vars_to_include != 'all':
        vars_to_include_all = dim_coord_names + vars_to_include
    varinfo_dict_copy = {**varinfo_dict,**{'dummy_key':0}}
    del varinfo_dict_copy['dummy_key']
    for varname in varinfo_dict_copy.keys():
        if varname in ['nb','nv']:
            # 2025-03-01 OpenDap fix: remove these dimension variables that seem to cause errors
            del varinfo_dict[varname]
            continue
        if ((vars_to_include != 'all') and (varname not in vars_to_include_all)):
            del varinfo_dict[varname]
        elif vars_to_omit != None:
            if ((vars_to_include == 'all') and (varname in vars_to_omit)):
                del varinfo_dict[varname]
        if ('GEOMETRY' not in ShortName) and (include_latlon_coords == False):
            # remove lon/lat coordinate variables    
            if varname in latlon_coord_names_in_llc90:
                del varinfo_dict[varname]
    
    
    ### complete Opendap URLs
    
    if netcdf4 == True:  
        url_append = '.dap.nc4?dap4.ce='
    else:
        url_append = '.dap.nc?dap4.ce='
    # append each variable to Opendap URL with subsetting
    for varname in varinfo_dict.keys():
        vardims_dict = varinfo_dict[varname]
        indices_str = ''
        for dimname in vardims_dict.keys():
            if dimname == 'time':
                bounds = vardims_dict[dimname]
                indices_str = indices_dim_append(indices_str,[bounds[0],bounds[1]+1,1],vardims_dict[dimname])
            if dimname in ['k','k_l','k_u']:
                indices_str = indices_dim_append(indices_str,k_isel,vardims_dict[dimname])
            if dimname == 'k_p1':
                k_p1_isel = [k_isel[0],k_isel[1]+1,k_isel[2]]
                indices_str = indices_dim_append(indices_str,k_p1_isel,vardims_dict[dimname])
            if dimname == 'tile':
                indices_str = indices_dim_append(indices_str,tile_isel,vardims_dict[dimname])
            if dimname in ['j','j_g']:
                indices_str = indices_dim_append(indices_str,j_isel,vardims_dict[dimname])
            if dimname in ['i','i_g']:
                indices_str = indices_dim_append(indices_str,i_isel,vardims_dict[dimname])
            if dimname == 'Z':
                indices_str = indices_dim_append(indices_str,Z_isel,vardims_dict[dimname])
            if dimname == 'latitude':
                indices_str = indices_dim_append(indices_str,latitude_isel,vardims_dict[dimname])
            if dimname == 'longitude':
                indices_str = indices_dim_append(indices_str,longitude_isel,vardims_dict[dimname])
            if dimname in 'nb':
                bounds = vardims_dict[dimname]
                indices_str = indices_dim_append(indices_str,[bounds[0],bounds[1]+1,1],bounds)
            if dimname in 'nv':
                bounds = vardims_dict[dimname]
                indices_str = indices_dim_append(indices_str,[bounds[0],bounds[1]+1,1],bounds)
        url_append += ('/'+varname+indices_str+';')
    
    url_append = url_append[0:(len(url_append)-1)]  # remove the extra ";" at the end.
    
    # convert URLs to encoded data URLs for download
    for url_ind,url in enumerate(grans_urls):
        grans_urls[url_ind] = encode_url(url)
    url_append = encode_url(url_append)
    
    
    ### Either output list of download URLs, or download files
    if download_or_list == 'list':
        with open(list_filename,'a') as f:
            for url in grans_urls:
                f.write(url+url_append+"\n")
    
        print("URL list written/appended to "+list_filename+".\n"\
              +"To download these files with wget,\n"\
              +"the bash shell script wget_download_fromlist.sh may be invoked, e.g.:\n\n"\
              +"bash ./wget_download_fromlist.sh -i "+list_filename+" \n"\
              +"-P ~/Downloads/ECCO_V4r4_PODAAC/"+ShortName+"/ \n"\
              +"-n "+subset_file_id+" -u username -p password")
    elif download_or_list == 'download':
        start = time.time()
        try:
            # use thread pool to download in parallel, with tqdm progress bar
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                results = list(tqdm(executor.map(download_wrapper, grans_urls, repeat(url_append),\
                                                 repeat(download_dir), repeat(subset_file_id),\
                                                 repeat(force_redownload), repeat(show_noredownload_msg)),\
                                    total=len(grans_urls), desc='DL Progress',\
                                    ascii=True, ncols=75, file=sys.stdout))
                downloaded_files = []
                download_sizes = []
                status_codes = []
                for result in results:
                    downloaded_files.append(result[0][0])
                    download_sizes.append(result[1])
                    status_codes.append(result[-1])
                total_download_size_in_bytes = np.sum(np.asarray(download_sizes))
                status_codes = np.asarray(status_codes)
        except:
            downloaded_files = []
            total_download_size_in_bytes = 0
            status_codes = np.array([]).astype('int32')
            for url in grans_urls:
                downloaded_file,download_size,status_code \
                  = download_wrapper(url,url_append,download_dir,subset_file_id,force_redownload,show_noredownload_msg)
                downloaded_files.append(downloaded_file[0])
                total_download_size_in_bytes += download_size
                status_codes = np.append(status_codes,status_code)
        
        end = time.time()
        total_time_download = end - start
        
        print('\n=====================================')
        print(f'total downloaded: {np.round(total_download_size_in_bytes/1e6,2)} Mb')
        print(f'avg download speed: {np.round(total_download_size_in_bytes/1e6/total_time_download,2)} Mb/s')
        print('Time spent = ' + str(total_time_download) + ' seconds')
        print('\n')
        
        # Display dates of granules that were not downloaded successfully
        status_codes_bad = (status_codes < 0).nonzero()[0]
        if len(status_codes_bad) > 0:
            datetimes_not_downloaded = []
            for datetime_ind in status_codes_bad:
                datetimes_not_downloaded.append(str(grans_datetimes[datetime_ind]))
            print('Granules from the following dates not downloaded successfully:')
            print(datetimes_not_downloaded)
            print('Please try downloading these granules again using times_to_include option.')

        
        if return_downloaded_files == True:
            if len(downloaded_files) == 1:
                # if only 1 file is downloaded, return a string of filename instead of a list
                downloaded_files = downloaded_files[0]
            return downloaded_files
        
    else:
        print("Error: Incorrect value of download_or_list.\n"\
              +"Please specify download_or_list = ('download' or 'list').")
