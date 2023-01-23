def ecco_podaac_download(ShortName,StartDate,EndDate,download_root_dir=None,n_workers=6,force_redownload=False):
    """
    This routine downloads ECCO datasets from PO.DAAC. It is adapted from the Jupyter notebooks created by Jack McNelis and Ian Fenty (https://github.com/ECCO-GROUP/ECCO-ACCESS/blob/master/PODAAC/Downloading_ECCO_datasets_from_PODAAC/README.md) and modified by Andrew Delman (https://ecco-v4-python-tutorial.readthedocs.io).
    
    Parameters
    ----------
    ShortName: the ShortName of the dataset (can be identified from https://search.earthdata.nasa.gov/search?fpj=ECCO, selecting the "i" information button and the ShortName will appear in a gray box in the upper-left corner)
    
    StartDate: the start of the time range to be downloaded, expressed in the format "YYYY-MM-DD"
    
    EndDate: the end of the time range to be downloaded, expressed in the format "YYYY-MM-DD"
    
    download_root_dir: path of the parent directory to download ECCO files
    
    n_workers: number of workers to use in concurrent downloads
    
    force_redownload: if True, existing files will be redownloaded and replaced; if False, existing files will not be replaced
    """
    
    
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
    
    # if no download directory specified, set directory under user's home directory
    if download_root_dir==None:
        user_home_dir = expanduser('~')
        download_root_dir = join(user_home_dir,'Downloads','ECCO_V4r4_PODAAC')
    
    # Predict the path of the netrc file depending on os/platform type.
    _netrc = join(expanduser('~'), "_netrc" if system()=="Windows" else ".netrc")
    
    ## Define Helper Subroutines
    
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
        response = requests.get(url="https://cmr.earthdata.nasa.gov/search/granules.csv", 
                                params=set_params(params),
                                headers=headers)
        return response, response.headers
    
    
    def get_granules(params: dict):
        response, headers = get_results(params=params)
        scroll = headers['CMR-Scroll-Id']
        hits = int(headers['CMR-Hits'])
        if hits==0:
            raise Exception("No granules matched your input parameters.")
        df = pd.read_csv(StringIO(response.text)) 
        while hits > df.index.size:
            response, _ = get_results(params=params, headers={'CMR-Scroll-Id': scroll})
            data = pd.read_csv(StringIO(response.text))
            df = pd.concat([df, data])
        return df
    
    ### Helper subroutine to gracefully download single files and avoids re-downloading if file already exists.
    # To force redownload of the file, pass **True** to the boolean argument *force* (default **False**)\n,
    def download_file(url: str, output_dir: str, force: bool=False):
        """url (str): the HTTPS url from which the file will download
        output_dir (str): the local path into which the file will download
        force (bool): download even if the file exists locally already
        """
        if not isdir(output_dir):
            raise Exception(f"Output directory doesnt exist! ({output_dir})")
        
        target_file = join(output_dir, basename(url))
        
        # if the file has already been downloaded, skip    
        if isfile(target_file) and force is False:
            print(f'\n{basename(url)} already exists, and force=False, not re-downloading')
            return 0
        
        with requests.get(url) as r:
            if not r.status_code // 100 == 2: 
                raise Exception(r.text)
                return 0
            else:
                with open(target_file, 'wb') as f:
                    total_size_in_bytes= int(r.headers.get('content-length', 0))
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk:
                            f.write(chunk)
    
                    return total_size_in_bytes
    
    ### Helper subroutine to download all urls in the list `dls`
    def download_files_concurrently(dls, download_dir, force=False):
        start_time = time.time()
    
        # use 3 threads for concurrent downloads
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
    
            # tqdm makes a cool progress bar
            results = list(tqdm(executor.map(download_file, dls, repeat(download_dir), repeat(force)), total=len(dls)))
        
            # add up the total downloaded file sizes
            total_download_size_in_bytes = np.sum(np.array(results))
            # calculate total time spent in the download
            total_time = time.time() - start_time
    
            print('\n=====================================')
            print(f'total downloaded: {np.round(total_download_size_in_bytes/1e6,2)} Mb')
            print(f'avg download speed: {np.round(total_download_size_in_bytes/1e6/total_time,2)} Mb/s')
    

    # define the directory where the downloaded files will be saved
    download_dir = Path(download_root_dir) / ShortName
    
    # create the download directory
    download_dir.mkdir(exist_ok = True, parents=True)
    
    print(f'created download directory {download_dir}')
    
    ## Log into Earthdata using your username and password
    
    # actually log in with this command:
    setup_earthdata_login_auth()
    
    # Query the NASA Common Metadata Repository to find the URL of every granule associated with the desired ECCO Dataset and date range of interest.
    
    # create a Python dictionary with our search criteria:  `ShortName` and `temporal`
    input_search_params = {'ShortName': ShortName,
                           'temporal': ",".join([StartDate, EndDate])}
    
    print(input_search_params)
    
    ### Query CMR for the desired ECCO Dataset
    
    # grans means 'granules', PO.DAAC's term for individual files in a dataset
    grans = get_granules(input_search_params)
    
    # grans.info()
    
    num_grans = len( grans['Granule UR'] )
    print (f'\nTotal number of matching granules: {num_grans}')
    
    
    ## Download the granules
    
    # convert the rows of the 'Online Access URLS' column to a Python list
    dls = grans['Online Access URLs'].tolist()
    
    try:
        # Attempt concurrent downloads, but if error arises switch to sequential downloads
        ### Method 1: Concurrent downloads
        
        # Define the maximum number of concurrent downloads (benefits typically taper off above 5-6)
        max_workers = 6
        
        # Force redownload (or not) depending on value of force_redownload
        download_files_concurrently(dls, download_dir, force_redownload)
        
    except:
        ### Method 2: Sequential Downloads
        
        # Download each URL sequentially in a for loop.
        total_download_size_in_bytes = 0
        start_time = time.time()
        
        # loop through all urls in dls
        for u in dls:
            u_name = u.split('/')[-1]
            print(f'downloading {u_name}')
            total_download_size_in_bytes += download_file(url=u, output_dir=download_dir, force=force_redownload)
        
        # calculate total time spent in the download
        total_time = time.time() - start_time
        
        print('\n=====================================')
        print(f'total downloaded: {np.round(total_download_size_in_bytes/1e6,2)} Mb')
        print(f'avg download speed: {np.round(total_download_size_in_bytes/1e6/total_time,2)} Mb/s')
