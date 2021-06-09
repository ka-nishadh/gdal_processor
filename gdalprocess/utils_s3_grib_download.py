import os
import re


import ntpath
import boto3


def path_leaf(path):
    """
    Get the name of a file without any extension from given path

    Parameters
    ----------
    path : file full path with extension
    
    Returns
    -------
    str
       filename in the path without extension

    """
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
    


s3_client = boto3.client('s3')


def foldercreator(path):
   """
    creates a folder

    Parameters
    ----------
    path : folder path
            
    Returns
    -------
    creates a folder
    """
   if not os.path.exists(path):
        os.makedirs(path)


def s3botodownloadfile(s3keylocation, bucket,localfolder):
    """
    Download file in S3

    Parameters
    ----------
    s3keylocation : key location in s3 bucket
    bucket : s3 bucket with target contents
    localfolder : folder path in local to download the file
    
    Returns
    -------
    None

    """
    filename=path_leaf(s3keylocation)
    #filename=configz.date+'_pfz_gen_logs.txt' 
    s3_client.download_file(bucket,s3keylocation, '%s/%s' % (localfolder,filename))
    download_filepath=localfolder+filename
    return download_filepath, filename


def get_byte_range(lines,searchstring_list):
    byte_range_list=[]
    for searchstring in searchstring_list:
        expr = re.compile(searchstring)
        byte_ranges = {}
        for n, line in enumerate(lines, start=1):
            # n is the line number (starting from 1) so that when we call for 
            # `lines[n]` it will give us the next line. (Clear as mud??)
    
            # Use the compiled regular expression to search the line
            if expr.search(line):   
                # aka, if the line contains the string we are looking for...
    
                # Get the beginning byte in the line we found
                parts = line.split(':')
                rangestart = int(parts[1])
    
                # Get the beginning byte in the next line...
                if n+1 < len(lines):
                    # ...if there is a next line
                    parts = lines[n].split(':')
                    rangeend = int(parts[1])
                else:
                    # ...if there isn't a next line, then go to the end of the file.
                    rangeend = ''
    
                # Store the byte-range string in our dictionary, 
                # and keep the line information too so we can refer back to it.
                byte_ranges['start_byte']=rangestart
                byte_ranges['end_byte']=rangeend
                byte_ranges['parse_line'] = line
                byte_range_list.append(byte_ranges)
    if len(byte_range_list) ==2:
        pbyte_range={}
        pbyte_range['start_byte']=byte_range_list[0]['start_byte']
        pbyte_range['end_byte']=byte_range_list[1]['end_byte']
        pbyte_range['parse_line'] = byte_range_list[1]['parse_line']
    else:
        pbyte_range=byte_range_list[0]
    return pbyte_range
            



def s3_gfs_download(params):
    s3keylocation=f'gfs.{params.startdate}/{params.run}/atmos/gfs.t{params.run}z.pgrb2.0p25.f{params.timestep}'
    s3keylocation_idx=f'{s3keylocation}.idx'
    download_filepath, filename=s3botodownloadfile(s3keylocation_idx, params.awsbucketname,params.localfolder)
    with open(download_filepath) as f:
        lines = [line.rstrip() for line in f]
        
    # Search 
    if isinstance(params.searchString, list):
        searchstring_list=params.searchString
        try:
            byte_ranges=get_byte_range(lines,searchstring_list)
        except Exception as e:
            pass
    else:
        searchstring_list=[params.searchString]
        try:
            byte_ranges=get_byte_range(lines,searchstring_list)
        except Exception as e:
            pass

    gribfilename=f'gfs.t{params.run}z.pgrb2.0p25.f{params.timestep}'
    try:
        obj = s3_client.get_object(
            Bucket=params.awsbucketname,
            Key=s3keylocation,
            Range='bytes={}-{}'.format(byte_ranges['start_byte'], byte_ranges['end_byte']))
        with open(params.localfolder+gribfilename, 'wb') as f:
            for chunk in obj['Body'].iter_chunks(chunk_size=4096):
                f.write(chunk)
        return gribfilename
    except Exception as e:
        pass


