import sys
sys.path.append('..')

import os
import pygrib
import numpy as np
from osgeo import gdal
from osgeo import osr
import glob

import metpy.calc as metpcalc
from metpy.units import units

import ntpath

from datetime import datetime, timedelta

from pathlib import Path

from gribprocess.utils_s3_grib_download import s3_gfs_download

def is_hour_between(start, end, now):
    is_between = False

    is_between |= start <= now <= end
    is_between |= end < start and (start <= now or now <= end)

    return is_between


def date_decide(date0):
    hr=int(date0.strftime('%H'))
    if is_hour_between(20, 2, hr):
        date=date0- timedelta(days=1)
        date_fmt=date.strftime('%Y%m%d')
        srn='18'
    elif is_hour_between(2, 8, hr):
        date=date0
        date_fmt=date.strftime('%Y%m%d')
        srn='00'
    elif is_hour_between(8, 14, hr):
        date=date0
        date_fmt=date.strftime('%Y%m%d')
        srn='06'
    else:
        date=date0
        date_fmt=date.strftime('%Y%m%d')
        srn='12'
    return date_fmt,srn



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



def gdal_tiff_creator(array,west,north,south,east,outputfile):
    #function from https://gis.stackexchange.com/a/37431
    xmin,ymin,xmax,ymax =west,south,east,north
    nrows,ncols = np.shape(array)
    xres = (xmax-xmin)/float(ncols)
    yres = (ymax-ymin)/float(nrows)
    geotransform=(xmin,xres,0,ymax,0, -yres)  
    output_raster = gdal.GetDriverByName('GTiff').Create(outputfile,ncols, nrows, 1 ,gdal.GDT_Float32)  # Open the file
    output_raster.SetGeoTransform(geotransform)  # Specify its coordinates
    srs = osr.SpatialReference()                 # Establish its coordinate encoding
    srs.ImportFromEPSG(4326)                     # This one specifies WGS84 lat long.
    output_raster.SetProjection( srs.ExportToWkt() )   # Exports the coordinate system 
    output_raster.GetRasterBand(1).WriteArray(array)   # Writes my array to the raster
    output_raster.FlushCache()
    
    
def getFeatures(gdf):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gdf.to_json())['features'][0]['geometry']]



def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def report_index_extent(rds, region_extent):
    array_y = np.asarray(rds['y'])
    array_x = np.asarray(rds['x'])
    idx = (np.abs(array_x - region_extent['ll_lon'])).argmin()
    west=array_x[idx]
    idx = (np.abs(array_y - region_extent['ur_lat'])).argmin()
    north=array_y[idx]
    idx = (np.abs(array_y - region_extent['ll_lat'])).argmin()
    south=array_y[idx]
    idx = (np.abs(array_x - region_extent['ur_lon'])).argmin()
    east=array_x[idx]
    output=(west,north,south,east)
    return output





def temp_grib_processor(grib_filename):
    myfile = pygrib.open(grib_filename)
    grb_obj = myfile.select(name='2 metre temperature')[0]
    forecast_time=grb_obj.validDate.strftime('%Y%m%d%H')
    gen_time=grb_obj.analDate
    fc_time=grb_obj.validDate
    lats, lons = grb_obj.latlons()
    array0=grb_obj.values
    array=array0-273.15
    south,north,west,east=lats.min(),lats.max(),lons.min(),lons.max()
    folderpath=os.path.dirname(grib_filename)
    outputfile=f'{folderpath}/gfs_{forecast_time}.tif'
    array1=np.round(array, 1)
    gdal_tiff_creator(array1,west,north,south,east,outputfile)
    #gfs_180=gfs_oneeightzero_conversion(outputfile)
    myfile.close()
    return outputfile,gen_time,fc_time



def prate_grib_processor(grib_filename):
    myfile = pygrib.open(grib_filename)
    grb_obj = myfile.select(name='Precipitation rate')[0]
    forecast_time=grb_obj.validDate.strftime('%Y%m%d%H')
    gen_time=grb_obj.analDate
    fc_time=grb_obj.validDate
    lats, lons = grb_obj.latlons()
    array0=grb_obj.values
    array=array0*3600
    south,north,west,east=lats.min(),lats.max(),lons.min(),lons.max()
    folderpath=os.path.dirname(grib_filename)
    outputfile=f'{folderpath}/gfs_{forecast_time}.tif'
    array1=np.round(array, 1)
    gdal_tiff_creator(array1,west,north,south,east,outputfile)
    #gfs_180=gfs_oneeightzero_conversion(outputfile)
    myfile.close()
    return outputfile,gen_time,fc_time


def wind_speed_grib_processor(grib_filename):
    myfile = pygrib.open(grib_filename)
    grb_obj = myfile.select(name='Wind speed (gust)')[0]
    forecast_time=grb_obj.validDate.strftime('%Y%m%d%H')
    gen_time=grb_obj.analDate
    fc_time=grb_obj.validDate
    lats, lons = grb_obj.latlons()
    array0=grb_obj.values
    array=array0
    south,north,west,east=lats.min(),lats.max(),lons.min(),lons.max()
    folderpath=os.path.dirname(grib_filename)
    outputfile=f'{folderpath}/gfs_{forecast_time}.tif'
    array1=np.round(array, 1)
    gdal_tiff_creator(array1,west,north,south,east,outputfile)
    #gfs_180=gfs_oneeightzero_conversion(outputfile)
    myfile.close()
    return outputfile,gen_time,fc_time



def calculate_winddirection(myfile):
    u_grb_obj = myfile.select(name='10 metre U wind component')[0]
    u_array0=u_grb_obj.values
    v_grb_obj = myfile.select(name='10 metre V wind component')[0]
    v_array0=v_grb_obj.values
    u_pdata3=u_array0* units.meter / units.second
    v_pdata3=v_array0* units.meter / units.second
    winddirection=metpcalc.wind_direction(u_pdata3,v_pdata3, convention='from')
    return winddirection.magnitude




def wind_direction_grib_processor(grib_filename):
    myfile = pygrib.open(grib_filename)
    grb_obj = myfile.select(name='10 metre U wind component')[0]
    forecast_time=grb_obj.validDate.strftime('%Y%m%d%H')
    gen_time=grb_obj.analDate
    fc_time=grb_obj.validDate
    lats, lons = grb_obj.latlons()
    array=calculate_winddirection(myfile)
    south,north,west,east=lats.min(),lats.max(),lons.min(),lons.max()
    folderpath=os.path.dirname(grib_filename)
    outputfile=f'{folderpath}/gfs_dircn_{forecast_time}.tif'
    array1=np.round(array, 0)
    gdal_tiff_creator(array1,west,north,south,east,outputfile)
    #gfs_180=gfs_oneeightzero_conversion(outputfile)
    myfile.close()
    return outputfile,gen_time,fc_time


def wind_vector_grib_processor(grib_filename,variable_name):
    myfile = pygrib.open(grib_filename)
    grb_obj = myfile.select(name=variable_name)[0]
    forecast_time=grb_obj.validDate.strftime('%Y%m%d%H')
    gen_time=grb_obj.analDate
    fc_time=grb_obj.validDate
    lats, lons = grb_obj.latlons()
    array0=grb_obj.values
    array=array0
    south,north,west,east=lats.min(),lats.max(),lons.min(),lons.max()
    folderpath=os.path.dirname(grib_filename)
    vec_var=variable_name.split(' ')[2]
    outputfile=f'{folderpath}/gfs_{vec_var}_{forecast_time}.tif'
    array1=np.round(array, 1)
    gdal_tiff_creator(array1,west,north,south,east,outputfile)
    #gfs_180=gfs_oneeightzero_conversion(outputfile)
    myfile.close()
    return outputfile,gen_time,fc_time


def cape_grib_processor(grib_filename):
    myfile = pygrib.open(grib_filename)
    grb_obj = myfile.select(name='Convective available potential energy')[0]
    forecast_time=grb_obj.validDate.strftime('%Y%m%d%H')
    gen_time=grb_obj.analDate
    fc_time=grb_obj.validDate
    lats, lons = grb_obj.latlons()
    array0=grb_obj.values
    array=array0
    south,north,west,east=lats.min(),lats.max(),lons.min(),lons.max()
    folderpath=os.path.dirname(grib_filename)
    outputfile=f'{folderpath}/gfs_{forecast_time}.tif'
    array1=np.round(array, 1)
    gdal_tiff_creator(array1,west,north,south,east,outputfile)
    #gfs_180=gfs_oneeightzero_conversion(outputfile)
    myfile.close()
    return outputfile,gen_time,fc_time



def variable_json_creator():
    temp={'varname':'temperature','s3_name':'temperature','label':'Temperature','shrt_name':'temp','unit':'Degree Celcius','gfs_level':'2_m','gfs_var':'tmp','gfs_fmt':'TMP','color':{'red':153,'green':0,'blue':0},'filter_string':':TMP:2 m above ground'}
    prate={'varname':'precipitation','s3_name':'precipitation','label':'Precipitation','shrt_name':'prate','unit':'mm.hr-1','gfs_level':'surface','gfs_var':'prate','gfs_fmt':'PRATE','color':{'red':0,'green':153,'blue':0},'filter_string':':PRATE:surface:'}
    wspd={'varname':'wind_speed','s3_name':'windspeed','label':'Wind speed','shrt_name':'wspd','unit':'m.s-1','gfs_level':'surface','gfs_var':'gust','gfs_fmt':'GUST','color':{'red':0,'green':153,'blue':0},'filter_string':':GUST:surface:'}
    wdir={'varname':'wind_direction','s3_name':'winddirection','label':'Wind direction','shrt_name':'wdir','unit':'angle in 360 degree','gfs_level':'10_m','gfs_var':'uvgrd','gfs_fmt':'UGRD:VGRD','color':{'red':0,'green':153,'blue':0},'filter_string':['UGRD:10 m above ground:','VGRD:10 m above ground:']}
    u_wcmpt={'varname':'wind_u_component','s3_name':'wind-vectors-uv','label':'wind_u_component','shrt_name':'u_cmpt','unit':'m.s-1','gfs_level':'10_m','gfs_var':'uvgrd','gfs_fmt':'UGRD:VGRD','color':{'red':0,'green':153,'blue':0},'filter_string':'UGRD:10 m above ground:'}
    v_wcmpt={'varname':'wind_v_component','s3_name':'wind-vectors-uv','label':'wind_v_component','shrt_name':'v_cmpt','unit':'m.s-1','gfs_level':'10_m','gfs_var':'uvgrd','gfs_fmt':'UGRD:VGRD','color':{'red':0,'green':153,'blue':0},'filter_string':'VGRD:10 m above ground:'}
    cape={'varname':'cape','s3_name':'lightning','label':'Lightning','shrt_name':'cape','unit':'cape_index_value','gfs_level':'surface','gfs_var':'cape','gfs_fmt':'CAPE','color':{'red':153,'green':0,'blue':0},'filter_string':':CAPE:surface:'}
    apcp={'varname':'precipitation','s3_name':'precipitation','label':'Precipitation','shrt_name':'apcp','unit':'mm.hr-1','gfs_level':'surface','gfs_var':'apcp','gfs_fmt':'APCP','color':{'red':0,'green':153,'blue':0},'filter_string':':APCP:surface:'}
    variable={}
    variable['temp']=[temp]
    variable['prate']=[prate]
    variable['wspd']=[wspd]
    variable['uvgrd']=[wdir,u_wcmpt,v_wcmpt]
    variable['cape']=[cape]
    variable['apcp']=[apcp]
    return variable

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


def gcs_uploader(params,gfs_180,bucket):
    filename=path_leaf(gfs_180)
    blob = bucket.blob(f'{params.startdate}{params.run}/{params.short_name}/{filename}')
    blob.upload_from_filename(gfs_180)


def export_temperature_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        #table_lastrow_id=return_table_last_id()
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/tmp/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        gfs_180,gen_time,fc_time=temp_grib_processor(grib_filename)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)        
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')


def export_precipitation_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/prate/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        gfs_180,gen_time,fc_time=prate_grib_processor(grib_filename)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')

            
#grib_array_extract_tiff(params,tiff_outputfile_3857,var_nctime,region,decrypted_credntials,databaseip,databasename)
def export_wind_speed_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        #table_lastrow_id=return_table_last_id()
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/gust/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        gfs_180,gen_time,fc_time=wind_speed_grib_processor(grib_filename)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)        
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')


def export_wind_direction_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        #table_lastrow_id=return_table_last_id()
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/uvgrd/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        gfs_180,gen_time,fc_time=wind_direction_grib_processor(grib_filename)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')


            

def export_wind_u_component_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        #table_lastrow_id=return_table_last_id()
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/uvgrd/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        var_name='10 metre U wind component'
        gfs_180,gen_time,fc_time=wind_vector_grib_processor(grib_filename,var_name)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')

def export_wind_v_component_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        #table_lastrow_id=return_table_last_id()
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/uvgrd/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        var_name='10 metre V wind component'
        gfs_180,gen_time,fc_time=wind_vector_grib_processor(grib_filename,var_name)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')


def export_cape_grib(params,bucket):
    fmt_fulltimestep=[str(i).zfill(3) for i in params.fulltimestep]
    for timestep in fmt_fulltimestep:
        #table_lastrow_id=return_table_last_id()
        grib_filepath=f'/home/gfs_weather_data/gfs_data/{params.startdate}{params.run}/cape/'
        grib_filename=f'{grib_filepath}gfs.t{params.run}z.pgrb2.0p25.f{timestep}'
        gfs_180,gen_time,fc_time=cape_grib_processor(grib_filename)
        if params.test_implmentation==1:
            pass
        else:
            gcs_uploader(params,gfs_180,bucket)
        print('completed run on '+params.varname)
        print(f'uploaded {params.varname} for time step {timestep}')
        
            



def pygrib_download_process(params,bucket):
    datefmt=params.startdate
    run=params.run
    var_name=params.short_name
    folderpath='/home/gfs_weather_data/gfs_data/'
    date_folder=folderpath+datefmt+run+'/'
    foldercreator(date_folder)
    variable_json=variable_json_creator()
    filter_json_list=variable_json[var_name]
    for filtr_json in filter_json_list:
        if filtr_json['shrt_name']=='temp':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_temperature_grib(params,bucket)
            #s3_temp_upload(params,mh_region_list)
        elif filtr_json['shrt_name']=='prate':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_precipitation_grib(params,bucket)
            #s3_prate_upload(params,mh_region_list)
        elif filtr_json['shrt_name']=='wspd':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_wind_speed_grib(params,bucket)
            #s3_wspd_upload(params,mh_region_list)
        elif filtr_json['shrt_name']=='wdir':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_wind_direction_grib(params,bucket)
            #s3_wdir_upload(params,mh_region_list)
        elif filtr_json['shrt_name']=='u_cmpt':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_wind_u_component_grib(params,bucket)
            #s3_uv_wind_upload(params,mh_region_list)
        elif filtr_json['shrt_name']=='v_cmpt':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_wind_v_component_grib(params,bucket)
        elif filtr_json['shrt_name']=='cape':
            var_folder=date_folder+filtr_json['gfs_var']+'/'
            foldercreator(var_folder)
            params.localfolder=var_folder
            params.searchString=filtr_json['filter_string']
            s3_gfs_download(params)
            export_cape_grib(params,bucket)
