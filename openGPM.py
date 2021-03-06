#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Created on 2021-05-27

Copyright (c) 2021 Erik Johansson

Author(s):
Erik Johansson <erik.johansson@lmd.ipsl.fr>
'''

import numpy as np
import pdb
from datetime import datetime, timedelta
import time
import netCDF4 as nc  # @UnresolvedImport
import h5py  # @UnresolvedImport
from scipy import spatial  # @UnresolvedImport
import glob
import sys
import os


import matplotlib  # @UnresolvedImport
matplotlib.use('Agg')
from matplotlib import pyplot as plt  # @UnresolvedImport




def convertDType(val):
#        name = val.dtype.name
    if val.dtype.names is None:
        testType = val.dtype
    else:
        testType = val.dtype[0]
    dtype = np.dtype(testType.type)
    promoType = np.promote_types(testType, testType) #@UndefinedVariable
    if dtype == promoType:
        if dtype == np.dtype('float32'):
            dtype = np.dtype('float')
        re = val.astype(dtype)
    else:
        re = val
    
    return re


def makeDataUsefull(data, nv = False, attrs = [], ntb = True):
    val = data[()]
    val = convertDType(val)
    if attrs == []:
        try:
            factor = data.attrs['factor']
        except:
            factor = 1
        try:
            offset = data.attrs['offset']
        except:
            offset = 0
        try:
            lv = data.attrs['valid_range'][0]
            hv = data.attrs['valid_range'][1]
        except:
            lv = np.nan
            hv = np.nan
        try:
            missing = data.attrs['missing']
        except:
            missing = np.nan
#         data.attrs['_FillValue']
    else:
        try:
            missing = attrs[0][attrs[1] + '.missing'][()][0][0]
        except:
            missing = np.nan
        try:
            lv = attrs[0][attrs[1] + '.valid_range'][()][0][0][0]
            hv = attrs[0][attrs[1] + '.valid_range'][()][0][0][1]
        except:
            lv = np.nan
            hv = np.nan
        factor = attrs[0][attrs[1] + '.factor'][()][0][0]
        offset = attrs[0][attrs[1] + '.offset'][()][0][0]
    if nv:
        # TODO: I remove error message
        np.seterr(invalid='ignore')
        val = np.where(val == missing, np.nan, val)
        if ntb:
            val = np.where(val < lv, np.nan, val)
            val = np.where(val > hv, np.nan, val)
        np.seterr(invalid='warn')
    val = val * 1. / factor
    val = val + offset
#     if attrs[1] in ['QR', 'TOACRE', 'BOACRE']:
#         pdb.set_trace()
    return val


def convertTime(cloudsat):
    dsec = time.mktime((1993,1,1,0,0,0,0,0,0)) - time.timezone
    start_sec1970 = cloudsat['TAI_start'] + dsec
#        start_time = time.gmtime(start_sec1970[0])
    cloudsat['sec1970'] = cloudsat['Profile_time'] + start_sec1970
    return cloudsat



def readFLXHRLidar(filename):
    retv = {'longitude': None, 'latitude': None, 'Height': None, 'TAI_start': None, \
            'UTC_start': None, 'elevation': None, 'Profile_time': None, 'sec1970': None, \
            'Albedo': None, 'FD': None, 'FD_NA': None, 'FD_NC': None, 'FU': None, \
            'FU_NA': None, 'FU_NC': None, 'QR': None, 'RH': None, 'TOACRE': None, \
            'BOACRE': None, 'FlagCounts': None, 'Land_Char': None, 'Scene_status': None, \
            'MeanOLR': None, 'MeanOSR': None, 'MeanQLW': None, 'MeanQSW': None, \
            'MeanSFCE': None, 'MeanSFCR': None, 'MeanSLR': None, 'MeanSSR': None, \
            'Meansolar': None, 'SigmaOLR': None, 'SigmaOSR': None, 'SigmaQLW': None, \
            'SigmaQSW': None, 'SigmaSFCE': None, 'SigmaSFCR': None, 'SigmaSLR': None, \
            'SigmaSSR': None, 'Sigmasolar': None, 'Solar_zenith_angle': None, 'Status': None}
    newName = {'Longitude':'longitude', 'Latitude': 'latitude', 'DEM_elevation': 'elevation'}
    
    try:
        h5file = h5py.File(filename, 'r')
    except:
        print('Stupied file')
        print(filename)
        return -1
    root = "2B-FLXHR-LIDAR"
    groups = ['Geolocation Fields', 'Data Fields']
    for group in groups:
        groupFields = h5file["%s/%s" % (root, group)]
        for dataset in groupFields.keys():
            if (dataset in retv.keys()) or (dataset in newName.keys()):
                if "Swath Attributes" not in h5file[root]:
                    print('No Attributes')
                    print(filename)
                    h5file.close()
                    return -1
                else:
                    ntb = True
                    if dataset in ['TOACRE', 'BOACRE']:
                        ntb = False
                    val = makeDataUsefull(groupFields[dataset], nv = True, attrs = [h5file["%s/%s" % (root, "Swath Attributes")], dataset], ntb = ntb)
                
                if dataset in newName.keys():
                    retv[newName[dataset]] = val
                else:
                    retv[dataset] = val
            else:
                continue

    retv = convertTime(retv)
    h5file.close()
    return retv


def readCLDLidar(filename, latlon = False):
    retv = {'CL': None, 'CLType': None, 'CLTop': None, 'CLTopF': None, 'CLBase': None}#, 'CF': None}
    h5file = h5py.File(filename, 'r')
    root = "2B-CLDCLASS-LIDAR"
    group = 'Data Fields'
    retv['CLType'] = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'CloudLayerType')], nv = True)
    retv['CLTop'] = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'CloudLayerTop')], nv = True)
    retv['CLTopF'] = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'CloudLayerTop')])
    retv['CLBase'] = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'CloudLayerBase')], nv = True)
    retv['CL'] = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'Cloudlayer')], nv = True)
    if latlon:
        group = 'Geolocation Fields'
        lat = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'Latitude')], nv = True)
        lon = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'Longitude')], nv = True)
        ele = makeDataUsefull(h5file["%s/%s/%s" % (root, group, 'DEM_elevation')], nv = True)
        retv.update({'latitude': lat})
        retv.update({'longitude': lon})
        retv.update({'elevation': ele})
       
    h5file.close() 
    return retv

def readGPM(filename):
    ncf = nc.Dataset(filename)
    ncf.variables.keys()
#         ncf.variables['Tb'], 'lon', 'lat', 'time'])
    ncTb = ncf.variables['Tb']
    Tb_d = ncf.variables['Tb'][:].data
    Tb_m = ncf.variables['Tb'][:].mask
    Tb_f = ncf.variables['Tb'][:].fill_value
     
    nclat = ncf.variables['lat']
    lat_d = ncf.variables['lat'][:].data
    lat_m = ncf.variables['lat'][:].mask
    lat_f = ncf.variables['lat'][:].fill_value
    
    nclon = ncf.variables['lon']
    lon_d = ncf.variables['lon'][:].data
    lon_m = ncf.variables['lon'][:].mask
    lon_f = ncf.variables['lon'][:].fill_value
    
    nctime = ncf.variables['time']
    time_d = ncf.variables['time'][:].data
    time_m = ncf.variables['time'][:].mask
    time_f = ncf.variables['time'][:].fill_value
    
    
    lats, lons = np.meshgrid(lat_d,lon_d, indexing='ij')
    t1 = time.gmtime(time_d[0]*24*3600)
    t2 = time.gmtime(time_d[1]*24*3600)
    ncf.close()
    retv = {'ncTb': ncTb, 'Tb_d': Tb_d, 'Tb_m': Tb_m, 'Tb_f': Tb_f, \
            'nclat': nclat, 'lat_d': lat_d, 'lat_m': lat_m, 'lat_f': lat_f, \
            'nclon': nclon, 'lon_d': lat_d, 'lon_m': lon_m, 'lon_f': lon_f, \
            'nctime': nctime, 'time_d': time_d, 'time_m': time_m, 'time_f': time_f, \
            'lats': lats, 'lons': lons, 't1': t1, 't2': t2}
    return retv


def cutCloudSatLat(cd, latc=91):
    lat = cd['latitude']
#     lon = cd['longitude']
    latm = (lat <= latc) & (lat >= (-1. * latc))
    for arname, data in cd.items():
        if 'start' in arname:
            continue
        if 'Mean' in arname:
            continue
        if 'Sigma' in arname:
            continue
        if arname in ['FlagCounts']:
            continue
        if data.ndim == 3:
            if data.shape[2] > 126:
                print('shape2')
                pdb.set_trace()
            if data[0].shape[0] == len(latm):
                temp = np.zeros([data.shape[0], latm.sum(), data.shape[2]]) 
                for i in range(data.shape[0]):
                    temp[i, :, :] = data[i, latm, :]
                cd[arname] = temp
            else:
                print('3Dim')
                pdb.set_trace()
        elif data.ndim == 2:
            if data.shape[0] == len(latm):
                temp = np.zeros([latm.sum(), data.shape[1]]) 
                temp[:, :] = data[latm, :]
                cd[arname] = temp
            elif data.shape[1] == len(latm):
                temp = np.zeros([data.shape[0], latm.sum()]) 
                temp[:, :] = data[:, latm]
                cd[arname] = temp
            else:
                print('2Dim')
                pdb.set_trace()
        else:
            if data.shape[0] == len(latm):
                temp = np.zeros([latm.sum()])
                temp[:] = data[latm]
                cd[arname] = temp
            else:
                print('1Dim')
                pdb.set_trace()

    return cd


def cutGPMLat(gpmd, latc=91):
    lats = gpmd['lats']
    lat = gpmd['lat_d']
#             lat = lats[:,0]
    latm = (lat <= latc) & (lat >= (-1. * latc))
#             map = (lats <= latc) & (lats >= (-1. * latc))
    for arname, data in gpmd.items():
        if 'nc' in arname:
            continue
        if '_' in arname:
            continue
        if arname in ['t1', 't2']:
            continue
        if data.ndim == 3:
            if data[0].shape == lats.shape:
                temp = np.zeros([data.shape[0], latm.sum(), lats.shape[1]]) 
                for i in range(data.shape[0]):
                    temp[i, :, :] = data[i, latm, :]
                gpmd[arname] = temp
            else:
                print('3Dim')
                pdb.set_trace()
        elif data.ndim == 2:
            if data.shape == lats.shape:
                temp = np.zeros([latm.sum(), lats.shape[1]]) 
                temp[:, :] = data[latm, :]
                gpmd[arname] = temp
            else:
                print('2Dim')
                pdb.set_trace()
        else:
            if data.shape == lat.shape:
                temp = np.zeros([latm.sum()])
                temp[:] = data[latm]
                gpmd[arname] = temp
            else:
                print('1Dim')
                pdb.set_trace()
                
    return gpmd


def loadTemps(tempnames):
    print('load %d files' %len(tempnames))
    for i, fn in enumerate(tempnames):
        lo = np.load(fn + '.npy', allow_pickle=True)[0]
        if i == 0:
            retv = lo
        else:
            for namn, val in retv.items():
                if (val.ndim == 1) or (val.shape[0] > val.shape[1]):
                    retv[namn] = np.concatenate((val, lo[namn]), axis=0)
                else:
                    retv[namn] = np.concatenate((val, lo[namn]), axis=1)
                   
    return retv['SW'], retv['SW_cle'], retv['SW_clo'], \
           retv['LW'], retv['LW_cle'], retv['LW_clo'], \
           retv['TB'], retv['CLL'], retv['rhSW'], \
           retv['tb'], retv['sw'], retv['lw'], retv['clL'], retv['rh_sw'], 


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--lat-bound", type=int, default=91,  
                        help="Latitude Boundary. Default=90")
#     parser.add_argument("-c1","--clim-max", type=int, default=500, 
#                         help="Max clim. Default=500")
#     parser.add_argument("-fn","--method-name", type=str, default='All', 
#                         help="name of the method used. Default=All")
#     parser.add_argument("-fm","--method-flag", action="store_true", default=False, 
#                         help="Do not show figure. Default=True (i.e. show figure)")
    parser.add_argument("-s","--show", action="store_false", default=True, 
                        help="Do not show figure. Default=True (i.e. show figure)")
    
    
    
    args = parser.parse_args()
    show = args.show
    
    tic_tot = time.time()
    homedir = os.getenv("HOME") 
    mainDir = '%s/Scratch' %homedir    #'/scratch/erikj' #'/scratchu/ejohansson' #'/scratch/erikj'
    mainGpmDir = '%s/Data/GPM' %mainDir
    mainFlxDir = '%s/Data/CloudSat/2B-FLXHR-LIDAR' %mainDir#'%s/Data/CloudSat/2B-FLXHR-LIDAR.v05.02_03' %mainDir
    tempDir = '%s/OpenGPM/TempFiles' %mainDir
    plotDir = '%s/OpenGPM/Plots' %mainDir
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    
    #: Some variables
    sv = 'v04'
    mainFlxDir = mainFlxDir + '.' + sv
    year = 2008
    month = 1
    day = 2
    hour_start = 8
    timediff = 15 * 60
    nrfiles = 500
    tempnames = []
    heighH = 20000
    lowH = 10000
    latbound = args.lat_bound
    loadtf = True
    height_lev = int((heighH - lowH) / 240.) + 1
    heights = np.arange(lowH, heighH + 240, 240)
    #: Datetime object
    datumST = datetime(year, month, day, hour_start)
    for i in range(nrfiles):
        print('File %i of %i' %(i, nrfiles))
        datum = datumST + timedelta(hours=(i))
        #: GPM        
        yday = datum.timetuple().tm_yday
        #: Dir where to find GPM data
        gpmDir = '%s/%d/%03d' %(mainGpmDir, year, yday)
        #: Filename of GPM data
        gpmname = '%s/merg_%d%02d%02d%02d_4km-pixel.nc4' %(gpmDir, datum.year, datum.month, datum.day, datum.hour)
        if latbound > 90:
            latbound_tn = 90
        else:
            latbound_tn = latbound
        tempname = '%s/%s_lat-%i' %(tempDir, os.path.basename(gpmname).replace('.nc4', ''), latbound_tn)
        tempnames.append(tempname)
        if (os.path.isfile(tempname + '.npy')) and (loadtf == True):
            continue
        #: CloudSat
        #: Incase the pm imediff puts the cloudsat on the other side of midnight
        #: Get the cloudsat from day before and after.
        #: Date minus 1 day. 
        datum_m1 = datum - timedelta(1)
        #: Date plus 1 day
        datum_p1 = datum + timedelta(1)
        #: Dir with cloudsat
        flxDir = '%s/%d/%d_%02d_%02d' %(mainFlxDir, datum.year, datum.year, datum.month, datum.day)
        flxDir_m1 = '%s/%d/%d_%02d_%02d' %(mainFlxDir, datum_m1.year, datum_m1.year, datum_m1.month, datum_m1.day)
        flxDir_p1 = '%s/%d/%d_%02d_%02d' %(mainFlxDir, datum_p1.year, datum_p1.year, datum_p1.month, datum_p1.day)
        #: All cloudsat files
        flxglob = glob.glob('%s/*.h5' %flxDir)
        flxglob = flxglob + glob.glob('%s/*.h5' %flxDir_m1)
        flxglob = flxglob + glob.glob('%s/*.h5' %flxDir_p1)
        flxglob.sort()
        #: Read the GPM file
        gpmdata = readGPM(gpmname)
        gpmdata = cutGPMLat(gpmdata, latbound)
        tic = time.time()
        #: Lets find the corresponding cloudsat files
        flxfiles = []
        cldfiles = []
        for flxF in flxglob:
            #: Date from cloudsat filename
            fb = os.path.basename(flxF)
            fd = os.path.dirname(flxF)
            y = fb[0:4]
            yd = fb[4:7]
            h = fb[7:9]
            m = fb[9:11]
            time_fn = datetime.strptime('%s %s %s %s' %(y, yd, h, m), '%Y %j %H %M')
            #: Make sure there are a Cloud file
            cldDir = fd.replace('2B-FLXHR-LIDAR', '2B-CLDCLASS-LIDAR')
            cldName = '%s/%s*2B-CLDCLASS-LIDAR*.h5' %(cldDir, fb.split('2B-FLXHR-LIDAR')[0])
            cldglob = glob.glob(cldName)
            if len(cldglob) == 0:
                continue
            else:
                cldF = cldglob[0]
            #: Get the difference between date and filename
            tsmin = datum - time_fn
            tsmax = time_fn - datum
            #: First control from file name
            #: Make a ruf selection so we do not need to spend time to open to many files
            if ((tsmin.days*24*3600 + tsmin.seconds) > (timediff + 2 * 3600)) or \
                ((tsmax.days*24*3600 + tsmax.seconds) > (2 * timediff + 0.5 * 3600)):
                continue
            #: Read cloudsat file
            liddata = readFLXHRLidar(flxF)
            liddata = cutCloudSatLat(liddata, latbound)
            #: Second control from time stamps
            #: Lets lock inside the file and remove (continu) if there is no data to be used.
            if ((np.abs(liddata['sec1970'] - gpmdata['time_d'][0]*24*3600) <= timediff).sum() == 0) and \
                ((np.abs(liddata['sec1970'] - gpmdata['time_d'][1]*24*3600) <= timediff).sum() == 0):
                continue
            #: Append files to be used
            flxfiles.append(flxF)
            cldfiles.append(cldF)
            print('way')
        print(time.time() - tic)
        if len(flxfiles) == 0:
            print('No matching CloudSat')
            #: Remove the tempname from tempnames
            if (tempname == np.asarray(tempnames)).any():
                tempnames.pop(np.where(tempname == np.asarray(tempnames))[0][0])
            continue
        tic = time.time()
        #: Create KDTree form gpm lats/lons
        tree = spatial.KDTree(list(zip(gpmdata['lats'].ravel(), gpmdata['lons'].ravel())))
        print(time.time() - tic)
        sw = None
        tb = None
        for i in [0]:#range(len(flxfiles)):
            #: Read the cloudsat data
            flxfile = flxfiles[i]
            cldfile = cldfiles[i]
            liddata = readFLXHRLidar(flxfile)
            liddata = cutCloudSatLat(liddata, latbound)
            clddata = readCLDLidar(cldfile, True)
            clddata = cutCloudSatLat(clddata, latbound)
    #         lid_latlon = np.zeros([2, len(liddata['latitude'])])
    #         lid_latlon[0, :] = liddata['latitude']
    #         lid_latlon[1, :] = liddata['longitude']
            #: There are 2 time stamps in each gpm file
            for ti in [0, 1]:
                #: Bothe cloudsat and gpm has time defined as time from 1970. Cloudsat in sec and gpm in days
                timeInd = np.abs(liddata['sec1970'] - gpmdata['time_d'][ti]*24*3600) <= timediff
                #: Different cloudsat files to different gpm time stamps
                if timeInd.sum() == 0:
                    continue
                if sw is None:
                    #: First file
                    sw = liddata['QR'][0,timeInd,:]
                    lw = liddata['QR'][1,timeInd,:]
                    height = liddata['Height'][timeInd,:]
                    rh_sw = liddata['RH'][0, timeInd]
                    clL = clddata['CL'][timeInd]
                    clT = clddata['CLTop'][timeInd,:]
                    clB = clddata['CLBase'][timeInd,:]
                else:
                    sw = np.concatenate((sw, liddata['QR'][0,timeInd,:]))
                    lw = np.concatenate((lw, liddata['QR'][1,timeInd,:]))
                    height = np.concatenate((height, liddata['Height'][timeInd,:]))
                    rh_sw = np.concatenate((rh_sw, liddata['RH'][0, timeInd]))
                    
                    clL = np.concatenate((clL, clddata['CL'][timeInd]))
                    clT = np.concatenate((clT, clddata['CLTop'][timeInd,:]))
                    clB = np.concatenate((clB, clddata['CLBase'][timeInd,:]))
                
                    
                #: cloudsat lat/lon needs to be in the format lists in list, [ [lat, lon], [lat, lon], ...]
                lid_latlon = []
                for i in range(len(liddata['latitude'][timeInd])):
                    lid_latlon.append([liddata['latitude'][timeInd][i],liddata['longitude'][timeInd][i]])
    
    #             lid_latlon = [ [liddata['latitude'][timeInd][i],liddata['longitude'][timeInd][i]] for i in range(len(liddata['latitude'][timeInd]))]
    #             pt = [[6, 30]]
    #             minp = []
    #             mina = []
    #             for i in range(len(lid_latlon)):
    #                 minp.append(np.argmin(np.sqrt(((gpmdata['lats'] - lid_latlon[i][0])**2) + ((gpmdata['lons'] - lid_latlon[i][1])**2))))
    #                 mina.append(np.unravel_index(minp[i], gpmdata['lats'].shape))
    #             latsrav = lats.ravel()
    #             lonsrav = lons.ravel()
    #              
    #             test=np.zeros([2, len(latsrav)])
    #             test[0,:] = latsrav
    #             test[1,:] = lonsrav
                
                #: Find the nearest neighburg
                #: https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
                minkded, minkdep = tree.query(lid_latlon)
                minkdea = np.unravel_index(minkdep, gpmdata['lats'].shape)
                if tb is None:
                    tb = gpmdata['Tb_d'][ti, minkdea[0], minkdea[1]]
                    gp = minkdep
                else:
                    tb = np.concatenate((tb, gpmdata['Tb_d'][ti, minkdea[0], minkdea[1]]))
                    gp = np.concatenate((gp, minkdep))
        
        #: Averaging
        #: heighest h in m - lowest h in m / 240 m. Add one to be sure to cover the whole wanted height
        if tb is None:
            tempnames.pop(-1)
            continue
        tic = time.time()
        SW = np.zeros([heights.shape[0], tb.shape[0]]) * np.nan
        LW = np.zeros([heights.shape[0], tb.shape[0]]) * np.nan
        SW_cle = np.zeros([heights.shape[0], tb.shape[0]]) * np.nan
        SW_clo = np.zeros([heights.shape[0], tb.shape[0]]) * np.nan
        LW_cle = np.zeros([heights.shape[0], tb.shape[0]]) * np.nan
        LW_clo = np.zeros([heights.shape[0], tb.shape[0]]) * np.nan
#         SW_cloudy
#         SW_cloudy_clear
        TB = np.zeros([tb.shape[0]])
        rhSW = np.zeros([tb.shape[0]])
        CLL = np.zeros([tb.shape[0]])
        t = 0
        for i in range(len(tb)):
            if i == 0:
    #             TB = [tb[i]]
    #             TB[0] = tb[0]
                sumSW = np.zeros((heights.shape))
                numSW = np.zeros((heights.shape))
                sumLW = np.zeros((heights.shape))
                numLW = np.zeros((heights.shape))
                
                sumSW_clo = np.zeros((heights.shape))
                numSW_clo = np.zeros((heights.shape))
                sumLW_clo = np.zeros((heights.shape))
                numLW_clo = np.zeros((heights.shape))
                
                sumSW_cle = np.zeros((heights.shape))
                numSW_cle = np.zeros((heights.shape))
                sumLW_cle = np.zeros((heights.shape))
                numLW_cle = np.zeros((heights.shape))
            else:
                if gp[i] != gp[i-1]:
                    TB[t] = tb[i-1]
                    rhSW[t] = rh_sw[i-1]
                    CLL[t] = clL[i-1]
#                     TB.append(tb[i])
                    SW[:, t] = sumSW / numSW
                    LW[:, t] = sumLW / numLW
                    
                    SW_clo[:, t] = sumSW_clo / numSW_clo
                    LW_clo[:, t] = sumLW_clo / numLW_clo
                    
                    SW_clo[:, t] = sumSW_clo / numSW_clo
                    LW_clo[:, t] = sumLW_clo / numLW_clo
                    
                    
                    t = t + 1
                    
                    sumSW = np.zeros((heights.shape))
                    numSW = np.zeros((heights.shape))
                    sumLW = np.zeros((heights.shape))
                    numLW = np.zeros((heights.shape))
                    
                    sumSW_clo = np.zeros((heights.shape))
                    numSW_clo = np.zeros((heights.shape))
                    sumLW_clo = np.zeros((heights.shape))
                    numLW_clo = np.zeros((heights.shape))
                    
                    sumSW_cle = np.zeros((heights.shape))
                    numSW_cle = np.zeros((heights.shape))
                    sumLW_cle = np.zeros((heights.shape))
                    numLW_cle = np.zeros((heights.shape))
                    
            if (tb[i] == -9999):
                continue
            
            for j in range(heights.shape[0]):
                hi = (height[i, :] >= heights[j]) & (height[i, :] < (heights[j] + 240))
                
                s = sw[i, hi]
                l = lw[i, hi]
    #             try:
                swn999ind = ~(s==-9.99)
                lwn999ind = ~(l==-9.99)
                s = s[swn999ind]
                l = l[swn999ind]
    #             if s != -9.99:
                sumSW[j] = sumSW[j] + s.sum()
                numSW[j] = numSW[j] + len(s)
                
    #             if s != -9.99:
                sumLW[j] = sumLW[j] + l.sum()
                numLW[j] = numLW[j] + len(s)
                
                clear = True
                if clL[i] != 0:
                    for k in range(int(clL[i])-1):
                        if not (~np.isnan(clB[i,:])).sum() == clL[i]:
                            print('check clL')
                            pdb.set_trace()
                        
                        if (height[i, hi][0] >= (clB[i, k] * 1000.)) and (height[i, hi][0] <= (clT[i, k] * 1000.)):
                            clear = False
                
                if clear:
                    sumSW_cle[j] = sumSW_cle[j] + s.sum()
                    numSW_cle[j] = numSW_cle[j] + len(s)
                    sumLW_cle[j] = sumLW_cle[j] + l.sum()
                    numLW_cle[j] = numLW_cle[j] + len(s)
                else:
                    sumSW_clo[j] = sumSW_clo[j] + s.sum()
                    numSW_clo[j] = numSW_clo[j] + len(s)
                    sumLW_clo[j] = sumLW_clo[j] + l.sum()
                    numLW_clo[j] = numLW_clo[j] + len(s)
                    
                        
                #: Find if pixel is cloudy
    #             except:
    #                 pdb.set_trace()
        TB[t] = tb[i]
        rhSW[t] = rh_sw[i]
        CLL[t] = clL[i]
    #     TB.append(tb[i])
        SW[:, t] = sumSW / numSW
        LW[:, t] = sumLW / numLW
        
        SW_cle[:, t] = sumSW_cle / numSW_cle
        LW_cle[:, t] = sumLW_cle / numLW_cle
        
        SW_clo[:, t] = sumSW_clo / numSW_clo
        LW_clo[:, t] = sumLW_clo / numLW_clo
        
        SW = SW[:, 0:t+1]
        LW = LW[:, 0:t+1]
        SW_cle = SW_cle[:, 0:t+1]
        LW_cle = LW_cle[:, 0:t+1]
        SW_clo = SW_clo[:, 0:t+1]
        LW_clo = LW_clo[:, 0:t+1]
        TB = TB[0:t+1]
        rhSW = rhSW[0:t+1]
        CLL = CLL[0:t+1]
        np.save(tempname, [{'SW': SW, 'SW_cle': SW_cle, 'SW_clo': SW_clo, \
                            'LW': LW, 'LW_cle': LW_cle, 'LW_clo': LW_clo, \
                            'TB': TB, 'CLL': CLL, 'rhSW': rhSW, \
                            'tb': tb, 'sw': sw, 'lw': lw, 'clL': clL, 'rh_sw': rh_sw }])
        print(time.time() - tic)
        print('save')
#         pdb.set_trace()
    
    SW, SW_cle, SW_clo, \
    LW, LW_cle, LW_clo, \
    TB, CLL, rhSW, \
    tb, sw, lw, clL, rh_sw = loadTemps(tempnames)

    avH = [] * height_lev
#     lev = 40
    
    tbMF = ~(tb == -9999)
    tbM = ~(TB == -9999)
    
    for cl in ['clo', 'cle', 'all']:
        if cl == 'clo':
            clMF = (clL != 0)
            clM = (CLL != 0)
        elif cl == 'cle':
            clMF = (clL == 0)
            clM = (CLL == 0)
        elif cl == 'all':
            clMF = np.ones(clL.shape).astype(bool)
            clM = np.ones(CLL.shape).astype(bool)
        
        for dn in ['d']:
            if dn == 'd':
                dnMF = rh_sw > 0 # day night mask full
                dnM = rhSW > 0
            for lev in [10, 20, 30, 40]:
                swMLev = ~(sw[:, lev] == -9.99)
                lwMLev = ~(lw[:, lev] == -9.99)
                useIndSw = tbMF & clMF & swMLev & dnMF
                useIndLw = tbMF & clMF & lwMLev & dnMF
                
                
                #: Just need edges
                fig = plt.figure()
                ax = fig.add_subplot(111)
                H, xedgesS, yedgesS, im = ax.hist2d(tb[useIndSw], sw[:, lev][useIndSw], bins=[50, 50])
                
                H, xedgesL, yedgesL, im = ax.hist2d(tb[useIndLw], lw[:, lev][useIndLw], bins=[50, 50], range= [[220, 300], [-1., 1.]])
    

                print('Nr of SW cases = %d' %useIndSw.sum())
    
#     xedgesS = np.asarray([209.  , 211.22, 213.44, 215.66, 217.88, 220.1 , 222.32, 224.54, 226.76, 228.98, 231.2 , 233.42, 235.64, 237.86, 240.08, 242.3 , 244.52, 246.74, 248.96, 251.18, 253.4 , 255.62, 257.84, 260.06, 262.28, 264.5 , 266.72, 268.94, 271.16, 273.38, 275.6 , 277.82, 280.04, 282.26, 284.48, 286.7 , 288.92, 291.14, 293.36, 295.58, 297.8 , 300.02, 302.24, 304.46, 306.68, 308.9 , 311.12, 313.34, 315.56, 317.78, 320.  ])
#     xedgesL = np.asarray([220. , 221.6, 223.2, 224.8, 226.4, 228. , 229.6, 231.2, 232.8, 234.4, 236. , 237.6, 239.2, 240.8, 242.4, 244. , 245.6, 247.2, 248.8, 250.4, 252. , 253.6, 255.2, 256.8, 258.4, 260. , 261.6, 263.2, 264.8, 266.4, 268. , 269.6, 271.2, 272.8, 274.4, 276. , 277.6, 279.2, 280.8, 282.4, 284. , 285.6, 287.2, 288.8, 290.4, 292. , 293.6, 295.2, 296.8, 298.4, 300. ])
    
                SW_mean = np.zeros((heights.shape[0], xedgesS.shape[0]-1))
                LW_mean = np.zeros((heights.shape[0], xedgesL.shape[0]-1))
                SW_medi = np.zeros((heights.shape[0], xedgesS.shape[0]-1))
                LW_medi = np.zeros((heights.shape[0], xedgesL.shape[0]-1))
                SW_std = np.zeros((heights.shape[0], xedgesS.shape[0]-1))
                LW_std = np.zeros((heights.shape[0], xedgesL.shape[0]-1))
                tbS_mean = np.zeros((xedgesL.shape[0]-1))
                tbL_mean = np.zeros((xedgesL.shape[0]-1))
#                 for i in range(len(xedgesS)-1):
#                     use = (TB >= xedgesS[i]) & (TB < xedgesS[i+1])
#                     SW_mean[:,i] = np.nanmean(SW[:,use], axis=1)
#                     SW_medi[:,i] = np.nanmedian(SW[:,use], axis=1)
#                     SW_std[:,i] = np.nanstd(SW[:,use], axis=1)
#                     tbS_mean[i] = (xedgesS[i] + xedgesS[i + 1]) /2.
#                 for i in range(len(xedgesL)-1):
#                     use = (TB >= xedgesL[i]) & (TB < xedgesL[i+1])
#                     LW_mean[:,i] = np.nanmean(LW[:,use], axis=1)
#                     LW_medi[:,i] = np.nanmedian(LW[:,use], axis=1)
#                     LW_std[:,i] = np.nanstd(LW[:,use], axis=1)
#                     tbL_mean[i] = (xedgesL[i] + xedgesL[i + 1]) /2.
                
                for i in range(len(xedgesS)-1):
                    use = (tb[useIndSw] >= xedgesS[i]) & (tb[useIndSw] < xedgesS[i+1])
                    SW_mean[:,i] = np.nanmean(sw[useIndSw, lev][use])#, axis=1)
                    SW_medi[:,i] = np.nanmedian(sw[useIndSw, lev][use])#, axis=1)
                    SW_std[:,i] = np.nanstd(sw[useIndSw, lev][use])#, axis=1)
                    tbS_mean[i] = (xedgesS[i] + xedgesS[i + 1]) /2.
                for i in range(len(xedgesL)-1):
                    use = (tb[useIndSw] >= xedgesL[i]) & (tb[useIndSw] < xedgesL[i+1])
                    LW_mean[:,i] = np.nanmean(lw[useIndLw, lev][use])#, axis=1)
                    LW_medi[:,i] = np.nanmedian(lw[useIndLw, lev][use])#, axis=1)
                    LW_std[:,i] = np.nanstd(lw[useIndLw, lev][use])#, axis=1)
                    tbL_mean[i] = (xedgesL[i] + xedgesL[i + 1]) /2.
                
                fig = plt.figure()
                fig.suptitle('Height = %d m, lev = %d, clear/cloudy = %s' %(heights[lev], lev, cl))
                ax = fig.add_subplot(1,2,1)
                H, xedgesS, yedgesS, im = ax.hist2d(tb[useIndSw], sw[:,lev][useIndSw], bins=[50, 50])#, range= [[220, 300], [0., 0.05]])#, range=rang, bins=bins, norm=LogNorm(), vmin=1., vmax=1e4, cmap=cmap)
                ax.plot(tbS_mean, SW_mean[lev, :], '-r')
                ax.axhline(np.nanmean(SW_mean[lev, :]), color='g', ls='-')
                ax.axhline(np.nanmean(SW_mean[lev, :]) + np.nanstd(SW_mean[lev, :]), color='g', ls='--')
                ax.axhline(np.nanmean(SW_mean[lev, :]) - np.nanstd(SW_mean[lev, :]), color='g', ls='--')
#                 ax.plot(tbS_mean, (SW_mean[lev, :] + SW_std[lev, :]), '--r')
#                 ax.plot(tbS_mean, (SW_mean[lev, :] - SW_std[lev, :]), '--r')
                ax.set_xlabel('GPM - Tb [K]')
                ax.set_ylabel('Cloudsat - SW [K/day]')
            #     ax.set_title('Height = %d m, lev = %d' %(height[useInd, lev][0], lev))
                
                
                print('Nr of LW cases = %d' %useIndLw.sum())
                ax = fig.add_subplot(1,2,2)
                H, xedgesL, yedgesL, im = ax.hist2d(tb[useIndLw], lw[:, lev][useIndLw], bins=[50, 50], range= [[220, 300], [-1., 1.]])#, range= [[220, 300], [0., 0.05]])#, range=rang, bins=bins, norm=LogNorm(), vmin=1., vmax=1e4, cmap=cmap)
                ax.plot(tbL_mean, LW_mean[lev, :], '-r')
                ax.plot(tbL_mean, (LW_mean[lev, :] + LW_std[lev, :]), '--r')
                ax.plot(tbL_mean, (LW_mean[lev, :] - LW_std[lev, :]), '--r')
                ax.axhline(np.nanmean(LW_mean[lev, :]), color='g', ls='-')
                ax.axhline(np.nanmean(LW_mean[lev, :]) + np.nanstd(LW_mean[lev, :]), color='g', ls='--')
                ax.axhline(np.nanmean(LW_mean[lev, :]) - np.nanstd(LW_mean[lev, :]), color='g', ls='--')
                ax.set_xlabel('GPM - Tb [K]')
                ax.set_ylabel('Cloudsat - LW [K/day]')
#                 ax.set_title('Height = %d m, lev = %d' %(heights[lev], lev))
    
                figname = '%s/histo_%s_%s_%d_l-%d' %(plotDir, cl, dn, lev, latbound_tn)
                fig.savefig(figname + '.png')
                if show:
                    fig.show()
    if show:
        pdb.set_trace()
    sys.exit()
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    im = ax.imshow(SW_mean[17:, :], origin='lower', aspect=(SW_mean[17:, :].shape[1]/SW_mean[17:, :].shape[0]))#, aspect=0.75)
    ax.set_title('Cloudsat - SW [K/day]')
    ax.set_ylabel('Height [m]')
    ax.set_xlabel('GPM - Tb [K]')
    yticks = [0, 10, 20]#, 30, 40]
    xticks = [0, 10, 20, 30, 40, 49]
    ax.set_yticks(yticks)
    ax.set_yticklabels(heights[17:][yticks].astype(int))#['10', '12', '15', '17', '20'])
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_xticklabels(tbS_mean[xticks].astype(int))
    
    
    cbar = fig.colorbar(im, orientation='horizontal')#, cax=cbar_ax, ticks=barticks)
    
    ax = fig.add_subplot(1,2,2)
    im = ax.imshow(LW_mean[17:, :], origin='lower', aspect=(LW_mean[17:, :].shape[1]/LW_mean[17:, :].shape[0]))#, aspect=0.75)
    ax.set_title('Cloudsat - LW [K/day]')
    ax.set_ylabel('Height [m]')
    ax.set_xlabel('GPM - Tb [K]')
    ax.set_yticks(yticks)
    ax.set_yticklabels(['', '', ''])
#     ax.set_yticklabels(['', '', '', '', ''])#['10', '12', '15', '17', '20'])
    ax.set_xticks([0, 10, 20, 30, 40, 50])
    ax.set_xticklabels(tbL_mean[xticks].astype(int))
    cbar = fig.colorbar(im, orientation='horizontal')#, cax=cbar_ax, ticks=barticks)
    figname = '%s/swlw' %plotDir
    fig.savefig(figname + '.png')
    fig.show()
    
    
    print(time.time() - tic_tot)
    pdb.set_trace()
#             minkded, minkdep = tree.query(pt)
    
    
    
    
    