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
import os
from tkinter.constants import SW
from matplotlib import pyplot as plt




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



if __name__ == '__main__':
    mainDir = '/scratch/erikj'
    mainGpmDir = '%s/Data/GPM' %mainDir
    mainClsDir = '%s/Data/CloudSat/2B-FLXHR-LIDAR.v05.02' %mainDir
    
    year = 2008
    month = 1
    day = 2
    hour_start = 1
    timediff = 10 * 60
    datum = datetime(year, month, day, hour_start)
    datum_m1 = datum - timedelta(1)
    datum_p1 = datum + timedelta(1)
    yday = datum.timetuple().tm_yday
    gpmDir = '%s/%d/%03d' %(mainGpmDir, year, yday)
    gpmname = '%s/merg_%d%02d%02d%02d_4km-pixel.nc4' %(gpmDir, year, month, day, hour_start)
    clsDir = '%s/%d/%d_%02d_%02d' %(mainClsDir, year, year, month, day)
    clsDir_m1 = '%s/%d/%d_%02d_%02d' %(mainClsDir, datum_m1.year, datum_m1.year, datum_m1.month, datum_m1.day)
    clsDir_p1 = '%s/%d/%d_%02d_%02d' %(mainClsDir, datum_p1.year, datum_p1.year, datum_p1.month, datum_p1.day)
    clsglob = glob.glob('%s/*.h5' %clsDir)
    clsglob = clsglob + glob.glob('%s/*.h5' %clsDir_m1)
    clsglob = clsglob + glob.glob('%s/*.h5' %clsDir_p1)
    clsglob.sort()

    gpmdata = readGPM(gpmname)
    tic = time.time()
    flxfiles = []
    for f in clsglob:
        fb = os.path.basename(f)
        y = fb[0:4]
        yd = fb[4:7]
        h = fb[7:9]
        m = fb[9:11]
        time_fn = datetime.strptime('%s %s %s %s' %(y, yd, h, m), '%Y %j %H %M')
        tsmin = datum - time_fn
        tsmax = time_fn - datum
        #: First control from file name
        if ((tsmin.days*24*3600 + tsmin.seconds) > (timediff + 2 * 3600)) or \
            ((tsmax.days*24*3600 + tsmax.seconds) > (2 * timediff + 0.5 * 3600)):
#             print('hmm')
            continue
        #: Read file
        liddata = readFLXHRLidar(f)
        #: Second control from time stamps
        if ((np.abs(liddata['sec1970'] - gpmdata['time_d'][0]*24*3600) <= timediff).sum() == 0) and \
            ((np.abs(liddata['sec1970'] - gpmdata['time_d'][1]*24*3600) <= timediff).sum() == 0):
#             print('ha')
#             pdb.set_trace()
            continue
        
        flxfiles.append(f)
        print('way')
    print(time.time() - tic)
    
    tree = spatial.KDTree(list(zip(gpmdata['lats'].ravel(), gpmdata['lons'].ravel())))
    sw = None
    tb = None
    for flxfile in flxfiles:
        liddata = readFLXHRLidar(flxfile)
#         lid_latlon = np.zeros([2, len(liddata['latitude'])])
#         lid_latlon[0, :] = liddata['latitude']
#         lid_latlon[1, :] = liddata['longitude']
        for ti in [0, 1]:
            
            timeInd = np.abs(liddata['sec1970'] - gpmdata['time_d'][ti]*24*3600) <= timediff
            if timeInd.sum() == 0:
                continue
            if sw is None:
                sw = liddata['QR'][0,timeInd,:]
                lw = liddata['QR'][1,timeInd,:]
                height = liddata['Height'][timeInd,:]
            else:
                sw = np.concatenate((sw, liddata['QR'][0,timeInd,:]))
                lw = np.concatenate((lw, liddata['QR'][1,timeInd,:]))
                height = np.concatenate((height, liddata['Height'][timeInd,:]))
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
        
        #     latsrav = lats.ravel()
        #     lonsrav = lons.ravel()
        #     
        #     test=np.zeros([2, len(latsrav)])
        #     test[0,:] = latsrav
        #     test[1,:] = lonsrav
            #: https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.spatial.KDTree.query.html#scipy.spatial.KDTree.query
            tic = time.time()
            minkded, minkdep = tree.query(lid_latlon)
            minkdea = np.unravel_index(minkdep, gpmdata['lats'].shape)
            if tb is None:
                tb = gpmdata['Tb_d'][ti, minkdea[0], minkdea[1]]
            else:
                tb = np.concatenate((tb, gpmdata['Tb_d'][ti, minkdea[0], minkdea[1]]))
            print(time.time() - tic)

    useInd = ~((tb == -9999) | (sw[:, 40] == -9.99))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    H, xedges, yedges, im = ax.hist2d(tb[useInd], sw[:,40][useInd])#, range=rang, bins=bins, norm=LogNorm(), vmin=1., vmax=1e4, cmap=cmap)
    fig.show()
    pdb.set_trace()
#             minkded, minkdep = tree.query(pt)
    
    
    
    
    