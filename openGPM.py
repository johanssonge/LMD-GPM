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
    tic_tot = time.time()
    mainDir = '/scratchu/ejohansson' #'/scratch/erikj'
    mainGpmDir = '%s/Data/GPM' %mainDir
    mainClsDir = '%s/Data/CloudSat/2B-FLXHR-LIDAR.v05.02' %mainDir
    tempDir = '%s/OpenGPM/TempFiles' %mainDir
    if not os.path.isdir(tempDir):
        os.makedirs(tempDir)
    
    #: Some variables
    year = 2008
    month = 1
    day = 2
    hour_start = 8
    timediff = 15 * 60
    nrfiles = 2
    tempnames = []
    heighH = 20000
    lowH = 10000
    height_lev = int((heighH - lowH) / 240.) + 1
    heights = np.arange(lowH, heighH + 240, 240)
    #: Datetime object
    datumST = datetime(year, month, day, hour_start)
    for i in range(nrfiles):
        datum = datumST + timedelta(hours=(i))
        #: GPM        
        yday = datum.timetuple().tm_yday
        #: Dir where to find GPM data
        gpmDir = '%s/%d/%03d' %(mainGpmDir, year, yday)
        #: Filename of GPM data
        gpmname = '%s/merg_%d%02d%02d%02d_4km-pixel.nc4' %(gpmDir, datum.year, datum.month, datum.day, datum.hour)
        
        tempname = '%s/%s' %(tempDir, os.path.basename(gpmname).replace('.nc4', ''))
        tempnames.append(tempname)
        if os.path.isfile(tempname + '.npy'):
            continue
        #: CloudSat
        #: Incase the pm imediff puts the cloudsat on the other side of midnight
        #: Get the cloudsat from day before and after.
        #: Date minus 1 day. 
        datum_m1 = datum - timedelta(1)
        #: Date plus 1 day
        datum_p1 = datum + timedelta(1)
        #: Dir with cloudsat
        clsDir = '%s/%d/%d_%02d_%02d' %(mainClsDir, year, year, month, day)
        clsDir_m1 = '%s/%d/%d_%02d_%02d' %(mainClsDir, datum_m1.year, datum_m1.year, datum_m1.month, datum_m1.day)
        clsDir_p1 = '%s/%d/%d_%02d_%02d' %(mainClsDir, datum_p1.year, datum_p1.year, datum_p1.month, datum_p1.day)
        #: All cloudsat files
        clsglob = glob.glob('%s/*.h5' %clsDir)
        clsglob = clsglob + glob.glob('%s/*.h5' %clsDir_m1)
        clsglob = clsglob + glob.glob('%s/*.h5' %clsDir_p1)
        clsglob.sort()
        
        #: Read the GPM file
        gpmdata = readGPM(gpmname)
        tic = time.time()
        #: Lets find the corresponding cloudsat files
        flxfiles = []
        for f in clsglob:
            #: Date from cloudsat filename
            fb = os.path.basename(f)
            y = fb[0:4]
            yd = fb[4:7]
            h = fb[7:9]
            m = fb[9:11]
            time_fn = datetime.strptime('%s %s %s %s' %(y, yd, h, m), '%Y %j %H %M')
            #: Get the difference between date and filename
            tsmin = datum - time_fn
            tsmax = time_fn - datum
            #: First control from file name
            #: Make a ruf selection so we do not need to spend time to open to many files
            if ((tsmin.days*24*3600 + tsmin.seconds) > (timediff + 2 * 3600)) or \
                ((tsmax.days*24*3600 + tsmax.seconds) > (2 * timediff + 0.5 * 3600)):
                continue
            #: Read cloudsat file
            liddata = readFLXHRLidar(f)
            #: Second control from time stamps
            #: Lets lock inside the file and remove (continu) if there is no data to be used.
            if ((np.abs(liddata['sec1970'] - gpmdata['time_d'][0]*24*3600) <= timediff).sum() == 0) and \
                ((np.abs(liddata['sec1970'] - gpmdata['time_d'][1]*24*3600) <= timediff).sum() == 0):
                continue
            #: Append files to be used
            flxfiles.append(f)
            print('way')
        print(time.time() - tic)
        
        tic = time.time()
        #: Create KDTree form gpm lats/lons
        tree = spatial.KDTree(list(zip(gpmdata['lats'].ravel(), gpmdata['lons'].ravel())))
        print(time.time() - tic)
        sw = None
        tb = None
        for flxfile in flxfiles:
            #: Read the cloudsat data
            liddata = readFLXHRLidar(flxfile)
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
                else:
                    sw = np.concatenate((sw, liddata['QR'][0,timeInd,:]))
                    lw = np.concatenate((lw, liddata['QR'][1,timeInd,:]))
                    height = np.concatenate((height, liddata['Height'][timeInd,:]))
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
        tic = time.time()
        SW = np.zeros([heights.shape[0], tb.shape[0]])
        LW = np.zeros([heights.shape[0], tb.shape[0]])
        tb_new = np.zeros([tb.shape[0]])
        t = 0
        for i in range(len(tb)):
            if i == 0:
    #             tb_new = [tb[i]]
    #             tb_new[0] = tb[0]
                sumSW = np.zeros((heights.shape))
                numSW = np.zeros((heights.shape))
                sumLW = np.zeros((heights.shape))
                numLW = np.zeros((heights.shape))
            else:
                if gp[i] != gp[i-1]:
                    tb_new[t] = tb[i-1]
#                     tb_new.append(tb[i])
                    SW[:, t] = sumSW / numSW
                    LW[:, t] = sumLW / numLW
                
                    t = t + 1
                    
                    sumSW = np.zeros((heights.shape))
                    numSW = np.zeros((heights.shape))
                    sumLW = np.zeros((heights.shape))
                    numLW = np.zeros((heights.shape))
                    
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
    #             except:
    #                 pdb.set_trace()
        tb_new[t] = tb[i]
    #     tb_new.append(tb[i])
        SW[:, t] = sumSW / numSW
        LW[:, t] = sumLW / numLW
        SW = SW[:, 0:t+1]
        LW = LW[:, 0:t+1]
        tb_new = tb_new[0:t+1]
        np.save(tempname, [{'SW': SW, 'LW': LW, 'tb': tb, 'tb_new': tb_new, 'sw': sw, 'lw': lw}])
        print(time.time() - tic)
        print('save')
        pdb.set_trace()
    
    def loadTemps(tempnames):
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
        return retv['SW'], retv['LW'], retv['tb'], retv['tb_new'], retv['sw'], retv['lw']
    print('load')
    SW, LW, tb, tb_new, sw, lw = loadTemps(tempnames)
    
    avH = [] * height_lev
    lev = 24
    useInd = ~((tb == -9999) | (sw[:, lev] == -9.99) | (sw[:, lev] == 0.))
    print('Nr of SW cases = %d' %useInd.sum())
    
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    H, xedgesS, yedgesS, im = ax.hist2d(tb[useInd], sw[:,lev][useInd], bins=[50, 50])#, range= [[220, 300], [0., 0.05]])#, range=rang, bins=bins, norm=LogNorm(), vmin=1., vmax=1e4, cmap=cmap)
    ax.set_xlabel('GPM - Tb [K]')
    ax.set_ylabel('Cloudsat - SW [K/day]')
#     ax.set_title('Height = %d m, lev = %d' %(height[useInd, lev][0], lev))
    
    
    useInd = ~((tb == -9999) | (lw[:, lev] == -9.99) | (lw[:, lev] == 0.))
    print('Nr of LW cases = %d' %useInd.sum())
    ax = fig.add_subplot(1,2,2)
    H, xedgesL, yedgesL, im = ax.hist2d(tb[useInd], lw[:,lev][useInd], bins=[50, 50], range= [[220, 300], [-1., 1.]])#, range= [[220, 300], [0., 0.05]])#, range=rang, bins=bins, norm=LogNorm(), vmin=1., vmax=1e4, cmap=cmap)
    ax.set_xlabel('GPM - Tb [K]')
    ax.set_ylabel('Cloudsat - LW [K/day]')
#     ax.set_title('Height = %d m, lev = %d' %(height[useInd, lev][0], lev))
    
    figname = 'test'
    fig.savefig(figname + '.png')
    fig.show()
    
    SW_mean = np.zeros((heights.shape[0], xedgesS.shape[0]-1))
    LW_mean = np.zeros((heights.shape[0], xedgesL.shape[0]-1))
    for i in range(len(xedgesS)-1):
        use = (tb_new >= xedgesS[i]) & (tb_new < xedgesS[i+1])
        SW_mean[:,i] = np.nanmean(SW[:,use], axis=1)
    for i in range(len(xedgesL)-1):
        use = (tb_new >= xedgesL[i]) & (tb_new < xedgesL[i+1])
        LW_mean[:,i] = np.nanmean(LW[:,use], axis=1)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,2,1)
    im = ax.imshow(SW_mean, origin='lower', aspect=(SW_mean.shape[1]/SW_mean.shape[0]))#, aspect=0.75)
#     ax.set_title('Cloudsat - SW [K/day]')
    ax.set_xlabel('GPM - Tb [K]')
    ax.set_ylabel('Height [m]')
    cbar = fig.colorbar(im, orientation='horizontal')#, cax=cbar_ax, ticks=barticks)
    
    ax = fig.add_subplot(1,2,2)
    im = ax.imshow(LW_mean, origin='lower', aspect=(LW_mean.shape[1]/LW_mean.shape[0]))#, aspect=0.75)
    ax.set_xlabel('GPM - Tb [K]')
    ax.set_ylabel('Height [m]')
    cbar = fig.colorbar(im, orientation='horizontal')#, cax=cbar_ax, ticks=barticks)
#     ax.set_title('Height = %d m, lev = %d' %(height[useInd, lev][0], lev))
    fig.show()
    
    
    print(time.time() - tic_tot)
    pdb.set_trace()
#             minkded, minkdep = tree.query(pt)
    
    
    
    
    