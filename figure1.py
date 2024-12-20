#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime, read
from obspy.core.inventory import read_inventory
import glob
import sys

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

from cartopy.feature import NaturalEarthFeature

# coast = NaturalEarthFeature(category='physical', scale='100m',
#                             facecolor=cfeature.COLORS['land'], name='coastline')

states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor=cfeature.COLORS['land'])

def setupmap(central_lon, central_lat,handle):
    #handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    handle.set_extent(extent)

    handle.add_feature(cfeature.LAND.with_scale('50m'), edgecolor='gray',facecolor=cfeature.COLORS['land'])
    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN.with_scale('50m'),facecolor=cfeature.COLORS['water'], edgecolor='gray'  )
    handle.add_feature(states_provinces, edgecolor='gray')
    handle.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor='gray',facecolor=cfeature.COLORS['land'] )
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)

    return handle
    
    
xmls = glob.glob('../magseis/*.xml')
for xml in xmls:
    if 'inv' not in vars():
        inv = read_inventory(xml)
    else:
        inv += read_inventory(xml)
        
latsTA, lonsTA, latsAK, lonsAK, latsNU, lonsNU = [] , [], [], [],[],[]
latsAdop, lonsAdop = [], []

AKupstas = ['SPIA', 'UNV', 'HOM','GAMB', 'COLD', 'MLY', 'CAPN', 'RC01','SSN', 'SKN', 'TRF', 'CHUM', 
        'CAST','DHY', 'KNK','KAI', 'MESA', 'PNL', 
        'PIN','MCAR', 'EYAK', 'BMR', 'PAX', 'NEA2', 'MLY', 'SCRK', 'RIDG','CUT','PPD', 'SCRK','SAW']
for net in inv:
    for sta in net:
        for chan in sta:
            try:
                coors = inv.get_coordinates(net.code + '.' + sta.code + '..LHZ')
                if sta.code in AKupstas:
                    lonsAdop.append(coors['longitude'])
                    latsAdop.append(coors['latitude'])
                elif net.code =='AK':
                    lonsAK.append(coors['longitude'])
                    latsAK.append(coors['latitude'])
                elif sta.code == 'EPYK':
                    lonsNU.append(coors['longitude'])
                    latsNU.append(coors['latitude'])
                elif (net.code == 'TA') and (sta.code[-1] =='K'):
                    lonsTA.append(coors['longitude'])
                    latsTA.append(coors['latitude'])
                else:
                    lonsNU.append(coors['longitude'])
                    latsNU.append(coors['latitude'])
            except:
                print('Problem: ' + sta.code)
                
                

boxcoords=[min(latsTA) -1., min(lonsTA)-1., max(latsTA) +1. , max(lonsTA) + 1.]
extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])            

fig= plt.figure(figsize=(14,14))
ax = plt.subplot(1,1,1, projection=ccrs.EquidistantConic(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
print('AK stations:' + str(len(latsAK)))
print('TA stations:' + str(len(latsTA)))
print('Upgraded stations:' + str(len(latsAdop)))
print('Canada stations:' + str(len(latsNU)))
sc = ax.plot(lonsAK, latsAK,markersize=10., transform=ccrs.Geodetic(), marker='o', linestyle= '', color='C1', label='AK Regional Network', zorder=50, markerfacecolor="None",markeredgewidth=3)
sc = ax.plot(lonsTA, latsTA,markersize=10., transform=ccrs.Geodetic(), marker='^', linestyle= '', color='C0', label='TA Alaska', zorder=50)

sc = ax.plot(lonsNU, latsNU,markersize=10., transform=ccrs.Geodetic(), marker='^', linestyle= '', color='C2', label='TA Alaska in Canada', zorder=50, markeredgewidth=3)
sc = ax.plot(lonsAdop, latsAdop,markersize=10., transform=ccrs.Geodetic(), marker='^', linestyle= '', color='C3', label='AK Upgraded Stations', zorder=50)

######## geomag  SIT DED CMO BRW          
glats = [57.0576,  64.8742, 71.3225]
glons = [-135.3273,   -147.8597, -156.6231]

sc = ax.plot(glons, glats, markersize=10.,transform=ccrs.Geodetic(), marker='s', color='C4',linestyle= '', label='USGS Geomag.', zorder=50)




fig.canvas.draw()
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
lats = [54, 57, 60, 63, 66, 69]
lons = [-170, -160, -150, -140, -130]
gl= ax.gridlines(draw_labels=True, xlocs=lons, ylocs=lats)
gl.xlabels_top = False
gl.ylabels_left = False
import matplotlib.ticker as mticker
gl.xlocator = mticker.FixedLocator(lons)
gl.ylocator = mticker.FixedLocator(lats)
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

plt.legend(ncol=5, loc=9, fontsize=14)
plt.savefig('Map_book.pdf', format='PDF', dpi=400)
plt.savefig('Map_book.jpg', format='JPEG', dpi=400)
plt.show()
