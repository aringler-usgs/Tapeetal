#!/usr/bin/env python
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from obspy.core import UTCDateTime
import utils
import pickle
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
import matplotlib.font_manager
from obspy.clients.fdsn import Client
import sys


path = '/home/aringler/inst_test/magseis/noise_data/'
net ='*'
chan ='LHZ'
client = Client('IRIS')
stime = UTCDateTime('2019-032T03:30:00')
etime = UTCDateTime('2019-032T04:30:00')
stime2 = UTCDateTime('2019-032T07:30:00')
etime2 = UTCDateTime('2019-032T08:30:00')
permin = 40.
permax = 400.



def setupmap(central_lon, central_lat,handle):
    #handle = plt.axes(projection=ccrs.AlbersEqualArea(central_lon, central_lat))
    handle.set_extent(extent)

    handle.add_feature(cfeature.LAND)
    handle.add_feature(cfeature.OCEAN)
    handle.add_feature(cfeature.COASTLINE)
    handle.add_feature(cfeature.BORDERS, linestyle=':')
    handle.add_feature(cfeature.LAKES, alpha=0.5)
    handle.add_feature(cfeature.RIVERS)
    handle.add_feature(cfeature.STATES, edgecolor='gray')
    return handle

inv = utils.get_dataless(net, chan, stime, etime, client)


def make_noise(inv, stime, etime):
    lats, lons, noise = [], [], []
    for net in inv:
        if net.code not in ['AK', 'TA']:
            continue
        for sta in net:
            for chans in sta:
                print(chans)
                try:
                    st = client.get_waveforms(net.code, sta.code, chans.location_code,
                            chans.code, stime, etime, attach_response = True)
                    st.detrend('constant')
                    st.merge(fill_value=0.)
                    st.remove_response()
                    noise.append(st[0].std())
                    lats.append(sta.latitude)
                    lons.append(sta.longitude)
                except:
                    pass
                    print('No data')
    return lats, lons, noise




# lats, lons, noise = make_noise(inv, stime, etime)

# f = open('Test_data_second','wb')
# pickle.dump([lats, lons, noise], f)
# f.close()

# lats, lons, noise2 = make_noise(inv, UTCDateTime('2019-032T07:30:00'), UTCDateTime('2019-032T08:30:00'))

# f = open('Test_data_first','wb')
# pickle.dump([lats, lons, noise2], f)
# f.close()

f = open('Test_data_first','rb')
lats, lons, noise = pickle.load(f)

f = open('Test_data_second','rb')
lats2, lons2, noise2 = pickle.load(f)

minval = np.mean(noise) - 2*np.std(noise)
maxval = np.mean(noise) + 2*np.std(noise)

lats = np.array(lats)
lons = np.array(lons)
noise = np.array(noise)
noise2 = np.array(noise2)



idxG = []
for idx, pair in enumerate(zip(lats2, lons2)):
    lat, lon = pair[0], pair[1]
    if (lat not in lats) and (lon not in lons):
        idxG.append(idx)
        
lats2 = np.delete(lats2, idxG)
lons2 = np.delete(lons2, idxG)
noise2 = np.delete(noise2, idxG)

print(lats2)
print(lons2)
print(noise2)

noise = 20.*np.log10(noise)
noise2 = 20.*np.log10(noise2)

diff = noise2 - noise
print(diff)
print(noise)


fig= plt.figure(figsize=(12,16))


boxcoords=[min(lats) -1., min(lons)-1., max(lats) +1. , max(lons) + 1.]
print(boxcoords)
extent=[boxcoords[1], boxcoords[3], boxcoords[0], boxcoords[2]]
central_lon = np.mean(extent[:2])
central_lat = np.mean(extent[2:])            
ax = plt.subplot(3,2,1, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = noise, transform=ccrs.PlateCarree() )
cbar = fig.colorbar(sc, orientation='horizontal', shrink=0.7)
cbar.ax.set_xlabel('Power (dB)', fontsize=18)
plt.text(-0.1, 1., '(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
ax = plt.subplot(3,2,3, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = noise2, transform=ccrs.PlateCarree())
cbar = fig.colorbar(sc, orientation='horizontal', shrink=0.7)
cbar.ax.set_xlabel('Power (dB)', fontsize=18)
plt.text(-0.1, 1., '(c)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
ax = plt.subplot(3,2,5, projection=ccrs.AlbersEqualArea(central_lon, central_lat))
ax = setupmap(central_lon, central_lat, ax)
sc = ax.scatter(lons, lats, c = diff, transform=ccrs.PlateCarree())
cbar = fig.colorbar(sc, orientation='horizontal', shrink=0.7)
cbar.ax.set_xlabel('Difference (dB)', fontsize=18)
plt.text(-0.1, 1., '(e)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)

st = client.get_waveforms("IU", "COLA", "*",
                           "LF*", stime, etime, attach_response = True)

print(st)
st.detrend('constant')
st.merge(fill_value=0)
ax = plt.subplot(3,2,2)
for tr in st:
    t=tr.times()/(24*60*60)
    ax.plot(tr.times()/(24*60*60), np.sqrt((tr.data/(4.1943*10))**2), label=tr.id, alpha=0.5)
plt.ylim((0., 800.))
plt.ylabel('RMS (nT)')
plt.xlim((min(t), max(t)))
plt.xlabel('Time (hr)')
plt.text(-0.15, 1., '(b)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)

st = client.get_waveforms("IU", "COLA", "*",
                           "LF*", stime2, etime2, attach_response = True)

print(st)
st.detrend('constant')
st.merge(fill_value=0)
ax = plt.subplot(3,2,4)
for tr in st:
    t=tr.times()/(24*60*60)
    ax.plot(tr.times()/(24*60*60), np.sqrt((tr.data/(4.1943*10))**2), label=tr.id.replace('.', ' '), alpha=0.5)
plt.ylim((0., 800.))
plt.ylabel('RMS (nT)')
plt.xlim((min(t), max(t)))
plt.xlabel('Time (hr)')
plt.text(-0.15, 1., '(d)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, fontsize=18)
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.3),
          fancybox=True, shadow=True, ncol=1, fontsize=18)









plt.savefig('Alaskamap.png', format='PNG', dpi=200)
plt.savefig('Alaskamap.pdf', format='PDF', dpi=400)
plt.show()
