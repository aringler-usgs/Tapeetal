#!/usr/bin/env python
import glob
import sys
import pickle
from PIL import Image
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime
from obspy.signal.filter import envelope
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

debug = True
#stime = UTCDateTime('2019-032T02:29:49')
#etime = UTCDateTime('2019-032T09:28:17')
#stime = UTCDateTime('2019-299T02:46:41')
#etime = UTCDateTime('2019-299T16:24:07')
stime = UTCDateTime('2019-02-28T00:00:00')
etime = UTCDateTime('2019-02-28T23:59:59')
f = open('20190228_photo.csv','r')
timescolor, r, g, b = [],[],[],[]
for line in f:
    line = line.replace(' ','')
    line = line.split(',')
    timescolor.append(float(line[0]))
    r.append(float(line[1]))
    g.append(float(line[2]))
    b.append(float(line[3]))

r = np.array(r)
g = np.array(g)
b = np.array(b)
timescolor = np.array(timescolor)
idx = [(timescolor > 5.) & (timescolor < 10)]

maxval = max([max(r[idx]), max(g[idx]), max(b[idx])])

r /= maxval
g /= maxval
b /= maxval

f.close()

#g =np.convolve(g, np.ones(10), 'same') / 10  
#r =np.convolve(r, np.ones(10), 'same') / 10 
#b =np.convolve(b, np.ones(10), 'same') / 10 
#3.8667 and 16.2167


net, sta, loc, chan = 'TA', 'POKR', '01', 'LHZ'
pmax = 1./300.
pmin = 1./100.

client = Client()
st = client.get_waveforms(net, sta, loc,
                            chan, stime-60*60, etime, attach_response = True)

st.detrend('constant')
st.merge(fill_value=0)
st.filter('bandpass', freqmin=pmax, freqmax=pmin, zerophase=True)
st.remove_response(output='ACC')
st.trim(stime,etime)


times2 = st[0].times()/(60*60) 
print(times2)
#times += stime.hour + stime.minute/(60.) + stime.second/(60.*60.)



times = (np.arange(len(r))/len(r))*(24) 


fig = plt.figure(1, figsize=(12,12))
plt.subplot(2,1,1)
data = st[0].data*10**9
#data /= max(data)

data_envelope = envelope(st[0].data*10**9)
plt.subplot(2,1,2)
plt.fill_between(times2, 0, data_envelope/max(data_envelope), label='Seismic', color='0.7')
plt.ylim((0., 1.2))
plt.plot(times, r, color='r', label='Red')
plt.plot(times, g, color='g', label='Green')
plt.plot(times, b, color = 'b', label='Blue')
plt.legend(loc=9, ncol=4)
plt.xlabel('Time (hr)')
plt.ylabel('Amplitude (Normalized)')
plt.ylim((0,1))
#3.8667 and 16.2167
plt.xlim((3.8867, 16.2167))
plt.text(2, 1, '(b)')
plt.subplot(2,1,1)
plt.plot(times2, data, label=st[0].id)

plt.subplot(2,1,1)
plt.xlabel('Time(hr)')
plt.ylabel('Acceleration ($nm/s^2$)')
plt.xlim((3.8867, 16.2167))
plt.text(2, 4, '(a)')
plt.legend(loc=9)
plt.savefig('Intensity_' + str(stime.month).zfill(2) + str(stime.day).zfill(2) + '_' + st[0].id + 'NEW.png', format='PNG')
plt.savefig('Intensity_' + str(stime.month).zfill(2) + str(stime.day).zfill(2) + '_' + st[0].id + 'NEW.pdf', format='PDF')
plt.show()