#!/usr/bin/env python
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from scipy.optimize import fmin
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

# This code generates figure 1 of the Ringler et al. magnetic field paper

stime = UTCDateTime('2019-02-02T00:00:00')
etime = stime + 1*24*60*60

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times')
#mpl.rc('text', usetex=True)
mpl.rc('font',size=20)

chans = ['LHZ']
locs = ['10', '00']
locs =['10']
fig = plt.figure(2, figsize=(14,10))
for idx, loc in enumerate(locs):
    plt.subplot(3,1,idx+1)
    for chan in chans:

        net, sta, loc, chan = 'TA', 'POKR', '*', 'LHZ'
        #net, sta = 'IU', 'COLA'

        pmax = 1./800.
        pmin = 1./100.

        client = Client()
        st = client.get_waveforms(net, sta, loc,
                                   chan, stime-4*60*60, etime, attach_response = True)
        
        # When using some TA sites you need to chop off a sample
        for tr in st:
            tr.data = tr.data[1:]
        stR = st.copy()
        stR.detrend('constant')
        st += client.get_waveforms("IU", "COLA", "*",
                                   "LF*", stime, etime, attach_response = True)

        
        st.detrend('constant')
        st.merge(fill_value=0)
        #st2 = st.select(location='40').copy()
        for tr in st.select(location='40'):
            tr.data /= 41.943
            tr.data /= 10**9

        # Magnetic field data is in T
        inv = client.get_stations(network=net, station=sta, location=loc,
                            channel=chan, starttime=stime,endtime=etime, level="response")


        
        st.select(channel=chan).remove_response(inventory=inv, output='ACC')
        # Filter and then estimate coefficents
        st2 = st.copy()
        
        st.filter('bandpass', freqmin=pmax, freqmax=pmin)
        st2.filter('highpass', freq=1./40.)
        st.trim(stime,etime)
        st2.trim(stime,etime)
        st.taper(0.05)
        st2.taper(0.05)
        
        #st2.filter('bandpass', freqmin=pmax, freqmax=pmin)
       


        tr = st.select(location=loc, channel=chan)[0]
        times = tr.times()/(60*60)

        #plt.subplot(,1,1)
        #plt.plot([ 8 + 11./60., 8+ 11./60], [-100, 100], 'k', linewidth=2)
        #plt.plot(times, st2[0].data*10**9,label = (tr.id.replace('IU.COLA','')).replace('.',' ') + ' High-pass 40 s')
        #plt.text(-3.7, 18, '(b)')
        #plt.ylim((-18,18))
        #plt.legend(loc=4, ncol=2)
        #plt.xlim((min(times), max(times)))
        #plt.ylabel('Acc. $(nm/s^2)$')
        plt.subplot(2,1,2)
        plt.plot([ 9 + 11./60., 9+ 11./60], [-100, 100], 'k', linewidth=2)
        plt.plot(times, tr.data*10**9,label = (tr.id.replace('IU.COLA','')).replace('.',' ') + ' Band-pass 100 to 800 s')
        #plt.plot(times, correct(sol)*10**9,label = 'Magnetic Field Corrected', alpha=0.5)
        plt.xlim((min(times), max(times)))
        plt.ylim((-55., 55.))
        #plt.plot([ 8 + 11./60., 8+ 11./60], [-100, 100], 'k', linewidth=2)
        #plt.text(-3.7, 55, '(a)')
        plt.legend(loc=4, ncol=2)
        plt.ylabel('Acc. $(nm/s^2)$')




plt.subplot(2,1,1)
#st2.normalize()
# Here we plot the magnetic field data
for tr in st.select(location='40'):
    plt.plot([ 9 + 11./60., 9+ 11./60], [-200, 200], 'k', linewidth=2)
    plt.plot(times, tr.data*10**9,label = tr.id.replace('IU.COLA.40.',''), alpha=0.5)
    plt.xlim((min(times), max(times)))
    plt.legend(loc=4, ncol=3)
    plt.ylim((-140., 140.))
    #plt.text(-3.7, 140, '(b)')
    plt.ylabel('Magnetic Field ($nT$)')
plt.subplot(2,1,2)
plt.xlabel('Time February 2, 2019 (UTC Hr)')
plt.savefig('figure1_new_other' +  '.png', format='png', dpi=400)
#plt.show()