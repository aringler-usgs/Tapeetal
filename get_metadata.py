#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy.core import UTCDateTime




stime = UTCDateTime('2020-001T00:00:00')
etime = stime + 1*24*60
client = Client('IRIS')
chan = 'LHZ'
netN = 'AK'


inv = client.get_stations(starttime=stime, endtime=etime, station="*",
                                  channel=chan, network=netN)

for net in inv:
    for sta in net:
        inv2 = client.get_stations(starttime=stime, endtime=etime, station=sta.code,
                                  channel=chan, network=netN, level="response")
        inv2.write(netN + '.' + sta.code + '.xml', format="STATIONXML")