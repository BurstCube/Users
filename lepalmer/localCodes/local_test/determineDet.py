#!/usr/bin/env python
# for triggerDict['Azimuth'][trigerDict['Energy']]:
#     for key = 'maxType':
#         if triggerDict[key] == 'counterXX':

import pylab as pl
import sys

import pickle 

f = open(sys.argv[1], 'rb')
mydict = pickle.load(f)
print(mydict)

AzimuthList = list(mydict.values())[0]
maxTypeList = list(mydict.values())[7]
print(AzimuthList)
AzimuthList = [float(i[0]) for i in AzimuthList]
GRB = {}
# print(AzimuthList)
maxDet1 = 0.0
minDet1 = 360.0
maxDet2 = 0.0
minDet2 = 360.0
maxDet3 = 0.0
minDet3 = 360.0
maxDet4 = 0.0
minDet4 = 360.0
count = -1
for i in maxTypeList:
    count += 1
    if i == 'counterXX':
        if AzimuthList[count] > maxDet1:
            maxDet1 == AzimuthList[count]
        if AzimuthList[count] < minDet1:
            minDet1 == AzimuthList[count]
    elif i == 'counterNXX':
        if AzimuthList[count] > maxDet2:
            maxDet2 == AzimuthList[count]
        if AzimuthList[count] < minDet2:
            minDet2 == AzimuthList[count]
    elif i == 'counterNXY':
        if AzimuthList[count] > maxDet3:
            maxDet3 == AzimuthList[count]
        if AzimuthList[count] < minDet3:
            minDet3 == AzimuthList[count]
    elif i == 'counterXY':
        if AzimuthList[count] > maxDet4:
            maxDet4 == AzimuthList[count]
        if AzimuthList[count] < minDet4:
            minDet4 == AzimuthList[count]

# energy = 200 ze = 15
count2 = -1
Azimuth = sys.argv[2]

for i in Azimuth:
    count2 += 1
    if sys.argv[2] >= minDet1 and angle <= maxDet1:
        GRB['detector'] = 'Det1'
    elif angle >= minDet2 and angle <= maxDet2:
        GRB['detector'] = 'Det2'
    elif angle >= minDet3 and angle <= maxDet3:
        GRB['detector'] = 'Det3'
    elif angle >= minDet4 and angle <= maxDet4:
        GRB['detector'] = 'Det4'
    GRB['Azimuth'][count2] = i


with open(sys.argv[3], 'wb') as f:
    pickle.dump(GRB, f, protocol=pickle.HIGHEST_PROTOCOL)

