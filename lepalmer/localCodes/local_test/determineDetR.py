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
GRB = {'Az': [], 'Detector': []}
# print(AzimuthList)
count2 = -1
# Az = [29.3, 128.99, 201.546, 300]
Az = sys.argv[2]
Detector = []
Detect1 = 'Detect1'
Detect2 = 'Detect2'
Detect3 ='Detect3'
Detect4 = 'Detect4'
maxDetect1 = 0.0
minDetect1 = 360.0
maxDetect2 = 0.0
minDetect2 = 360.0
maxDetect3 = 0.0
minDetect3 = 360.0
maxDetect4 = 0.0
minDetect4 = 360.0
count = -1

for i in maxTypeList:
    count += 1
#     print(AzimuthList[count])
#     print(i)
    if i == 'counterXX':
        if AzimuthList[count] > maxDetect1 and AzimuthList[count] != 0 and AzimuthList[count] != 360:
            maxDetect1 = float(AzimuthList[count])
        if AzimuthList[count] < minDetect1:
            minDetect1 = float(AzimuthList[count])
    elif i == 'counterNXX':
        if AzimuthList[count] > maxDetect2:
            maxDetect2 = float(AzimuthList[count])
        if AzimuthList[count] < minDetect2:
            minDetect2 = float(AzimuthList[count])
    elif i == 'counterNXY':
        if AzimuthList[count] > maxDetect3:
            maxDetect3 = float(AzimuthList[count])
        if AzimuthList[count] < minDetect3:
            minDetect3 = float(AzimuthList[count])
    elif i == 'counterXY':
        if AzimuthList[count] > maxDetect4:
            maxDetect4 = float(AzimuthList[count])
        if AzimuthList[count] < minDetect4:
            minDetect4 = float(AzimuthList[count])
            
GRB['max1'] = maxDetect1
GRB['max2'] = maxDetect2
GRB['max3'] = maxDetect3
GRB['max4'] = maxDetect4

GRB['min1'] = minDetect1
GRB['min2'] = minDetect2
GRB['min3'] = minDetect3
GRB['min4'] = minDetect4

# energy = 200 ze = 15
for angle in Az:
    count2 += 1
    if angle >= minDetect1 and angle <= maxDetect1:            
        Detector.append(Detect1)
    elif angle >= minDetect2 and angle <= maxDetect2:            
        Detector.append(Detect2)
# GRB['detector'] = 'Det2'
    elif angle >= minDetect3 and angle <= maxDetect3:
        Detector.append(Detect3)
#         GRB['detector'] = 'Det3'
    elif angle >= minDetect4 and angle <= maxDetect4:
        Detector.append(Detect4)
#         GRB['detector'] = 'Det4'
#     GRB['Azimuth'] = angle
GRB['Az'].append(Az)
GRB['Detector'].append(Detector)

print("**", minDetect1, minDetect2, minDetect3, minDetect4, "**")
print("**", maxDetect1, maxDetect2, maxDetect3, maxDetect4, "**")

print('**', GRB, "**")
with open('sys.argv[3]', 'wb') as f:
    pickle.dump(GRB, f, protocol=pickle.HIGHEST_PROTOCOL)
