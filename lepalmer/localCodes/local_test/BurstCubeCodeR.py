#!/usr/bin/env python
from BurstCube.simGenerator import configurator
conf = configurator('config.yaml')

from BurstCube.bcSim import simFiles
sfs = simFiles('config.yaml')

# print(sfs.sims[0].simFile)
# type(sfs.sims[0].simFile)# fileN = 'FIN_1000.000keV_15.00ze_0.00az.inc1.id1.sim'

area = sfs.calculateAeff()

simulation = 0; detectOne = 0; detectTwo = 0; detectThree = 0; detectFour = 0
percentOne = 0; percentTwo = 0; percentThree = 0; percentFour = 0
# detectOneTot = 0
# detectTwoTot = 0
# detectThreeTot = 0
# detectFourTot = 0
import sys
detectTot = 0
triggerDict = {}
# f = open("dict.txt","w")
doubTriggerDict = {}
triggerSim = []
triggerSimNum = []
triggerSimCount = []
triggerSimEnergy = []
triggerSimList = []
counterSim = []
Det1 = []
Det2 = []
Det3 = []
Det4 = []
spercentOne = 0
spercentTwo = 0
spercentThree = 0
spercentFour = 0
simType = ''
# count = 0

import numpy as np
# triggerData = np.zeros(len(self.sims)), dtype={'names': ['keV', 'ze', 'az', 'counterXX', 'counterNXX', 'counterNXY', 'counterXY'], 'formats': ['float32', 'float32', 'float32','float32', 'float32', 'float32', 'float32']}

for file in sfs.sims:
    simulation = simulation + 1
#     count += 1
#     print(file.simFile)
    import re
    counterXX = 0
    counterNXX = 0
    counterNXY = 0
    counterXY = 0
    count = 0 
    d = 0
    prevLine = 'random'
    counter = 0

    with open(file.simFile) as sim:

#         print(file.simFile)
#         triggerDict = {'trig': [0], 'd' : 'sample'}

        for i, line in enumerate(sim):

            if line.startswith('ID'):
                ID = line.split()

                for i in range(len(ID)):
                    trigger = ID[1]

            if line.startswith('HTsim'):
                data = line.split(";")
                data = [w.strip(' ') for w in data]
#                 counter = counter + 1
        
                data[1] = float(data[1])
                data[2] = float(data[2])
#                 print(data[1])
#                 print("HELLO")
#                 print(data[2])
#                 print("GOODBYE")
#                 for i in range(len(data)):
#                 triggerDict.setdefault('triggers', [])
#                 triggerDict.setdefault('d', [])

                
                if (abs(5.52500 - data[1]) <= 0.05) and  (abs(5.52500 - data[2]) <= .05):
                    counterXX += 1
                    
                    d = 1
#                     triggerDict['trig'].append(counterXX)
#                     triggerDict['d'].append(d)
                    
                elif (abs(-5.52500 - data[1]) <= 0.05) and  (abs(5.52500 - data[2]) <= .05):
                    counterNXX += 1
                    d = 2
#                     triggerDict['trig'].append(counterNXX)
#                     triggerDict['d'].append(d)
                    
                elif (abs(-5.52500 - data[1]) <= 0.05) and  (abs(-5.52500 - data[2]) <= .05):
                    counterNXY += 1
                    d = 3
#                     triggerDict['trig'].append(counterNXY)
#                     triggerDict['d'].append(d)
                    
                elif (abs(5.52500 - data[1]) <= 0.05) and  (abs(-5.52500 - data[2]) <= .05):
                    counterXY += 1
                    d = 4
#                     triggerDict['trig'].append(counterXY)
#                     triggerDict['d'].append(d)
                    
#                 detectOne = counterXX + detectOne
#                 detectTwo = counterNXX + detectTwo
#                 detectThree = counterNXY + detectThree
#                 detectFour = counterXY + detectFour
#                 print(detectOne)
#                 print(detectTwo)
#                 print(detectThree)
#                 print(detectFour)
#                 print('*****')
                counter += 1   
                
            if prevLine.startswith('HTsim') and line.startswith('HTsim'):
                sfo = re.split('[_ . //]', file.simFile)
#                 doubTriggerDict['energy'] = sfo[4]
#                 doubTriggerDict['zenith'] = sfo[6]
#                 doubTriggerDict['azimuth'] = sfo[8]
#                 doubTriggerDict['trigNum'] = trigger
    
            prevLine = line
            if (counterXX > counterNXX) and (counterXX > counterNXY) and (counterXX > counterXY):
                triggerSimMax = counterXX
                triggerSimType = 'counterXX'
#                 triggerSimType = 'counterXX'
            elif (counterNXX > counterXX) and (counterNXX > counterNXY) and (counterNXX > counterXY):
                triggerSimMax = counterNXX
                triggerSimType = 'counterNXX'
#                 triggerSimType = 'counterNXX'
            elif (counterNXY > counterXX) and (counterNXY > counterNXX) and (counterNXY > counterXY):
                triggerSimMax = counterNXY
                triggerSimType = 'counterNXY'
#                 triggerSimType = 'counterNXY'

            elif (counterXY > counterXX) and (counterXY > counterNXX) and (counterXY > counterNXY):
                triggerSimMax = counterXY
                triggerSimType = 'counterXY'
#                 print(counterXX, '**')
#                 print(counterNXY, '***')
#                 triggerSimType = 'counterXY'
            count += 1

            sfo = re.split('[/ _]', file.simFile)
#     triggerSimList = [counterXX, counterNXX, counterNXY, counterXY]
#     print(triggerSimList)
#     triggerSimMax = max(triggerSimList)

    
#     for i in triggerSimList:
#         if triggerSimMax == triggerSimList[0]:
#             triggerSimType = 'counterNXX'
#         if triggerSimMax == triggerSimList[1]:
#             triggerSimType = 'counterNXX'
#         if triggerSimMax == triggerSimList[2]:
#             triggerSimType = 'counterNXY'
#         if triggerSimMax == triggerSimList[3]:
#             triggerSimType = 'counterXY'
                
    triggerSimCount.append(triggerSimType)
    triggerSimNum.append(triggerSimMax)
    triggerSim.append(triggerSimMax)
    triggerSim.append(triggerSimType)
#     triggerSimEnergy.append(sfo[12])
    counterSim.append(triggerSimType)
    spercentOne = counterXX / counter * 100
    spercentTwo = counterNXX / counter * 100        
    spercentThree = counterNXY / counter * 100
    spercentFour = counterXY / counter * 100
#     print(counter)
#     print(counter, spercentOne)
    Det1.append([spercentOne])
    Det2.append([spercentTwo])
    Det3.append([spercentThree])
    Det4.append([spercentFour])
    
    
#             print('***', trigger)

#             vars = re.split('[_ .]', file.simFile)
#             print("----------------------------------------------------------------------------")
#             print('Simulation', simulation, '--', trigger, ': Energy ({} keV), Azimuth ({} degrees), Zenith ({} degrees)'.format(vars[3], vars[5], vars[7]))
#             print("----------------------------------------------------------------------------")

#             print('Detector 1 has', counterXX , 'hits.')
#             print('Detector 2 has', counterNXX , 'hits.')
#             print('Detector 3 has', counterNXY , 'hits.')
#             print('Detector 4 has', counterXY , 'hits.')
        
#     print('**', count)
    detectOne = counterXX + detectOne
    detectTwo = counterNXX + detectTwo
    detectThree = counterNXY + detectThree
    detectFour = counterXY + detectFour
# print(len(triggerSimMax), '****')
# print(triggerSim)
detectTot = detectOne + detectTwo + detectThree + detectFour
# print(detectOne, detectTot)
percentOne = detectOne / detectTot * 100
percentTwo = detectTwo / detectTot * 100
percentThree = detectThree / detectTot * 100
percentFour = detectFour / detectTot * 100

'''Important, but not printed for saving purposes'''
print('**', percentOne, "**", percentTwo, '**', percentThree, "**", percentFour, "**")

import operator
'''Important, but not printed for saving purposes'''
print(sorted(enumerate(triggerSimNum), key=operator.itemgetter(1)))

Energy = list(area['keV'])
'''Important, but not printed for saving purposes'''

Azimuth = list(area['az'])
Zenith = list(area['ze'])

Zenith1 = []
for element in Zenith:
    Zenith1.append([element])
Zenith = Zenith1

Azimuth1 = []
for element in Azimuth:
    Azimuth1.append([element])
Azimuth = Azimuth1

Energy1 = []
for element in Energy:
    Energy1.append([element])
Energy = Energy1


triggerDict['Azimuth'] = Azimuth
triggerDict['Zenith'] = Zenith
triggerDict['Energy'] = Energy

c = 0
for i in range(len(triggerSimCount)):
    print('Index:', c, 'Value: ', triggerSimCount[i], 'Azimuth:', Azimuth[i], 'Det1:', Det1[i], 'Det2:', Det2[i], 'Det3:', Det3[i], 'Det4:', Det4[i])
    c = c + 1

triggerDict['Detector1'] = Det1
triggerDict['Detector2'] = Det2
triggerDict['Detector3'] = Det3
triggerDict['Detector4'] = Det4
triggerDict['maxType'] = counterSim
import pickle
with open(sys.argv[1], 'wb') as f:
    pickle.dump(triggerDict, f, protocol=pickle.HIGHEST_PROTOCOL)

# with open(sys.argv[1], 'rb') as f:
#    b = pickle.load(f)
# with open(sys.argv[1], 'w') as fp:
#   for key, value in triggerDict.items():
#      fp.write('%s:%s\n' % (key, value))


