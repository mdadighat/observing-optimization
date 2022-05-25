# -*- coding: utf-8 -*-
"""

@author: Michelle Dadighat

ACP Target List Checker/Organizer
"""

import re
import pandas as pd
import sqlite3
from sqlite3 import Error
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astroplan import Observer, FixedTarget, AirmassConstraint, AtNightConstraint, MoonSeparationConstraint
from astroplan import is_observable

################################################################
# Constants
################################################################
observatoryLat = 31.940556
observatoryLong = -110.258056
observatoryAlt = 1000 # in meters
utcOffset = -7 # hours off from  UTC, Mountain Standard Time
date = Time(["2018-01-05"])
imageOverhead = 60
focusTime = 180

acpPlanHeader = "#CHILL -20.0\n#AUTOFOCUS\n\n"

observatory = Observer(longitude=observatoryLong*u.deg, latitude=observatoryLat*u.deg,
                       elevation=observatoryAlt*u.m, name="Moka", 
                       timezone="America/Phoenix")



################################################################
# Target
# 
# A target extracted from the ACP plan - relevant directives 
# are stored separately, others are stored just in order to be
# read out later.
################################################################
class Target:
    name = ""
    ra = []
    dec = []
    pa = 0
    filters = []
    binning = []
    counts =[]
    intervals = []
    additionalKeywords = []
    
    def __init__(self, name, ra, dec, pa, filters, binning, counts, intervals, 
                 additionalKeywords):
        self.name = name
        self.ra = ra
        self.dec = dec
        self.pa = pa
        self.filters = filters
        self.binning = binning
        self.counts = counts
        self.intervals = intervals
    
class ObservingItem:
    name = ""
    coord = []
    airmass = 0
    integration = 0
    
    def __init__(self, name, coord, airmass, integration):
        self.name = name
        self.coord = coord
        self.airmass = airmass
        self.integration = integration

def getAstronomicalTwilightEvening(observatory, date):
    astronomicalTwilight = observatory.twilight_evening_astronomical(date, which="next")
    atFormatted = Time(astronomicalTwilight, out_subfmt='date_hm')
    return atFormatted

def getAstronomicalTwilightMorning(observatory, date):
    astronomicalTwilight = observatory.twilight_morning_astronomical(date, which="next")
    atFormatted = Time(astronomicalTwilight, out_subfmt='date_hm')
    return atFormatted

def createDatabaseConnection(dbFileName):
    try:
        connection =sqlite3.connect(dbFileName)
        return connection
    except Error as e:
        print(e)
        
    return None

def getAllTargets(connection):
    df = pd.read_sql_query("SELECT * FROM targets", connection)

    return df

def getExposuresForTarget(targetId, connection):
    df = pd.read_sql_query("SELECT * from exposures WHERE target_id=" + str(targetId), connection)
    
    return df

def getTargetId(targetName, connection):
    df = pd.read_sql_query("SELECT id from targets WHERE name='" + targetName + "'", connection)
    
    return df

def getPositionAngle(targetName, connection):
    df = pd.read_sql_query("SELECT position_angle from targets WHERE name='" + targetName + "'", connection)
    
    return df["position_angle"][0]

################################################################
# readACPFile(filename)
# param: File name of original ACP Plan
# returns: A list of targets with associated information
################################################################
def readACPFile(filename):
    targetList = [] # how to store targets?
    non_decimal = re.compile(r'[^\d.]+')

    
    with open(filename) as fp:
        name = ""
        ra = 0
        dec = 0
        pa = 0
        filters = []
        binning = []
        counts =[]
        intervals = []
        workingLines = []
        additionalKeywords = []
        
        for line in fp:
            #process line
            if(line[0] ==  '#' or line[0] == ';' or line[0] == '\n'):
               workingLines.append(line)
               
               if("#filter" in line.lower()):
                  filters = line[8:-1].split(',')
               elif("#binning" in line.lower()):
                    binning = line[9:-1].split(',')
               elif("#count" in line.lower()):
                    counts = line[7:-1].split(',')
               elif("#interval" in line.lower()):
                    intervals = line[10:-1].split(',')
               elif("#posang" in line.lower()):
                    pa = float(line[8:])
               else:
                   additionalKeywords.append(line)
            
            else:
                #we have found the first target
                workingLines.append(line)
                targetLine = line.split("\t")
                name = targetLine[0]
                raSplit = non_decimal.sub(' ', targetLine[1]).strip().split(' ')
                ra = map(float, raSplit)
                decSplit = non_decimal.sub(' ', targetLine[2][:-1]).strip().split(' ')
                dec = map(float, decSplit)
                newTarget = Target(name, ra, dec, pa, filters, binning, counts,
                                   intervals, additionalKeywords)
                targetList.append(newTarget)
                
    return targetList

def getAirmass(ra, dec, time):
    observatory = EarthLocation(lat=observatoryLat*u.deg, lon=observatoryLong*u.deg,
                                height=observatoryAlt*u.m)

    targetCoord = SkyCoord(ra=ra, dec=dec)

    targetCoordAltAz = targetCoord.transform_to(AltAz(obstime=time,location=observatory))  
    return targetCoordAltAz.secz


def isObservable(table, key_colnames):
    if table['is_observable'] == True:
        return True
    return False

def getBestTarget(targets, time):
    bestAirmass = 1.4
    bestTarget = 0
    
    for target in targets:
        airmass = getAirmass(target.ra, target.dec, time)
        if airmass < bestAirmass:
            bestAirmass = airmass
            bestTarget = target
    return (bestTarget,bestAirmass)

def recheckPlan(observingPlan, connection):
    observingTime = startTime
    for item in observingPlan:
        airmass = getAirmass(item.coord.ra, item.coord.dec, observingTime)
        if  airmass > 1.4:
            print("Plan failure - " + item.name + " " +str(airmass) + " at " + str(observingTime.iso))
            return (False, item)
        else:
            item.airmass = airmass
            targetid = getTargetId(item.name,connection)
            exposures = getExposuresForTarget(targetid['id'][0], connection)
            
            integrationTime = imageOverhead
            for x, exposure in exposures.iterrows():
                integrationTime += exposure['interval'] * exposure['count']
            time_delta = TimeDelta(integrationTime, format='sec')
            
            observingTime = Time(observingTime + time_delta)
        
    return (True, observingPlan)

def getExposureInfo(targetId, connection):
    df = pd.read_sql_query("SELECT * from exposures WHERE target_id=" + str(targetId), connection)
    
    filterString = "#FILTER "
    intervalString = "#INTERVAL "
    binningString = "#BINNING "
    countString = "#COUNT "
    for index, row in df.iterrows():

        filterString += (row["filter"] + ",")
        intervalString += (str(row["interval"]) + ",")
        binningString += (str(row["binning"]) + ",")
        countString += (str(row["count"]) + ",")

    filterString = filterString[:-1] 
    intervalString = intervalString[:-1] 
    binningString = binningString[:-1] 
    countString = countString[:-1]
    
    exposureInfo = filterString + '\n' + binningString + '\n' + countString + '\n' + intervalString + '\n'
    
    return exposureInfo

def printToACPPlan(observingPlan, connection):
    
    acpTxtFile = open("plan.txt", "w")
    sunset = getAstronomicalTwilightEvening(observatory, date).to_datetime()
    sunrise = getAstronomicalTwilightMorning(observatory, date).to_datetime()
    acpTxtFile.write("#WAITUNTIL 1, " + sunset.strftime('%m-%d-%Y %H:%M') + '\n')
    acpTxtFile.write("#SHUTDOWNAT " + sunrise.strftime('%m-%d-%Y %H:%M') + '\n')
    acpTxtFile.write(acpPlanHeader)
    
    for target in observingPlan:
        targetId = getTargetId(target.name, connection)['id'][0]
        exposureInfo = getExposureInfo(targetId, connection)
        
        acpTxtFile.write("#WAITAIRMASS 1.4, 5\n")
        acpTxtFile.write(str(exposureInfo))
        acpTxtFile.write("#POSANG " + str(getPositionAngle(target.name, connection)) + "\n")
        acpTxtFile.write(target.name + "\t" + target.coord.to_string('hmsdms').replace(' ', '\t'))
        acpTxtFile.write('\n\n')
    
    acpTxtFile.close()

#Read in plan targets
#readACPFile(sys.argv[1])
focus_delta = TimeDelta(focusTime, format='sec')
startTime = getAstronomicalTwilightEvening(observatory, date) + focus_delta
endTime = getAstronomicalTwilightMorning(observatory, date)
connection = createDatabaseConnection("targets.db")
targets = getAllTargets(connection)

targetTable = Table.from_pandas(targets)
targetTable.remove_column("position_angle")
targetTable.remove_column("id")

#fixedTargets = [FixedTarget(coord=SkyCoord(ra=ra*u.hour, dec=dec*u.deg), name=name)
#           for name, ra, dec in targetTable]

fixedTargets = [] 
for row in targetTable:
    decCoord  = re.sub(r"(\W\d\d)\s",r"\1d ", str(row["dec"] ))
    stringCoord = str(row["ra"]) + " " + decCoord.replace(' ', "")
    stringCoord = stringCoord.replace("'", 'm')
    stringCoord = stringCoord.replace('"', 's')
   
    fixedTarget = FixedTarget(coord=SkyCoord(stringCoord), name = row["name"])
    fixedTargets.append(fixedTarget)

time_range = Time([startTime, endTime])

#check overall target suitability for the night - airmass and twilight
constraints = [AirmassConstraint(1.4), AtNightConstraint.twilight_astronomical()]
ever_observable = is_observable(constraints, observatory, fixedTargets, time_range=time_range)

observabilityTable = Table()

observabilityTable['name'] = [target.name for target in fixedTargets]
observabilityTable['coord'] = [target.coord for target in fixedTargets]
observabilityTable['is_observable'] = ever_observable
#print(observabilityTable)

#reduce target list to those that are observable during the night
nightlyTargets = []
notObservable = []
obsTableGrouped = observabilityTable.group_by('is_observable')

for key, group in zip(obsTableGrouped.groups.keys, obsTableGrouped.groups):

    if key['is_observable'] == True:
        nightlyTargets = group
    else:
        notObservable = group

print("Nightly targets: ")  
for i in nightlyTargets:
    print(i['name']) 
print("\n\n")         
print("Not observable: ")
for i in notObservable:
    print(i['name'])
observingList = []

for target in fixedTargets:
    if target.name in nightlyTargets['name']:
        observingList.append(target)

#sort remaining targets by visibility at current observing time

observingTime = startTime
currentConstraints = [AirmassConstraint(1.4), MoonSeparationConstraint(min=20*u.deg)]
observing_delta = TimeDelta(600.0, format='sec')
time_delta = TimeDelta(600.0, format='sec')
end_time = Time(observingTime+time_delta)

observing_range=Time([observingTime.iso, end_time.iso])

observingPlan = []

while (len(observingList) > 0):
    
    print("...")
    currently_observable = is_observable(currentConstraints, observatory, 
                                         observingList, time_range=observing_range)
    
    #print(observing_range)
    currentTable = Table()

    currentTable['name'] = [target.name for target in observingList]
    currentTable['coord'] = [target.coord for target in observingList]
    currentTable['is_observable'] = currently_observable

    #possibleTargetsTable = currentTable.group_by('is_observable')#.groups[1]
    #filtered = possibleTargetsTable.groups.filter(isObservable)

    usableTargets = []
    for row in currentTable:
        if row['is_observable']:
           # print(row['name'], row['is_observable'])
            fTarget = FixedTarget(coord=row['coord'], name=row['name'])
            usableTargets.append(fTarget)
        else:
            continue
    #for target in usableTargets:
      #  print(target)
    
    #print(observingTime)
    if usableTargets:
        bestTarget =  getBestTarget(usableTargets, observingTime)
        
        #print(bestTarget[0])
        
        for item in observingList:
            if item.name == bestTarget[0].name:
                 observingList.remove(item)
        
        #increase observing check time by exposure length times a factor to account
        #for slewing, plate solving, etc.
        id = getTargetId(bestTarget[0].name,connection)
        exposures = getExposuresForTarget(id['id'][0], connection)
        integrationTime = imageOverhead
        for x, exposure in exposures.iterrows():
            integrationTime += exposure['interval'] * exposure['count']
        
        print("integration time: " + str(integrationTime))
        time_delta = TimeDelta(integrationTime, format='sec')
        print(integrationTime, id['id'][0])
        
        observingTime = Time(observingTime + time_delta)
        end_time = Time(observingTime + observing_delta)
        
        observing_range=Time([observingTime.iso, end_time.iso])
        
        planItem = ObservingItem(bestTarget[0].name, bestTarget[0].coord,
                                 bestTarget[1], integrationTime)
        observingPlan.append(planItem) 
    else:
        print('Current list:')
        timeTest = startTime
        for planItem in observingPlan:
            item_time_delta = TimeDelta(planItem.integration, format='sec')
            timeTest = timeTest + item_time_delta
            print(planItem.name, planItem.airmass, str(timeTest.iso))
        print('\nUnusable targets for round 1:')
        for target in observingList:
            print(target.name)
        break

#recheck list
        
checkTime = 0
fixed = []
print("\n")
copyPlan = observingPlan
for recheck_item in observingList:
    print('rechecking ' + recheck_item.name)
    
    for plan_item in observingPlan:

        delta_check = TimeDelta(plan_item.integration, format="sec")
        if not checkTime:
            checkTime = startTime + delta_check
        else:
            checkTime = checkTime + delta_check
            
        print("Checking at " + str(checkTime.iso))
        itemAirmass = getAirmass(recheck_item.coord.ra, recheck_item.coord.dec, checkTime)
        
        if (itemAirmass < 1.4 and itemAirmass > 1):
            print("works after " + plan_item.name +" "+ str(itemAirmass))
            index = observingPlan.index(plan_item)
            
            id = getTargetId(bestTarget[0].name,connection)
            exposures = getExposuresForTarget(id['id'][0], connection)
            itemIntegration = imageOverhead
            for x, exposure in exposures.iterrows():
                itemIntegration += exposure['interval'] * exposure['count']
            observingItem = ObservingItem(recheck_item.name, recheck_item.coord,
                                 itemAirmass, itemIntegration)
            copyPlan.insert(index+1, observingItem)
            workable = recheckPlan(copyPlan, connection)
            if not workable[0]:
                
                print("doesn't work after " + plan_item.name + " - breaks plan")
                print("attempting reshuffle")
                #take item that doesn't work, move it back a space, and recheck again
                reshuffleItem = workable[1]
                reshuffleIndex = copyPlan.index(reshuffleItem)
                success = False

                
                for idx in range(0,(reshuffleIndex - 1)):
                    copyPlan[reshuffleIndex-idx], copyPlan[reshuffleIndex-idx-1] = copyPlan[reshuffleIndex-idx-1], copyPlan[reshuffleIndex-idx]
                    reshuffledWorkable = recheckPlan(copyPlan, connection)
                    if not reshuffledWorkable[0]:
                        copyPlan[reshuffleIndex-idx-1], copyPlan[reshuffleIndex-idx] = copyPlan[reshuffleIndex-idx], copyPlan[reshuffleIndex-idx-1]
                    else:
                        success = True
                        break
                if not success:
                    copyPlan.remove(observingItem)
                else:
                    fixed.append(recheck_item)
                    print("on to next item...")
                    break
                #keep going through the observing plan
            else:
                fixed.append(recheck_item)
                print("on to next item...")
                
                break

        else:
            print("doesn't work after " + plan_item.name + " " +str(itemAirmass))

    print('\n\n')
    checkTime = 0
    
obsTime = startTime
finalPlan = recheckPlan(copyPlan, connection)[1]
print("Final Results:\n\n")
for line in finalPlan:
    delta = TimeDelta(line.integration, format="sec")
    print(delta)
    print(line.name, line.airmass, str(obsTime.iso))
    obsTime = obsTime + delta
    
for line in observingList:
    if not any(item.name == line.name for item in copyPlan):
        print("failed to place: " + line.name)

printToACPPlan(finalPlan, connection)















