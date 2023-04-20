# -*- coding: utf-8 -*-
"""
@author: Michelle Dadighat

ACP Target List Checker/Organizer

This script is designed to generate an ACP plan file of the best targets 
to fill a given night.  You can either set it to read in an existing file and update
it for the current date, or read in all the desired observing targets from the DB and
create a plan from that. 

Optimal targets are defined by lowest airmass, although a moon separation constraint is
also used when determining if a target is considered observable for the night.
 
Current configuration - reading in targets from DB to create a full plan. 
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
observatoryAlt = 1000  # in meters
utcOffset = -7  # hours off from  UTC, Mountain Standard Time

# Date for plan to be created - hardcoded for now since this isn't 
# called from anywhere else
date = Time(["2022-01-05"])

# in seconds
imageOverhead = 60  # average amount of time to take an exposure not included elsewhere
focusTime = 180 # on average, < 3.5 FWHM

recheckAll = True

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
    counts = []
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


################################################################
# Observing Item
#
# similar to a target, but with just the basic information needed 
# to see where the target fits in the plan
################################################################
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
    """ getAstronomicalTwilightEvening
    Gets the earliest possible time for the targets to be observable

    Parameters
    ----------
    observatory : astroplan.Observer
        Location at which to determine twilight
    date : astropy.Time
        Date for which to determine twilight

    Returns
    -------
    atFormatted : astropy.Time
        astronomical twilight formatted yyyy-mm-dd hh:mm
    """
    astronomicalTwilight = observatory.twilight_evening_astronomical(
        date, which="next")
    atFormatted = Time(astronomicalTwilight,
                       format='iso', out_subfmt='date_hm')
    return atFormatted


def getAstronomicalTwilightMorning(observatory, date):
    """ getAstronomicalTwilightMorning
    Gets the latest possible time for the targets to be observable

    Parameters
    ----------
    observatory : astroplan.Observer
        Location at which to determine twilight
    date : astropy.Time
        Date for which to determine twilight

    Returns
    -------
    atFormatted : astropy.Time
        astronomical twilight formatted yyyy-mm-dd hh:mm
    """
    astronomicalTwilight = observatory.twilight_morning_astronomical(
        date, which="next")
    atFormatted = Time(astronomicalTwilight,
                       format='iso', out_subfmt='date_hm')
    return atFormatted


def createDatabaseConnection(dbFileName):
    """ createDatabaseConnection
    Connects to sqlite database of all desired targets, even those out of season

    Parameters
    ----------
    dbFileName : string
        Name/path of database to connect

    Returns
    -------
    connection : sqlite3.Connection
        Connection object used to load targets
    """
    try:
        connection = sqlite3.connect(dbFileName)
        return connection
    except Error as e:
        print(e)

    return None


def getAllTargets(connection):
    """ getAllTargets
    Returns a list of all targets from the database

    Parameters
    ----------
    connection : sqlite3.Connection
        Connection object used to load targets

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing all targets from the database
    """
    df = pd.read_sql_query("SELECT * FROM targets", connection)

    return df


def getExposuresForTarget(targetId, connection):
    """ getExposuresForTarget
    Returns a list of all desired exposures for the given target (time, filters, etc)

    Parameters
    ----------
    targetId : int
        Unique target id from the database

    connection : sqlite3.Connection
        Connection object used to load targets

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing all exposures needed in a plan instance for a given target
    """
    df = pd.read_sql_query(
        "SELECT * from exposures WHERE target_id=" + str(targetId), connection)

    return df


def getTargetId(targetName, connection):
    """ getTargetId
    Returns unique id for target of a given name

    Parameters
    ----------
    targetName : string
        Name of target for which to obtain id

    connection : sqlite3.Connection
        Connection object used to load targets

    Returns
    -------
    targetId : int
        Unique target id from the database for a given target
    """
    df = pd.read_sql_query(
        "SELECT target_id from targets WHERE name='" + targetName + "'", connection)
    
    targetId = df['target_id'][0]

    return targetId


def getPositionAngle(targetName, connection):
    """ getPositionAngle
    Returns the camera position angle for the given target

    Parameters
    ----------
    targetName : string
        Name of target for which to obtain position angle

    connection : sqlite3.Connection
        Connection object used to load targets

    Returns
    -------
    posAngle : int
        Unique target id from the database for a given target
    """
    df = pd.read_sql_query(
        "SELECT position_angle from targets WHERE name='" + targetName + "'", connection)
    posAngle = df["position_angle"][0]
    return posAngle


def readACPFile(filename):
    """ readACPFile
    A list of targets from a given ACP plan file with associated information

    Parameters
    ----------
    filename : string
        File name of original ACP Plan

    Returns
    -------
    targetList : list[Target]
        all targets with exposure info as read in from the file
    """
    targetList = []  # how to store targets?
    non_decimal = re.compile(r'[^\d.]+')

    with open(filename) as fp:
        name = ""
        ra = 0
        dec = 0
        pa = 0
        filters = []
        binning = []
        counts = []
        intervals = []
        workingLines = []
        additionalKeywords = []

        for line in fp:
            # process line
            if(line[0] == '#' or line[0] == ';' or line[0] == '\n'):
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
                # we have found the first target
                workingLines.append(line)
                targetLine = line.split("\t")
                name = targetLine[0]
                raSplit = non_decimal.sub(
                    ' ', targetLine[1]).strip().split(' ')
                ra = map(float, raSplit)
                decSplit = non_decimal.sub(
                    ' ', targetLine[2][:-1]).strip().split(' ')
                dec = map(float, decSplit)
                newTarget = Target(name, ra, dec, pa, filters, binning, counts,
                                   intervals, additionalKeywords)
                targetList.append(newTarget)

    return targetList


def getAirmass(ra, dec, time):
    """ getAirmass
    Get the airmass for a given location at a given time

    Parameters
    ----------
    ra : string
        target right ascension
    dec : string
        target declination
    time : string
        time to get airmass for

    Returns
    -------
    airmass : float
        airmass of the target
    """
    observatory = EarthLocation(lat=observatoryLat*u.deg, lon=observatoryLong*u.deg,
                                height=observatoryAlt*u.m)

    targetCoord = SkyCoord(ra=ra, dec=dec)

    targetCoordAltAz = targetCoord.transform_to(
        AltAz(obstime=time, location=observatory))
    return targetCoordAltAz.secz


def getBestTarget(targets, time):
    """ getBestTarget
    Finds the best target in the list for the given time, as determined by airmass

    Parameters
    ----------
    targets : list[Target]
        all observable targets
    time : astropy.Time
        time to check for the target
        

    Returns
    -------
    bestTarget : tuple(astroplan.FixedTarget, astropy.Time)
        ideal target at the given time
    """
    bestAirmass = 1.4
    bestTarget = 0

    for target in targets:
        airmass = getAirmass(target.ra, target.dec, time)
        if airmass < bestAirmass:
            bestAirmass = airmass
            bestTarget = target
    return (bestTarget, bestAirmass)


def recheckPlan(observingPlan, connection):
    """ recheckPlan
    Uses the total exposure time for each target to determine if 
    all targets meet the minimum airmass requirement where they were
    placed in the plan.

    Parameters
    ----------
    observingPlan : list[Target]
        all observable targets
    connection : sqlite3.Connection
        connection to observing database with targets and exposure info
        

    Returns
    -------
    checked : tuple(boolean, list[ObservingItem])
        Boolean is True if the plan is valid, False if a target in the 
        plan ends up with an airmass > 1.4 at the time it was placed.
        In the case of a failed plan, the ObservingItem is the target 
        that needs to be moved, otherwise it is just the original plan
    """
    observingTime = startTime
    for item in observingPlan:
        airmass = getAirmass(item.coord.ra, item.coord.dec, observingTime)
        if airmass > 1.4:
            #print("Plan failure - " + item.name + " " +
            #      str(airmass) + " at " + str(observingTime.iso))
            return (False, item)
        else:
            item.airmass = airmass
            targetid = getTargetId(item.name, connection)
            exposures = getExposuresForTarget(targetid, connection)

            integrationTime = imageOverhead
            for x, exposure in exposures.iterrows():
                integrationTime += exposure['interval'] * exposure['count']
            time_delta = TimeDelta(integrationTime, format='sec')

            observingTime = Time(observingTime + time_delta)

    return (True, observingPlan)


def getExposureInfo(targetId, connection):
    """ getExposureInfo
    Finds all exposure info for a given target (filters, count,
    intervals and binning)

    Parameters
    ----------
    targetId : int
        database id for the target for which to find exposure info
    connection : sqlite3.Connection
        connection to observing database with targets and exposure info
        

    Returns
    -------
    exposureInfo : string
        Formatted string of exposure info for a given target ready to print
        to an ACP Plan file
    """
    df = pd.read_sql_query(
        "SELECT * from exposures WHERE target_id=" + str(targetId), connection)

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

    exposureInfo = filterString + '\n' + binningString + \
        '\n' + countString + '\n' + intervalString + '\n'

    return exposureInfo


def printToACPPlan(observingPlan, connection):
    """ printToACPPlan
    Prints the final observing plan to a text file in ACP Plan format

    Parameters
    ----------
    observingPlan : list[Target]
        the final list of targets to observe
    connection : sqlite3.Connection
        connection to observing database with targets and exposure info
        

    Returns
    -------
    none
    """

    acpTxtFile = open("plan.txt", "w")
    sunset = getAstronomicalTwilightEvening(observatory, date).to_datetime()
    sunrise = getAstronomicalTwilightMorning(observatory, date).to_datetime()
    acpTxtFile.write("#WAITUNTIL 1, " +
                     sunset.strftime('%m-%d-%Y %H:%M') + '\n')
    acpTxtFile.write("#SHUTDOWNAT " +
                     sunrise.strftime('%m-%d-%Y %H:%M') + '\n')
    acpTxtFile.write(acpPlanHeader)

    for target in observingPlan:
        targetId = getTargetId(target.name, connection)
        exposureInfo = getExposureInfo(targetId, connection)

        acpTxtFile.write("#WAITAIRMASS 1.4, 5\n")
        acpTxtFile.write(str(exposureInfo))
        acpTxtFile.write(
            "#POSANG " + str(getPositionAngle(target.name, connection)) + "\n")
        acpTxtFile.write(target.name + "\t" +
                         target.coord.to_string('hmsdms').replace(' ', '\t'))
        acpTxtFile.write('\n\n')

    acpTxtFile.close()


###############################
# Fixing/Creating the plan
###############################

# Read in plan targets
# Option 1 - Read in from current ACP plan, with intent to modify for current date
# readACPFile(sys.argv[1])

#Option 2 - Read targets in from DB, recreate plan from scratch
focus_delta = TimeDelta(focusTime, format='sec')
startTime = getAstronomicalTwilightEvening(observatory, date) + focus_delta
endTime = getAstronomicalTwilightMorning(observatory, date)
connection = createDatabaseConnection("targets.db")
targets = getAllTargets(connection)

targetTable = Table.from_pandas(targets)
targetTable.remove_column("position_angle")
#targetTable.remove_column("id")

# fixedTargets = [FixedTarget(coord=SkyCoord(ra=ra*u.hour, dec=dec*u.deg), name=name)
#           for name, ra, dec in targetTable]

fixedTargets = []
for row in targetTable:
    decCoord = re.sub(r"(\W\d\d)\s", r"\1d ", str(row["dec"]))
    stringCoord = str(row["ra"]) + " " + decCoord.replace(' ', "")
    stringCoord = stringCoord.replace("'", 'm')
    stringCoord = stringCoord.replace('"', 's')

    fixedTarget = FixedTarget(coord=SkyCoord(stringCoord), name=row["name"])
    fixedTargets.append(fixedTarget)

time_range = Time([startTime, endTime])

# check overall target suitability for the night - airmass and twilight
constraints = [AirmassConstraint(
    1.4), AtNightConstraint.twilight_astronomical()]
ever_observable = is_observable(
    constraints, observatory, fixedTargets, time_range=time_range)

observabilityTable = Table()

observabilityTable['name'] = [target.name for target in fixedTargets]
observabilityTable['coord'] = [target.coord for target in fixedTargets]
observabilityTable['is_observable'] = ever_observable
# print(observabilityTable)

# reduce target list to those that are observable during the night
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

# sort remaining targets by visibility at current observing time

observingTime = startTime
currentConstraints = [AirmassConstraint(
    1.4), MoonSeparationConstraint(min=20*u.deg)]
observing_delta = TimeDelta(600.0, format='sec')
time_delta = TimeDelta(600.0, format='sec')
end_time = Time(observingTime+time_delta)

observing_range = Time([observingTime.iso, end_time.iso])

observingPlan = []

# go through the list of targets and add them to the plan according to 
# best airmass at the time they would be observed, this accounts for the
# different exposure times for each target
while (len(observingList) > 0):

    print("...")
    currently_observable = is_observable(currentConstraints, observatory,
                                         observingList, time_range=observing_range)

    # print(observing_range)
    currentTable = Table()

    currentTable['name'] = [target.name for target in observingList]
    currentTable['coord'] = [target.coord for target in observingList]
    currentTable['is_observable'] = currently_observable

    usableTargets = []
    for row in currentTable:
        if row['is_observable']:
           # print(row['name'], row['is_observable'])
            fTarget = FixedTarget(coord=row['coord'], name=row['name'])
            usableTargets.append(fTarget)
        else:
            continue

    if usableTargets:
        bestTarget = getBestTarget(usableTargets, observingTime)

        #for item in observingList:
        #    if item.name == bestTarget[0].name:
        #        observingList.remove(item)
        observingList = [x for x in observingList if x.name is not bestTarget[0].name]

        # increase observing check time by exposure length times a factor to account
        # for slewing, plate solving, etc.
        id = getTargetId(bestTarget[0].name, connection)
        exposures = getExposuresForTarget(id, connection)
        integrationTime = imageOverhead
        for x, exposure in exposures.iterrows():
            integrationTime += exposure['interval'] * exposure['count']

        intTime = "integration time: " + str(integrationTime)
        time_delta = TimeDelta(integrationTime, format='sec')
        print(intTime, id)

        observingTime = Time(observingTime + time_delta)
        end_time = Time(observingTime + observing_delta)

        observing_range = Time([observingTime.iso, end_time.iso])

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

# Recheck list of unused targets, optional
# This takes a while to run so only use if plan is not urgently needed
copyPlan = observingPlan
if (recheckAll):
    checkTime = 0
    fixed = []
    print("\n")
    for recheck_item in observingList:
        print('rechecking ' + recheck_item.name)

        for plan_item in observingPlan:

            delta_check = TimeDelta(plan_item.integration, format="sec")
            if not checkTime:
                checkTime = startTime + delta_check
            else:
                checkTime = checkTime + delta_check

            print("Checking at " + str(checkTime.iso))
            itemAirmass = getAirmass(
                recheck_item.coord.ra, recheck_item.coord.dec, checkTime)

            if (itemAirmass < 1.4 and itemAirmass > 1):
                print(recheck_item.name + " possible after " + plan_item.name + " " + str(itemAirmass))
                index = observingPlan.index(plan_item)

                id = getTargetId(bestTarget[0].name, connection)
                exposures = getExposuresForTarget(id, connection)
                itemIntegration = imageOverhead
                for x, exposure in exposures.iterrows():
                    itemIntegration += exposure['interval'] * exposure['count']
                observingItem = ObservingItem(recheck_item.name, recheck_item.coord,
                                            itemAirmass, itemIntegration)
                copyPlan.insert(index+1, observingItem)

                #verify new item doesn't break the rest of the plan
                workable = recheckPlan(copyPlan, connection)
                if not workable[0]:

                    print(observingItem.name + " doesn't work after " + plan_item.name + " - breaks plan")
                    print("attempting reshuffle")
                    # take item that doesn't work, move it back a space, and recheck again
                    reshuffleItem = workable[1]
                    reshuffleIndex = copyPlan.index(reshuffleItem)
                    success = False

                    for idx in range(0, (reshuffleIndex - 1)):
                        copyPlan[reshuffleIndex-idx], copyPlan[reshuffleIndex-idx -
                                                            1] = copyPlan[reshuffleIndex-idx-1], copyPlan[reshuffleIndex-idx]
                        reshuffledWorkable = recheckPlan(copyPlan, connection)
                        if not reshuffledWorkable[0]:
                            copyPlan[reshuffleIndex-idx-1], copyPlan[reshuffleIndex -
                                                                    idx] = copyPlan[reshuffleIndex-idx], copyPlan[reshuffleIndex-idx-1]
                        else:
                            success = True
                            break
                    if not success:
                        for item in copyPlan:
                            if item.name == observingItem.name:
                                copyPlan.remove(item)
                                print("Removing " + item.name)
                        #copyPlan.remove(observingItem)
                    else:
                        fixed.append(recheck_item)
                        print("on to next item...")
                        break
                    # keep going through the observing plan
                else:
                    fixed.append(recheck_item)
                    print("on to next item...")

                    break

            else:
                print("doesn't work after " + plan_item.name + " " + str(itemAirmass))

        print('\n\n')
        checkTime = 0

# Final results for review
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
