# -*- coding: utf-8 -*-
"""
@author: M. Dadighat
"""

import pandas as pd
import sqlite3
import numpy as np
import random
import time
import multiprocessing as mp
from operator import itemgetter
from math import ceil
from sqlite3 import Error
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astroplan import (Observer, AirmassConstraint,
                       AtNightConstraint, MoonSeparationConstraint)
from astroplan import is_observable
from astroplan import download_IERS_A


# Genetic algorithm to calculate observing schedules

# =============================================================================
#
# Constants
#
# =============================================================================
observatory_lat = 31.940556
observatory_long = -110.258056
observatory_alt = 1000  # in meters
observatory = Observer(longitude=observatory_long*u.deg,
                       latitude=observatory_lat*u.deg,
                       elevation=observatory_alt*u.m, name="Moka",
                       timezone="America/Phoenix")
overhead_per_target = 60  # in seconds
focus_delta = TimeDelta(180.0, format='sec')  # auto-focus time
observing_delta = TimeDelta(600.0, format='sec')
plan_date = Time(["2018-02-09"])
db_file_name = "targets.db"
acp_plan_header = "#CHILL -20.0\n#AUTOFOCUS\n\n"

initial_population_size = 50
maximum_airmass = 1.4
minimum_moon_sep = 20  # in degrees


# =============================================================================
#
# Utility functions
#
# =============================================================================
def getAstronomicalTwilightEvening(observatory, date):
    """ getAstronomicalTwilightEvening

    Uses the Astroplan module to calculate the time of astronomical twilight at 
    the beginning of the observing night - returns an astropy Time in UTC 
    formatted to include the date, hours, and minutes
    """
    astronomicalTwilight = observatory.twilight_evening_astronomical(
        date, which="next")
    atFormatted = Time(astronomicalTwilight, out_subfmt='date_hm')
    return atFormatted


def getAstronomicalTwilightMorning(observatory, date):
    """ getAstronomicalTwilightMorning

    Uses the Astroplan module to calculate the time of astronomical twilight at 
    the end of the observing night - returns an astropy Time in UTC 
    formatted to include the date, hours, and minutes
    """
    astronomicalTwilight = observatory.twilight_morning_astronomical(
        date, which="next")
    atFormatted = Time(astronomicalTwilight, out_subfmt='date_hm')
    return atFormatted


def get_airmass(target_coord, time):
    """getAirmass

    Calculates the airmass of a target at the given time - uses a globally 
    defined observer location. Time should be in astropy Time format
    """

    observatory = EarthLocation(lat=observatory_lat*u.deg,
                                lon=observatory_long*u.deg,
                                height=observatory_alt*u.m)

    target_alt_az = target_coord.transform_to(AltAz(obstime=time,
                                                    location=observatory))
    return target_alt_az.secz


def get_normalized_airmass(airmass):
    normalized_airmass = (airmass - maximum_airmass) / (1.0 - maximum_airmass)
    return normalized_airmass

# =============================================================================
#
# Fitness function
#
# =============================================================================


def evaluate_fitness(observing_schedule, total_time, total_targets):
    """ evaluate_fitness

    This evaluates the average fitness of a schedule taking into account:
        - total number of targets
        - total observing time (percentage of available time)
        - average target distance from moon
        - average target airmass

    It will be expanded to include preferred cadence per target and whether
    rising or setting targets are preferred (seasonally)

    Parameters
    ----------
    observing_schedule : list or `~pandas.DataFrame`...TBD
        Observing schedule of targets for a night

    total_time : int
        Total time available for observing for the night, in seconds

    total_targets : int
        Total number of targets that are theoretically observable in a night 

    Returns
    -------
    fitness_score : float
        The final fitness score for an observing plan
    """

    fitness_score = 0

    # Percentage weighting for each of the criteria
    time_weighting = .50
    count_weighting = .40
    airmass_weighting = .10
    moon_distance_weighting = 0

    # Scores will be out of 100 for each part individually

    # Time score
    time_score = 0
    observing_time = observing_schedule['total_interval'].sum()
    observing_time += focus_delta.value
    time_fraction = observing_time / total_time.value
    time_score = time_fraction * 100.0
    print('Observing time: ' + str(observing_time / 60.0) + ' min')
    print('Total time: ', str(total_time / 60.0) + ' min')
    print('Time fraction: ', str(time_fraction))
    print('Time score: ' + str(time_score))

    # Count score
    count_score = 0
    schedule_targets = len(observing_schedule.index)
    target_fraction = schedule_targets / float(total_targets)
    count_score = target_fraction * 100.0
    print('Schedule target count: ' + str(schedule_targets))
    print('Total possible targets: ', str(float(total_targets)))
    print('Count score: ' + str(count_score))

    # Airmass score
    # an airmass of 1.00 is perfect so it is 100%, 1.4 is the maximum acceptable
    # airmass so it is 0%
    airmass_score = 0
    airmass_interim_total = 0
    for value in observing_schedule['airmass']:
        airmass_interim_total += get_normalized_airmass(value)
    airmass_score = (airmass_interim_total /
                     float(len(observing_schedule.index))) * 100
    print('Average airmass score: ' + str(airmass_score))

    # Moon distance score
    moon_distance_score = 0

    fitness_score = (time_weighting * time_score
                     + count_weighting * count_score
                     + airmass_weighting * airmass_score
                     + moon_distance_weighting * moon_distance_score)
    print(fitness_score)
    return fitness_score

# =============================================================================
#
# Split
#
# =============================================================================


def split_schedule(schedule):
    number_of_targets = len(schedule.index)
    head_index = ceil(number_of_targets / 2.0)
    tail_index = number_of_targets - head_index

    head = schedule.head(head_index)
    tail = schedule.tail(tail_index)
#    print('Head row count: ' + str(len(head.index)))
#    print(head)
#    print('Tail row count: ' + str(len(tail.index)))
#    print(tail)
    return (head, tail)


# =============================================================================
#
# Refresh
#
# =============================================================================
def refresh_schedule(schedule):
    constraints = [AirmassConstraint(maximum_airmass),
                   AtNightConstraint.twilight_astronomical(),
                   MoonSeparationConstraint(min=minimum_moon_sep*u.deg)]
    dusk = getAstronomicalTwilightEvening(observatory, plan_date) + focus_delta
    #dawn = getAstronomicalTwilightMorning(observatory, plan_date)

    observing_slot_beginning = dusk
    observing_slot_end = observing_slot_beginning + observing_delta
    observing_slot = Time([observing_slot_beginning, observing_slot_end],
                          out_subfmt='date_hm')
    observable_at_new_spot = is_observable(constraints, observatory,
                                           schedule['coord'].tolist(),
                                           time_range=observing_slot)
    schedule['observable'] = pd.Series(observable_at_new_spot,
                                       index=schedule.index)
    schedule = schedule[schedule.observable]
    next_time = observing_slot_beginning
    next_delta = TimeDelta(0.0, format='sec')
    for index, row in schedule.iterrows():
        row['airmass'] = get_airmass(row['coord'], (next_time + next_delta))
        next_delta = TimeDelta(row['total_interval'], format='sec')
    return schedule


# =============================================================================
#
# Mutate
#
# =============================================================================
def mutate(schedule_to_mutate):
    """ mutate

    """
    print('Mutating...')
    # find the last item in the list
    #full_targets = data_frame.copy(deep=True)
    return False


# =============================================================================
#
# Reproduce
#
# =============================================================================
def reproduce(parents, total_time, total_targets):
    """ reproduce
    Uses a population of parent schedules to create a new generation of 
    potential schedules, as well as mutating them

    Parameters
    ----------
    parents : list
        List of tuples containing parent schedules and their fitness scores

    Returns
    -------
    children : list
        List tuples containing children schedules and their fitness scores
    """
    print('Reproducing...')
    children = []

    random.shuffle(parents)
    for index in range(0, (ceil(len(parents)/2))):
        parent_A = split_schedule(parents[index][0])
        parent_B = split_schedule(parents[index + 1][0])

        parent_A_head = parent_A[0]
        parent_B_head = parent_B[0]

#        print('Parent A head:')
#        print(parent_A_head)
#        print('Parent B head:')
#        print(parent_B_head)

        parent_A_tail = parent_A[1]
#        print('Parent A tail:')
#        print(parent_A_tail)
        parent_B_tail = parent_B[1]
#        print('Parent B tail:')
#        print(parent_B_tail)

        child_C = pd.concat([parent_A_head, parent_B_tail]
                            ).drop_duplicates().reset_index(drop=True)
        #parent_A_head.append(parent_B_tail, ignore_index=True)
#        child_C.reset_index(inplace=True, drop=True)
        child_D = parent_B_head.append(parent_A_tail, ignore_index=True)
        child_D.reset_index(inplace=True, drop=True)

        child_C = mutate(child_C)
        child_D = mutate(child_D)

        child_C = refresh_schedule(child_C)
        child_D = refresh_schedule(child_D)

#        print('Child C:')
#        print(child_C)

       # print(child_D)

        child_C = (child_C, evaluate_fitness(child_C, total_time,
                   total_targets))
        child_D = ((child_D, evaluate_fitness(child_D, total_time,
                   total_targets)))
        children.append(child_C)
        children.append(child_D)
    return children


# =============================================================================
#
# Get data
#
# =============================================================================
def get_all_targets():
    """ getAllTargets

    Gets all available target info from sqlite db - not all targets will 
    be used and additional sorting will be performed

    Parameters
    ----------
    none

    Returns
    -------
    target_df_fixed : list
        A list containing all targets from the database, with coordinate info
        sanitized and total exposure time per target calculated
    """

    try:
        connection = sqlite3.connect(db_file_name)
        target_df = pd.read_sql_query("SELECT * FROM targets", connection)
        exposure_df = df = pd.read_sql_query("SELECT * from exposures",
                                             connection)

        df = target_df.merge(exposure_df, how='outer')

        # Use target exposure's count, count, interval and overhead to
        # calculate total target time
        del df['id']
        df['total_interval'] = df['count']*df['interval']
        target_df_fixed = df.groupby(['target_id', 'name', 'ra', 'dec',
                                      'position_angle'],
                                     as_index=False)['total_interval'].sum()
        target_df_fixed['total_interval'] = target_df_fixed['total_interval']
        + overhead_per_target

        # Clean up dec string to make it compatible with SkyCoord constructor
        target_df_fixed['dec'] = target_df_fixed['dec'].str.replace(r'(\W\d\d)\s',
                                                                    '\\1d ')
        target_df_fixed['dec'] = target_df_fixed['dec'].str.replace("'", 'm')
        target_df_fixed['dec'] = target_df_fixed['dec'].str.replace('"', 's')

        # Replace RA/Dec columns with SkyCoord for use with constraints
        target_df_fixed['coord'] = SkyCoord(target_df_fixed['ra'],
                                            target_df_fixed['dec'])

        # TODO - delete ra/dec columns?

        return target_df_fixed
    except Error as e:
        print(e)
        return False


# =============================================================================
#
# Create schedule
#
# =============================================================================
def create_schedule(target_list, start_time, end_time):
    """ create_schedule

    Creates an inital schedule that fulfills all basic constraints both for the
    night and for each target. Targets added to the list are randomly selected
    from a pool of targets that are suitable at the given time

    Parameters
    ----------
    target_list : `~pandas.DataFrame`
        List of all targets possible during the night

    start_time : `~astropy.time.Time`
        Earliest possibile time for observations to begin - will be 
        astronomical twilight at dusk unless changed

    end_time : `~astropy.time.Time`
        Latest possibile time for observations to occur - will be 
        astronomical twilight at dawn unless changed


    Returns
    -------
    observing_schedule : list
        List of targets in observing schedule. Currently returning name, 
        airmass, and time - potentially need to change return format to be 
        more useful later depending on input needed to the reproduce function
    """

    print('Starting new schedule...')
    observing_schedule = pd.DataFrame(data=None, columns=target_list.columns)
    # create constraints - per target suitability
    current_constraints = [AirmassConstraint(maximum_airmass),
                           AtNightConstraint.twilight_astronomical(),
                           MoonSeparationConstraint(min=minimum_moon_sep*u.deg)]

    # start with start time, full list of targets
    observing_slot_beginning = start_time
    observing_slot_end = observing_slot_beginning + observing_delta
    observing_slot = Time([observing_slot_beginning, observing_slot_end],
                          out_subfmt='date_hm')
    target_pool = target_list.copy(deep=True)
    target_pool.reset_index(inplace=True, drop=True)
    print(target_pool)

    # repeat until list of available targets at the given time is empty or
    # none are available
    while not target_pool.empty:
        # apply constraints to list of targets
        observing_slot = Time([observing_slot_beginning, observing_slot_end],
                              out_subfmt='date_hm')
        observable = is_observable(current_constraints, observatory,
                                   target_pool['coord'].tolist(),
                                   time_range=observing_slot)
        target_pool['observable'] = observable
        target_pool['airmass'] = np.vectorize(get_airmass)(target_pool['coord'],
                                                           observing_slot_beginning)
        print(target_pool)

        # randomly pick target from list of possibilities
        possibilities = target_pool[target_pool.observable]
        print('Number of observable targets: '
              + str(len(possibilities.index)))
        if not possibilities.empty:
            # add target to list and remove from total targets available
            target_index = random.choice(possibilities.index.values)

            observing_schedule = observing_schedule.append(
                target_pool.loc[target_index])
            exposure_time = TimeDelta(target_pool['total_interval'][target_index],
                                      format='sec')
            observing_slot_beginning += exposure_time
            observing_slot_end = observing_slot_beginning + observing_delta
            target_pool.drop(target_pool.index[target_index], inplace=True)
            target_pool.reset_index(inplace=True, drop=True)

        else:
            break
    observing_schedule.reset_index(inplace=True, drop=True)
    return observing_schedule


def getExposureInfo(targetId, connection):
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


def print_acp_plan(plan_to_print):
    connection = sqlite3.connect(db_file_name)
    acpTxtFile = open("plan.txt", "w")
    sunset = getAstronomicalTwilightEvening(
        observatory, plan_date).to_datetime()
    sunrise = getAstronomicalTwilightMorning(
        observatory, plan_date).to_datetime()
    acpTxtFile.write("#WAITUNTIL 1, " +
                     sunset.strftime('%m-%d-%Y %H:%M') + '\n')
    acpTxtFile.write("#SHUTDOWNAT " +
                     sunrise.strftime('%m-%d-%Y %H:%M') + '\n')
    acpTxtFile.write(acp_plan_header)

    for index, target in plan_to_print.iterrows():
        print(target['target_id'])
        exposureInfo = getExposureInfo(target['target_id'], connection)

        acpTxtFile.write("#WAITAIRMASS 1.4, 5\n")
        acpTxtFile.write(str(exposureInfo))
        acpTxtFile.write("#POSANG " + str(target['position_angle']) + "\n")
        acpTxtFile.write(
            target['name'] + "\t" + target['coord'].to_string('hmsdms').replace(' ', '\t'))
        acpTxtFile.write('\n\n')

    acpTxtFile.close()


# =============================================================================
#
# Main algorithm
#
# =============================================================================

# 1. Get data from data source and remove any targets that can't be observed
#    at all during the night

if __name__ == '__main__':
    download_IERS_A()
    data_frame = get_all_targets()

    start_time = getAstronomicalTwilightEvening(
        observatory, plan_date) + focus_delta
    end_time = getAstronomicalTwilightMorning(observatory, plan_date)
    observing_range = Time([start_time, end_time], out_subfmt='date_hm')

    print("Time range is : " + str(observing_range.iso) + "\n\n")

    nightly_constraints = [AirmassConstraint(maximum_airmass),
                           AtNightConstraint.twilight_astronomical()]
    ever_observable = is_observable(nightly_constraints, observatory,
                                    data_frame['coord'].tolist(),
                                    time_range=observing_range)
    data_frame['observable'] = ever_observable
    data_frame = data_frame[data_frame.observable]

    print('Available targets for the night: ')
    print(data_frame['name'])

    # 2. Generate initial schedule set
    timer_start = time.time()

    schedule_population = []
    cpu_count = mp.cpu_count()
    print(str(cpu_count) + ' CPUs')
    pool = mp.Pool(processes=cpu_count)
    results = []
    for i in range(0, initial_population_size):
        #schedule = create_schedule(data_frame, start_time, end_time)
        # schedule_population.append(schedule)
        print(i)
        pool.apply_async(create_schedule, (data_frame, start_time, end_time,),
                         callback=results.append)
    for r in results:
        r.wait()
    pool.close()
    pool.join()

    total_night_time = TimeDelta(end_time - start_time, format='sec')

    total_target_count = len(data_frame)

    population = []

    print('\n\n')
    print('Initial schedules: ')
    for item in results:
        print('Initial schedule:')
        print('Number of targets: ' + str(len(item)))
        item_fitness = evaluate_fitness(item, total_night_time,
                                        total_target_count)
        print('Fitness score: ' + str(item_fitness))
        population.append((item, item_fitness))

        print('\n\n')
    print("%f seconds to calculate schedules" % (time.time() - timer_start))

    maximum_fitness = 70
    generation_maximum = max(population, key=itemgetter(1))[1]
    winner = max(population, key=itemgetter(1))[0]  # 0
    print(generation_maximum)
#    while generation_maximum < 70:
#
#        # 3. Reproduce and mutate
#        current_generation = population
#        current_generation.sort(key=itemgetter(1))
#        #print(current_generation)
#        current_generation = current_generation[:3]
#        for i in range(0, 1):
#            bonus_schedule = create_schedule(data_frame, start_time, end_time)
#            bonus_fitness = evaluate_fitness(bonus_schedule, total_night_time,
#                                        total_target_count)
#            current_generation.append((bonus_schedule, bonus_fitness))
#        next_generation = reproduce(current_generation, total_night_time,
#                                    total_target_count)
#
#        # 4. Evaluate
#        generation_maximum = max(next_generation, key=itemgetter(1))[1]
#        winner = max(next_generation, key=itemgetter(1))[0]
#        # 5. Next generation if needed
#

    print(winner)
    print_acp_plan(winner)
