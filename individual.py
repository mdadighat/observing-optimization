# -*- coding: utf-8 -*-
"""

@author: mmpkb8
"""
import pandas as pd
import sqlite3
from sqlite3 import Error
import random
from astroplan import download_IERS_A
from astropy.coordinates import SkyCoord
from astroplan import (Observer, AirmassConstraint,
                       AtNightConstraint)
from astropy import units as u
from astropy.time import Time, TimeDelta
from astroplan import is_observable


db_file_name = "targets.db"
overhead_per_target = 60  # in seconds
initial_population_size = 20
observatory_lat = 31.940556
observatory_long = -110.258056
observatory_alt = 1000  # in meters
observatory = Observer(longitude=observatory_long*u.deg,
                       latitude=observatory_lat*u.deg,
                       elevation=observatory_alt*u.m, name="Moka",
                       timezone="America/Phoenix")
focus_delta = TimeDelta(180.0, format='sec')  # auto-focus time
maximum_airmass = 1.4
plan_date = Time(["2018-01-13"])


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


start_time = getAstronomicalTwilightEvening(
    observatory, plan_date) + focus_delta
end_time = getAstronomicalTwilightMorning(observatory, plan_date)
night_time_range = Time([start_time, end_time], out_subfmt='date_hm')

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


def confirm_all_at_night(schedule):
    valid_schedule = False
   # interim_schedule = schedule.copy(deep=True)
   # airmass_constraint = [AirmassConstraint(maximum_airmass)]
   # observable = is_observable(airmass_constraint, observatory,
   #                                 data_frame['coord'].tolist(),
   #                                 time_range=observing_range)
   # interim_schedule['observable'] = observable
   # interim_schedule = interim_schedule[interim_schedule.observable]
   # valid_schedule =  not any(interim_schedule.observable == False)

   # start at dusk
   # start_time

   # confirm observation is at night

    return valid_schedule


def confirm_maximum_airmass(schedule):
    return False


def confirm_minimum_moon_distance(schedule):
    return False


def calculate_observing_time(schedule):
    return False


def calculate_violations(chromosome):
    violations = []
    # 1 - everything occurs at night
    violations.append(confirm_all_at_night(chromosome.genotype))

    # 2 - everything is less than 1.4 airmass
    violations.append(confirm_maximum_airmass(chromosome.genotype))

    # 3 - everything is minimum distance from moon
    violations.append(confirm_minimum_moon_distance(chromosome.genotype))

    return violations

# holds the schedule information


class Chromosome:
    def __init__(self, genotype):
        self.genotype = genotype


class Observation:
    def __init__(self, name, coord, exposure_time):
        self.name = name
        self.coord = coord
        self.exposure_time = exposure_time


class Individual:
    def __init__(self, chromosome_seed):
        self.chromosome = generate_chromosome(chromosome_seed)
        self.ranking = 0
        self.crowded_distance = 0
        self.target_count = len(self.chromosome.genotype)
        self.observing_time = calculate_observing_time(self.chromosome)
        self.violations = calculate_violations(self.chromosome)

    def __lt__(self, other):
        return crowded_comparison(self, other) == False

    def __gt__(self, other):
        return crowded_comparison(self, other) == True
   # def __eq__(self, other):
    #    return crowded_comparison(self, other) == True

    def __le__(self, other):
        return crowded_comparison(self, other) == False

    def __ge__(self, other):
        return crowded_comparison(self, other) == True

    def __ne__(self, other):
        return crowded_comparison(self, other) == False

    def __key(self):
        return (self.chromosome.genotype)

    def __eq__(x, y):
        return x.__key() == y.__key()


def crowded_comparison(individual_A, individual_B):

    rank_A = individual_A.ranking
    rank_B = individual_B.ranking
    distance_A = individual_A.crowded_distance
    distance_B = individual_B.crowded_distance

    if (rank_A < rank_B) or ((rank_A == rank_B) and (distance_A > distance_B)):
        return True
    return False

# max parameter is the maximum number of targets that are possible to schedule -
# not supporting non-seqential multiple observations of a single target yet


def generate_chromosome(seed_values):
    genotype = []
    number_of_possibilities = len(seed_values.index)
    genotype_length = random.randint(0, number_of_possibilities)
    print('Chromosome length: ' + str(genotype_length))
    for i in range(0, genotype_length):
        name = seed_values.sample(1)['name']
        coord = seed_values.sample(1)['coord']
        exposure_time = seed_values.sample(1)['total_interval']
        genotype.append(Observation(name, coord, exposure_time))
    # wtf is a valid chromosome for a schedule? how much to limit

    chromosome = Chromosome(genotype)
    return chromosome


def generate_parent_population(chromosome_seed):
    parent_population = []

    for i in range(0, initial_population_size):
        individual = Individual(chromosome_seed)
        parent_population.append(individual)

    return parent_population


class ProblemDefinition:

    name = 'Optimal schedule'

    # total number of targets
    def f1(self, individual):
        return individual.target_count

    # total observing time
    def f2(self, individual):
        return individual.observing_time

    # airmass (avg.)
    def f3(self, individual):
        return False


if __name__ == '__main__':
    download_IERS_A()
    data_frame = get_all_targets()
    print('Initial data frame:')
    print(data_frame)

    start_time = getAstronomicalTwilightEvening(
        observatory, plan_date) + focus_delta
    end_time = getAstronomicalTwilightMorning(observatory, plan_date)
    observing_range = Time([start_time, end_time], out_subfmt='date_hm')

    print("Time range is : " + str(observing_range.iso) + "\n\n")

    nightly_constraints = [AirmassConstraint(maximum_airmass),
                           AtNightConstraint.twilight_astronomical()]
    ever_observable = is_observable(nightly_constraints, observatory,
                                    data_frame['coord'].tolist(),
                                    time_range=night_time_range)
    data_frame['observable'] = ever_observable
    data_frame = data_frame[data_frame.observable]

    population = generate_parent_population(data_frame)
    # for thingy in population:
    # print(thingy.chromosome)
    # 1 print('\n\n')
