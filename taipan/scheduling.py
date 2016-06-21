import numpy as np
import datetime
from collections import OrderedDict

import os
import sys
import time
import logging
import pickle

import ephem

# ______________________________________________________________________________
#
# CONSTANTS

# define twilight values for when observing can start/must end

SOLAR_TWILIGHT_HORIZON = np.radians(-10.)      # radians everywhere
LUNAR_TWILIGHT_HORIZON = np.radians(0.)        # radians everywhere

SOLAR_HORIZON = SOLAR_TWILIGHT_HORIZON - np.radians(0.24)
LUNAR_HORIZON = LUNAR_TWILIGHT_HORIZON - np.radians(0.27)

SUN = ephem.Sun()
MOON = ephem.Moon()

# define and establish emph.Observer object for UKST-Taipan
#  all from https://en.wikipedia.org/wiki/Siding_Spring_Observatory
UKST_LATITUDE = -31.272231
UKST_LONGITUDE = +149.071233 
UKST_ELEVATION = 1165  # metres

UKST_TELESCOPE = ephem.Observer()
UKST_TELESCOPE.lat = np.radians(UKST_LATITUDE)   # radians everywhere
UKST_TELESCOPE.lon = np.radians(UKST_LONGITUDE)  # radians everywhere
UKST_TELESCOPE.elevation = UKST_ELEVATION
UKST_TELESCOPE.temp = 10.

SECONDS_PER_DAY = 86400.

ALMANAC_RESOLUTION_MAX = 60. * 4.


# ______________________________________________________________________________
# CLASS DEFINITIONS
# ______________________________________________________________________________


class Almanac(object):
    """
    Object which stores observability information for a specific point on the
    sky.
    """

    # Define class attributes
    _ra = None
    _dec = None
    _start_date = None
    _end_date = None
    _airmass = None
    _resolution = None
    _minimum_airmass = None
    _observer = None

    # Setters & getters
    @property
    def ra(self):
        """
        RA for this almanac
        """
        return self._ra

    @ra.setter
    def ra(self, r):
        if r is None:
            self._ra = None
            return
        r = float(r)
        if r < 0. or r >= 360.:
            raise ValueError('Almanac must have 0 <= RA < 360')
        self._ra = r

    @property
    def dec(self):
        """
        RA for this almanac
        """
        return self._dec

    @dec.setter
    def dec(self, d):
        if d is None:
            self._dec = None
            return
        d = float(d)
        if d < -90. or d > 90.:
            raise ValueError('Almanac must have -90 <= RA <= 90')
        self._dec = d

    @property
    def start_date(self):
        """
        Almanac start date as a datetime date object
        """
        return self._start_date

    @start_date.setter
    def start_date(self, d):
        if not isinstance(d, datetime.date):
            raise ValueError("start_date must be an instance of "
                             "datetime.datetime.date")
        if self.end_date is not None:
            if d >= self.end_date:
                raise ValueError("Requested start_date is after the existing "
                                 "end_date for this almanac")
        self._start_date = d

    @property
    def end_date(self):
        """
        Almanac start date as a datetime date object
        """
        return self._end_date

    @end_date.setter
    def end_date(self, d):
        if not isinstance(d, datetime.date):
            raise ValueError("end_date must be an instance of "
                             "datetime.datetime.date")
        if self.end_date is not None:
            if d >= self.end_date:
                raise ValueError("Requested end_date is after the existing "
                                 "start_date for this almanac")
        self._end_date = d

    @property
    def airmass(self):
        """
        An OrderedDict of the type:
        date_j2000: airmass
        An OrderedDict is used to keep the entries in date (i.e. key) order
        """
        return self._airmass

    @airmass.setter
    def airmass(self, a):
        if not isinstance(a, dict):
            raise ValueError("airmass must be a dictionary")
        self._airmass = a

    @property
    def minimum_airmass(self):
        """
        The minimum airmass that this almanac will consider 'observable'
        """
        return self._minimum_airmass

    @minimum_airmass.setter
    def minimum_airmass(self, a):
        if a is None:
            self._minimum_airmass = None
            return
        a = float(a)
        if a < 0:
            raise ValueError('Minimum airmass must be > 0')
        if a > 100:
            raise ValueError('Minimum airmass must be < 100 '
                             '(max. practical value is ~ 38; 99 is used '
                             'as a special value by this module')
        self._minimum_airmass = a

    @property
    def resolution(self):
        """
        Resolution of almanac in minutes
        """
        return self._resolution

    @resolution.setter
    def resolution(self, r):
        r = float(r)
        if r <= 0:
            raise ValueError('Resolution must be > 0 mins')
        if r < float(ALMANAC_RESOLUTION_MAX):
            raise ValueError('Resolution must be < %d mins' %
                             (ALMANAC_RESOLUTION_MAX, ))

    @property
    def observer(self):
        """
        ephem observer object
        """
        return self._observer

    @observer.setter
    def observer(self, o):
        if not isinstance(o, ephem.Observer):
            raise ValueError('observer must be an instance of the Observer '
                             'class from pyephem')
        self._observer = o

    # Helper functions
    def compute_end_date(self, observing_period):
        """
        Compute the end date for this almanac and set it
        """
        delta = datetime.timedelta(observing_period)
        end_date = self.start_date + delta
        self.end_date = end_date

    def get_observing_period(self):
        return (self.end_date -
                self.start_date).total_seconds() / SECONDS_PER_DAY

    def generate_file_name(self):
        filename = 'almanac_R%3.1f_D%2.1f_start%s_end%s_res%3f.amn' % (
            self.ra, self.dec, self.start_date.strftime('%y%m%d'),
            self.end_date.strftime('%y%m%d'), self.resolution,
        )
        return filename

    # Initialization
    def __init__(self, ra, dec, start_date, end_date=None,
                 observing_period=None, observer=UKST_TELESCOPE,
                 minimum_airmass=2.0, resolution=15.,
                 populate=True):
        if end_date is None and observing_period is None:
            raise ValueError("Must specify either one of end_date or "
                             "observing_period")

        self.ra = ra
        self.dec = dec
        self.start_date = start_date
        if end_date:
            self.end_date = end_date
        else:
            self.compute_end_date(self, observing_period)
        self.airmass = {}
        self.observer = observer
        self.minimum_airmass = minimum_airmass
        self.resolution = resolution

        # See if an almanac with these properties already exists in the PWD
        if populate:
            load_success = self.load()
            if not load_success:
                # Couldn't find an almanac to load, so we need to
                # populate this one
                self.calculate_airmass()

        return

    # Save & read from disk
    # Uses pickle to seralize objects
    def save(self, filename=None, filepath='./'):
        if filepath[-1] != '/':
            raise ValueError('filepath must end with /')
        if filename is None:
            filename = self.generate_file_name()
        with open('%s%s' % (filepath, filename, )) as fileobj:
            pickle.dump(self, fileobj)
        return

    def load(self, filepath='./'):
        """
        This function will attempt to load an almanac based on the file name of
        almanacs available in the directory specified by filepath. This will be
        based on the RA, Dec, starttime, endtime and resolution information
        encoded within the filename. If successful, the function returns True.
        If a filename matching the calling
        almanac is not found, then the function will exit and return False.
        """
        if filepath[-1] != '/':
            raise ValueError('filepath must end with /')
        files_present = os.listdir(filepath)
        if self.generate_file_name() not in files_present:
            return False

        with open('%s%s' % (filepath, filename,)) as fileobj:
            file_almanac = pickle.load(self, fileobj)

        self.ra = file_almanac.ra
        self.dec = file_almanac.dec
        self.start_date = file_almanac.start_date
        self.end_date = file_almanac.end_date
        self.observer = file_almanac.observer
        self.resolution = file_almanac.resolution
        self.airmass = file_almanac.airmass
        self.minimum_airmass = file_almanac.airmass

        return True

    # Computation functions
    def generate_almanac_bruteforce(self, full_output=False):
        # Define almanac-dependent horizons
        observable = np.arcsin(1. / self.minimum_airmass)

        # initialise ephem fixed body object for pointing
        target = ephem.FixedBody()
        target._ra = np.radians(self.ra)  # radians everywhere
        target._dec = np.radians(self.dec)  # radians everywhere

        # Calculate the time grid
        if self.airmass is None or len(self.airmass) == 0:
            self.observer.date = self.start_date
            observing_period = (
                self.end_date - self.start_date
                               ).total_seconds() / SECONDS_PER_DAY
            dates_j2000 = (self.observer.date +
                           np.arange(0, observing_period,
                                     self.resolution / 1440.))
        else:
            dates_j2000 = sorted(self.airmass.keys())

        # Perform the computation
        time_total = 0.
        sol_alt, lun_alt, target_alt, dark_time = [], [], [], []
        for d in dates_j2000:
            self.observer.date = d
            SUN.compute(self.observer)
            MOON.compute(self.observer)
            target.compute(self.observer)

            can_observe = (SUN.alt <
                           SOLAR_HORIZON) and (MOON.alt < LUNAR_HORIZON)
            if full_output:
                sol_alt.append(SUN.alt)
                lun_alt.append(MOON.alt)
                target_alt.append(target.alt)
                dark_time.append(can_observe)
            else:
                if can_observe and (target.alt > observable):
                    time_total += self.resolution

        if full_output:
            # convert results to numpy arrays
            sol_alt = np.array(sol_alt)
            lun_alt = np.array(lun_alt)
            target_alt = np.array(target_alt)
            dark_time = np.array(dark_time)
            return dates_j2000, sol_alt, lun_alt, target_alt, dark_time
        else:
            return time_total / 60.

    def calculate_airmass(self):
        logging.debug('Computing airmass values for almanac at %3.1f, %2.1f '
                      'from %s to %s' %
                      (self.ra, self.dec, self.start_date.strftime('%y-%m-%d'),
                       self.end_date.strftime('%y-%m-%d'),
                       ))
        dates, sun, moon, target, dark_time = self.generate_almanac_bruteforce(
            full_output=True)

        airmass_values = np.clip(np.where(target > np.radians(10.),
                                          1./np.sin(target), 99.), 0., 9.)

        if self.airmass is None:
            self.airmass = {}

        for i in range(len(dates)):
            self.airmass[dates[i]] = airmass_values[i]

        return


class DarkAlmanac(Almanac):
    """
    Subclass of Almanac, which holds information on what times are 'dark'
    Holds no RA, Dec information
    """

    # Setters and getters
    @property
    def ra(self):
        return self._ra

    @ra.setter
    def ra(self, r):
        r = float(r)
        if abs(r - 0.) > 1e-5:
            raise ValueError('Dark almanac must have RA = 0')

    @property
    def dec(self):
        return self._dec

    @dec.setter
    def dec(self, r):
        r = float(r)
        if abs(r - 0.) > 1e-5:
            raise ValueError('Dark almanac must have DEC = 0')

    @property
    def minimum_airmass(self):
        return self._minimum_airmass

    @minimum_airmass.setter
    def minimum_airmass(self, a):
        if a is not None:
            raise ValueError('Dark almanac does not require a minimum_airmass'
                             ' value')

    # Create a new class attribute, dark_time, which simply aliases
    # airmass
    # Furthermore, we'll override the Almanac getter & setter for airmass
    # so that it doesn't appear as a valid class attribute
    # However, the hidden value _airmass will still be used
    @property
    def airmass(self):
        raise AttributeError('Dark almanacs do not have an airmass catalogue -'
                             ' use dark_time instead')

    @airmass.setter
    def airmass(self, a):
        # Fail silent
        return
        # if a is not None or a != {}:
        #     raise AttributeError('Dark almanacs do not have an airmass '
        #                          'catalogue - use dark_time instead')

    @property
    def dark_time(self):
        """
        A catalogue of whether time should be considered 'dark' or not
        """
        return self._airmass

    @dark_time.setter
    def dark_time(self, d):
        if not isinstance(d, dict):
            raise ValueError('dark_time must be a dictionary of datetimes '
                             'and corresponding Boolean values')
        self._airmass = d

    # Class functions
    def generate_file_name(self):
        filename = 'darkalmanac_start%s_end%s_res%3f.amn' % (
            self.start_date.strftime('%y%m%d'),
            self.end_date.strftime('%y%m%d'), self.resolution,
        )
        return filename

    # Initialization
    def __init__(self, start_date, end_date=None,
                 observing_period=None, observer=UKST_TELESCOPE,
                 resolution=15., populate=True):
        # 'super' the Almanac __init__ method, but do NOT attempt to
        # populate the DarkAlmanac from this method (uses a different
        # method)
        logging.debug('super-ing standard almanac creation')
        super(DarkAlmanac, self).__init__(0., 0., start_date, end_date=end_date,
                                          observing_period=observing_period,
                                          observer=observer,
                                          minimum_airmass=None,
                                          resolution=resolution, populate=False)

        logging.debug('Looking for existing dark almanac file')
        # See if an almanac with these properties already exists in the PWD
        if populate:
            load_success = self.load()
            if not load_success:
                # Couldn't find an almanac to load, so we need to
                # populate this one
                self.create_dark_almanac()

    def create_dark_almanac(self):
        """
        Perform the necessary calculations to populate this almanac
        """
        logging.debug('Creating dark almanac')
        logging.debug('Calculating from standard almanac functions...')
        dates, sun, moon, target, is_dark_time = \
            self.generate_almanac_bruteforce(full_output=True)

        times_per_day = 1440. / self.resolution
        if abs(times_per_day - int(times_per_day)) > 1e-5:
            raise ValueError('Dark almanac resolution muyst divide into a day '
                             'with no remainder (i.e. resolution must be a '
                             'divisor of 1440).')

        if self.dark_time is None:
            self.dark_time = {}

        logging.debug('Populating dark time dict')
        for i in range(len(dates)):
            self.dark_time[dates[i]] = is_dark_time[i]

        logging.debug('Dark almanac created!')

        return



# ______________________________________________________________________________


def airmass(ra, dec, date_j2000=None, observer=UKST_TELESCOPE):
    # initialise ephem observer object with specified datetime
    if date_j2000 is not None:
        observer.date = date_j2000
    # initialise ephem target object with specified coordinates
    target = ephem.FixedBody()             # create ephem target object
    target._ra = np.radians(ra)          # radians everywhere
    target._dec = np.radians(dec)        # radians everywhere
    target.compute(observer)             # intialise ephem target object

    # return airmass if target is accessible with UKST
    if np.degrees(target.alt) > 5:
        return 1./np.sin(target.alt)
    else:
        return 99.

# ______________________________________________________________________________


def next_dark_period(observer=UKST_TELESCOPE):
    _date = observer.date
    # determine next sun and moon rise and set times from now.
    sunrise, sunset = next_sunrise(observer), next_sunset(observer)
    moonrise, moonset = next_moonrise(observer), next_moonset(observer)
    # are the sun and moon up now? (will they set before they rise?)
    is_sun_up = sunset < sunrise
    is_moon_up = moonset < moonrise
    while (is_sun_up or is_moon_up):
        until = observer.date + 1./1440 # advance at least 1 min into the future
        # advance time to when sun/moon have set (whichever is furthest ahead)
        if is_sun_up:
            until = max(until,  sunset)
        if is_moon_up:
            until = max(until, moonset)
        observer.date = until
        # update rise/set times as necessary
        if sunrise <= until:
            sunrise = next_sunrise(observer)
        if sunset <= until:
            sunset = next_sunset(observer)
        if moonrise <= until:
            moonrise = next_moonrise(observer)
        if moonset <= until:
            moonset = next_moonset(observer)
        # how about now? are the sun and moon up now?
        is_sun_up = sunset < sunrise
        is_moon_up = moonset < moonrise
    darkstart = observer.date
    # dark time ends when either the sun or the moon rises
    darkend = min(sunrise, moonrise)
    logging.debug('next dark period: %s --> %s' % (darkstart, darkend, ))
    # return observer.date to its original value
    observer.date = _date
    return darkstart, darkend        

# ______________________________________________________________________________


def hours_observable(ra, dec, observer=UKST_TELESCOPE,
                     start_date=None, end_date=None,
                     observing_period=365.25 / 2.,
                     minimum_airmass=2.0, verbose=False):

    # Input checking
    if end_date is None and observing_period is None:
        raise ValueError('Either end_date or observing_period must be defined')

    if start_date is None:
        start_date = utc_to_j2000() # Set start time to now

    # determine end of observing period if necessary
    if end_date is None:
        observer.date = start_date  # pyephem accepts datetime object
        end_date = observer.date + observing_period
        observer.date = end_date

    if end_date <= start_date:
        raise ValueError('The end date is before the start date! Refine either '
                         'your passed end_date or observing_period')
        
    # initialise ephem fixed body object for target
    target = ephem.FixedBody()
    target._ra = np.radians(ra)    # radians everywhere
    target._dec = np.radians(dec)  # radians everywhere

    # define angle corresponding to minimum_airmass and twilight_horizon
    observable = np.arcsin(1. / minimum_airmass)  # radians everywhere

    # hack around weird pyephem behavior near target's zenith
    # Precisely what is this doing?
    zenith = np.abs(np.abs(observer.lat) + np.radians(dec - 90.))
    observable = np.clip(observable, 0., zenith - np.radians(.1))
    
    def next_rise(obs):
        obs.horizon = observable
        target.compute(obs)
        return obs.next_rising(target)

    def next_set(obs):
        obs.horizon = observable
        target.compute(obs)
        return obs.next_setting(target)
        
    # from start of observing, compute next rise/set dates 
    observer.date = start_date    # pyephem accepts datetime object
    sunrise = next_sunrise(observer)
    sunset = next_sunset(observer)
    moonrise = next_moonrise(observer)
    moonset = next_moonset(observer)

    # given its Dec, check whether target actually ever sets, or is circumpolar
    target.compute(observer)
    always_up_declination = -np.pi/2. + np.abs(observer.lat) - observable
    target_always_up = dec < np.degrees(always_up_declination)

    # and then compute target rise/set times if necessary
    if target_always_up:
        targetrise, targetset = observer.date, end_date
    else:
        targetrise = next_rise(observer)
        targetset = next_set(observer)

    # and now step forward through the observing period
    hours_total = 0.
    while observer.date < end_date:
        
        # if sun rises before it sets then it is nighttime (good)
        is_night_time = sunrise < sunset 
        # if moon rises before it sets then it is darktime (good)
        is_dark_time = moonrise < moonset
        # if target sets before it rises then it is visible (good)
        is_target_up = (targetset < targetrise)  # or target_always_up

        # if target is observable record how long it is observable for
        if is_night_time and is_dark_time and is_target_up:

            # determine when target becomes unobservable
            until = min(sunrise, moonrise, targetset, end_date)

            # record amount of observable time in hours
            hours_tonight = (until - observer.date) * 24.
            # add this to the total observable time in the full period
            hours_total += hours_tonight

            if verbose:
                print (observer.date, '-->', until,
                       '(%.3f hr)' % hours_tonight)
        else:
            # if target is not observable, determine when it might next be
            # advance by a minimum of one minute to sidestep pyephem oddnesses
            until = observer.date + 1./1440
            if not is_night_time:
                until = max(until, sunset)
            if not is_dark_time:
                until = max(until, moonset)
            if not is_target_up:
                until = max(until, targetrise)

        if observer.date == until:
            pass
            # This code block triggers if until is the same as the osberver date
            # So, that's in the event that either:
            # - The target is observable, but the until unobservable date is
            #   somehow the same;
            # - The target is unobservable, but somehow is magically observable
            #   at the current date
            # print 'oh no!', observer.date
            # print
            # until = observer.date
            # print until, float(until)
            # print is_night_time; until = max(until, sunset)
            # print until, float(until)
            # print is_dark_time;  until = max(until, moonset)
            # print until, float( until )
            # print is_target_up; until = max(until, targetrise)
            # print until, float(until)
            # print
            # print 'currently', observer.date, float(observer.date)
            # print 'sunrise', sunrise, sunrise - observer.date
            # print 'sunset', sunset, sunset - observer.date
            # print 'is_night_time', is_night_time, sunset - observer.date
            # print 'moonrise', moonrise, moonrise - observer.date
            # print 'moonset', moonset, moonset - observer.date
            # print 'is_dark_time', is_dark_time, moonset-observer.date
            # print 'targetrise', targetrise, targetrise - observer.date
            # print 'targetset', targetset, targetset - observer.date
            # print 'is_target_up', is_target_up, targetrise - observer.date
            # print 'checks', targetset <= until, targetset - until, \
            #     targetset - next_set()
            # print 'until', until, until - observer.date
            # print 'target_always_up', target_always_up
#            time.sleep( 10. )

        # advance time to next change in circumstances
        observer.date = min(until, end_date)
        
        # update any/all event times that need updating
        if sunrise <= until:
            sunrise = next_sunrise(observer=observer)
        if sunset <= until:
            sunset = next_sunset(observer=observer)
        if moonrise <= until:
            moonrise = next_moonrise(observer=observer)
        if moonset <= until:
            moonset = next_moonset(observer=observer)
        if targetrise <= until and not target_always_up:
            targetrise = next_rise()
        if targetset <= until and not target_always_up:
            targetset = next_set()

    # all done!
    return hours_total
        
        
# ______________________________________________________________________________

def hours_observable_bruteforce(ra, dec, observer=UKST_TELESCOPE,
                                start_date=None, end_date=None,
                                observing_period=365.25 / 2.,
                                minimum_airmass=2.0, verbose=0,
                                resolution=15., full_output=False,
                                dates_j2000=None):

    # TODO: Add input checking

    # define horizons for observability computation
    observable = np.arcsin(1./ minimum_airmass) # radians everywhere
    solar_horizon = SOLAR_TWILIGHT_HORIZON - np.radians(0.24)
    lunar_horizon = LUNAR_TWILIGHT_HORIZON - np.radians(0.27)
    # remember that the sun and moon both have a size! (0.24 and 0.27 deg)
    
    # initialise ephem fixed body object for pointing
    target = ephem.FixedBody()
    target._ra = np.radians(ra)    # radians everywhere
    target._dec = np.radians(dec)  # radians everywhere

    # instantiate ephem sun/moon objects
    sun, moon = ephem.Sun(), ephem.Moon()

    # define time grid if not given 
    if dates_j2000 is None:
        # initialise ephem observer object for UKST @ SSO
        observer.date = start_date

        if end_date is None:
            end_date = observer.date + observing_period
        else :
            observing_period = end_date - observer.date
        dates_j2000 = (observer.date +
                       np.arange(0, observing_period, resolution / 1440.))

    # and finally the actual computation
    time_total = 0.
    sol_alt, lun_alt, target_alt, dark_time = [], [], [], []
    for d in dates_j2000.ravel():
        observer.date = d
        sun.compute(observer)
        moon.compute(observer)
        target.compute(observer)

        can_observe = (sun.alt < solar_horizon) and (moon.alt < lunar_horizon)
        if full_output:
            sol_alt.append(sun.alt)
            lun_alt.append(moon.alt)
            target_alt.append(target.alt)
            dark_time.append(can_observe)
        else :
            if can_observe and (target.alt > observable):
                time_total += resolution
            
    if full_output :
        # convert results to numpy arrays
        sol_alt = np.array(sol_alt).reshape(dates_j2000.shape)
        lun_alt = np.array(lun_alt).reshape(dates_j2000.shape)
        target_alt=np.array(target_alt).reshape(dates_j2000.shape)
        dark_time=np.array(dark_time).reshape(dates_j2000.shape)
        return dates_j2000, sol_alt, lun_alt, target_alt, dark_time
    else:
        return time_total/60.


# ______________________________________________________________________________

def hours_observable_from_almanac(almanac_filename,
                                  dates_j2000, is_dark_time,
                                  observer=UKST_TELESCOPE,
                                  ra=None, dec=None,
                                  start_date_j2000=None, end_date_j2000=None,
                                  observing_period=365.25 / 2.,
                                  minimum_airmass=2.0):
    # load almanac if it exists
    if os.path.exists(almanac_filename):
        with np.load(almanac_filename) as save_dict:
            ddates, airmass = save_dict['dates_J2000'], save_dict['airmass']
        # check that almanac follows same date grid as input is_dark_time
        if not np.allclose(dates_j2000, ddates):
            logging.info('''\n\nBADNESS!
            dates found in almanac with filename %s do not match those given
            to hours_observable_from_almanac. this almanac must be remade.
            existing almanac is being removed and reconstructed.
            ''' % almanac_filename)
            os.remove(almanac_filename)
            # discard this almanac if its not right

    # create almanac if necessary
    if not os.path.exists(almanac_filename):
        if ra is None or dec is None:
            raise RuntimeError('''BADNESS!
            almanac with filename %s does not exist.
            call hours_observable_from_almanac using ra, dec keywords to
            create a new almanac.  dying gracelessly now.
            ''' % almanac_filename)
        airmass = create_almanac(ra, dec, dates_j2000, almanac_filename)

    if end_date_j2000 is None:
        end_date_j2000 = observer.date + observing_period

    # identify when target is observable during observing period
    fair_game = is_dark_time and \
                (start_date_j2000 <= dates_j2000) and \
                (dates_j2000 <= end_date_j2000)

#    better_entries = np.sum( fair_game & ( airmass <= minimum_airmass ) )
    better_entries = np.interp(minimum_airmass,
                               np.sort(airmass[fair_game]),
                               np.arange(fair_game.sum()))
    # this gives the number of dark time entries within the specified time
    # period where the airmass in the almanac is less than the minimum.

    # use dates grid to determine resolution in hours
    first_date = dates_j2000[np.unravel_index(0, dates_j2000.shape)]
    last_date = dates_j2000[np.unravel_index(dates_j2000.size - 1,
                                              dates_j2000.shape)]
    resolution_in_hours = (last_date - first_date) * 24. / dates_j2000.size
    better_hours = better_entries * resolution_in_hours

    return better_hours
        
    
# ______________________________________________________________________________

def create_dark_almanac(dark_almanac_filename,
                        start_date_j2000=None, end_date_j2000=None,
                        observing_period=int(4 * 365.25),
                        resolution=15.):

    logging.info('creating dark almanac %s ...' % (dark_almanac_filename), )
    # sys.stdout.flush()
    
    dates, sun, moon, target, is_dark_time = hours_observable_bruteforce(
        ra=0., dec=0., 
        start_date=start_date_j2000, end_date=end_date_j2000,
        observing_period=observing_period,
        resolution=resolution, full_output=True)

    solar_horizon = SOLAR_TWILIGHT_HORIZON - np.radians(0.24)
    lunar_horizon = LUNAR_TWILIGHT_HORIZON - np.radians(0.27)
    # remember that the sun and moon both have a size! (0.24 and 0.27 deg)

    times_per_day = 1440. / resolution
    if not times_per_day == int(times_per_day):
        logging.info('''BADNESS!
        resolution_in_min must divide a day into an integer no pieces.
        in other words, resolution_in_min must be a divisor of 1440.
        valid values include 1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60, etc.
        the recommended value is 15.''')
    shape = (dates.size/times_per_day, times_per_day)

    dates_J2000 = dates.reshape(shape)
    sun, moon = sun.reshape(shape), moon.reshape(shape)
    target = target.reshape(shape)
    is_dark_time = is_dark_time.reshape(shape)

    save_dict = { 'dates_J2000' : dates_J2000, 'is_dark_time' : is_dark_time }
    np.savez(dark_almanac_filename, **save_dict)
    logging.info('done.')
    
    return dates_J2000, is_dark_time
    
    
# ______________________________________________________________________________

def create_almanac(ra, dec, dates_j2000, almanac_filename):
    logging.info('creating almanac' % (almanac_filename, ))
    # sys.stdout.flush()
    
    dates, sun, moon, target, dark_time = hours_observable_bruteforce(
        ra, dec, dates_j2000=dates_j2000, full_output=True)

    airmass = np.clip(np.where(target > np.radians(10.),
                               1./np.sin(target), 99.), 0., 9.)

    save_dict = {'dates_J2000': dates_j2000, 'airmass': airmass}
    np.savez(almanac_filename, **save_dict)
    logging.info('done.')

    return airmass

        
# ______________________________________________________________________________
#
# CONVENIENCE FUNCTIONS FOR RISE AND SET TIMES

# ephem does everything internally using J2000 dates

def utc_to_j2000(datetime_utc=datetime.datetime.utcnow()):
    return ephem.julian_date(datetime_utc) - 2415020.


def next_sunrise(observer=UKST_TELESCOPE, horizon=SOLAR_TWILIGHT_HORIZON):
    observer.horizon = horizon             # specify horizon for rise/set
    target = ephem.Sun()                   # set target as Sol
    target.compute(observer)             # intialise ephem target object
    return observer.next_rising(target)  # return next rise date


def next_sunset(observer=UKST_TELESCOPE, horizon=SOLAR_TWILIGHT_HORIZON):
    observer.horizon = horizon             # specify horizon for rise/set
    target = ephem.Sun()                   # set target as Sol
    target.compute(observer)             # intialise ephem target object
    return observer.next_setting(target) # return next set date


def next_moonrise(observer=UKST_TELESCOPE, horizon=LUNAR_TWILIGHT_HORIZON):
    observer.horizon = horizon             # specify horizon for rise/set
    target = ephem.Moon()                  # set target as the moon
    target.compute(observer)             # intialise ephem target object
    return observer.next_rising(target)  # return next rise date


def next_moonset(observer=UKST_TELESCOPE, horizon=LUNAR_TWILIGHT_HORIZON):
    observer.horizon = horizon             # specify horizon for rise/set
    target = ephem.Moon()                  # set target as the moon
    target.compute(observer)             # intialise ephem target object
    return observer.next_setting(target)  # return next set date


def next_riseset(RA, Dec, observer=UKST_TELESCOPE, horizon_in_degrees=0.):
    observer.horizon = np.radians(horizon_in_degrees) # radians everywhere
    target = ephem.FixedBody()             # create ephem target object
    target._ra = np.radians(RA)          # radians everywhere
    target._dec = np.radians(Dec)        # radians everywhere
    target.compute(observer)             # intialise ephem target object
    risedate = observer.next_rising(target)
    setdate = observer.next_setting(target)
    return risedate, setdate               # return rise and set dates

# ______________________________________________________________________________
