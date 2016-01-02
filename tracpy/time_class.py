"""
Object for keeping track of time in TracPy
"""

import time


class Time(object):
    """
    All times except for starttime and time_at_last_call are time differences
    """

    def __init__(self):
        '''
        Initialize instance of tracpy timer.
        '''

        self.starttime = time.time()  # start time of simulation

        self.times = {}  # initialized dictionary of time flags

        self.total = 0.  # will keep track of total time for simulation

        # This will get reset each time Time is called for calculating time
        # differences
        self.time_at_last_call = time.time()

    def addtime(self, name):
        """
        Add a dt of time spent on a section of code to a string name. If used
        before, it will be summed, otherwise it will start from 0.
        """

        # Add dt to time of section of code `name`. Default is 0, in case
        # it hasn't been used yet.
        # amount of time elapsed since last Time call
        dt = (time.time() - self.time_at_last_call)
        # update the time flag called
        self.times[name] = self.times.get(name, 0) + dt

        # Update total time too
        self.total += dt

        # Update to current time
        self.time_at_last_call = time.time()

    def write(self):
        """Write out all available times."""

        print "Total run time: %f (seconds)" % self.total
        print "---------------------------------------------"
        print "Time spent on:"

        for key in sorted(self.times.keys()):
            print "\t%s \t\t%4.4f (%4.4f%%)" % \
                     (key, self.times[key], (self.times[key]/self.total)*100)

        print "============================================="
