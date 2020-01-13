# __init__.py: Top level .py file for python stationresponse interface
# Copyright (C) 2011
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id: __init__.py 33126 2015-12-11 19:09:42Z dijkema $
from __future__ import absolute_import, print_function, division

from . import _stationresponse

class stationresponse(object):
    """
    The Python interface to the LOFAR station beam model.
    """

    def __init__ (self, msname, inverse = False, useElementResponse = True,
        useArrayFactor = True, useChanFreq = False):
        """Create a stationresponse object that can be used to evaluate the
        LOFAR beam for the given Measurement Set.

        The Measurement Set defines the station and dipole positions, the phase
        center, and the channel frequencies (and reference frequency) for which
        the LOFAR beam will be evaluated.

        `msname`
          Name of the Measurement Set.
        `inverse`
          Compute the inverse of the LOFAR beam (default False).
        `useElementResponse`
          Include the effect of the dual dipole (element) beam (default True).
        `useArrayFactor`
          Include the effect of the station and tile array factor (default
          True).
        `useChanFreq`
          Compute the phase shift for the station beamformer using the channel
          frequency instead of the subband reference frequency. This option
          should be enabled for Measurement Sets that contain multiple subbands
          compressed to single channels inside a single spectral window
          (default: False).

        For example::

        import pyrap.tables
        import lofar.stationresponse
        response = lofar.stationresponse.stationresponse('test.MS')

        # Iterate over all time stamps in the Measurement Set and compute the
        # beam Jones matrix for station 0, channel 0.
        ms = pyrap.tables.table('test.MS')
        for subtable in ms.iter('TIME'):
            time = subtable.getcell("TIME", 0)
            print time, response.evaluateChannel(time, 0, 0)
        """
        self._response = _stationresponse.StationResponse(msname, inverse,
          useElementResponse, useArrayFactor, useChanFreq)

    def version (self, type='other'):
        """Show the software version."""
        return self._response.version (type)

    def setRefDelay (self, ra, dec):
        """Set the reference direction used by the station beamformer. By
        default, DELAY_DIR of field 0 is used.

        `ra`
          Right ascension (in radians, J2000)
        `dec`
          Declination (in radians, J2000)
        """
        self._response.setRefDelay(ra, dec)

    def getRefDelay (self, time):
        """Get the reference direction used by the station beamformer.
        Returns an ITRF vector in meters (numpy array of 3 floats).
 
        `time`
          Time at which to evaluate the direction
        """
        return self._response.getRefDelay(time)

    def setRefTile (self, ra, dec):
        """Set the reference direction used by the analog tile beamformer
        (relevant for HBA observations only). By default, LOFAR_TILE_BEAM_DIR
        of field 0 is used. If not present, DELAY_DIR of field 0 is used
        instead.

        `ra`
          Right ascension (in radians, J2000)
        `dec`
          Declination (in radians, J2000)
        """
        self._response.setRefTile(ra, dec)

    def getRefTile (self, time):
        """Get the reference direction used by the analog tile beamformer
        (relevant for HBA observations only).
        Returns an ITRF vector in meters (numpy array of 3 floats).
 
        `time`
          Time at which to evaluate the direction
        """
        return self._response.getRefTile(time)

    def setDirection (self, ra, dec):
        """Set the direction of interest (can be and often will be different
        from the pointing). By default, PHASE_DIR of field 0 is used.

        `ra`
          Right ascension (in radians, J2000)
        `dec`
          Declination (in radians, J2000)
        """
        self._response.setDirection(ra, dec)

    def getDirection (self, time):
        """Get the direction of interest.
        Returns an ITRF vector in meters (numpy array of 3 floats).

        `time`
          Time at which to evaluate the direction
        """
        return self._response.getDirection(time)

    def evaluate (self, time):
        """Compute the beam Jones matrix for all stations and channels at the
        given time. The result is returned as a 4-dim complex numpy array with
        shape: no. of stations x no. of channels x 2 x 2.

        `time`
          Time (MJD in seconds)
        """
        return self._response.evaluate0(time)

    def evaluateStation (self, time, station):
        """Compute the beam Jones matrix for all channels at the given time for
        the given station. The result is returned as a 3-dim complex numpy array
        with shape: no. of channels x 2 x 2.

        `time`
          Time (MJD in seconds).
        `station`
          Station number (as in the ANTENNA table of the Measurement Set).
        """
        return self._response.evaluate1(time, station)

    def evaluateChannel (self, time, station, channel):
        """Compute the beam Jones matrix for the given time, station, and
        channel. The result is returned as a 2-dim complex numpy array with
        shape: 2 x 2.

        `time`
          Time (MJD in seconds).
        `station`
          Station number (as defined in the ANTENNA table of the Measurement
          Set).
        `channel`
          Channel number (as in the SPECTRAL_WINDOW table of the Measurement
          Set).
        """
        return self._response.evaluate2(time, station, channel)
    
    def evaluateFreq (self, time, station, freq):
        """Compute the beam Jones matrix for the given time, station, and
        frequency. The result is returned as a 2-dim complex numpy array with
        shape: 2 x 2.

        `time`
          Time (MJD in seconds).
        `station`
          Station number (as defined in the ANTENNA table of the Measurement
          Set).
        `frequency`
          Frequency to compute beam at (in Hz)
        """
        return self._response.evaluate3(time, station, freq)

    def evaluateFreqITRF (self, time, station, freq, direction, station0, tile0):
        """Compute the beam Jones matrix for the given time, station, and
        frequency, with the given ITRF directions.
        The result is returned as a 2-dim complex numpy array with
        shape: 2 x 2.

        `time`
          Time (MJD in seconds).
        `station`
          Station number (as defined in the ANTENNA table of the Measurement
          Set).
        `frequency`
          Frequency to compute beam at (in Hz)
        `direction`
          ITRF direction to compute beam at (numpy array with 3 floats)
        `station0`
          ITRF direction of the station beamformer (numpy array with 3 floats)
        `tile0`
          ITRF direction of the tile beamformer (numpy array with 3 floats)
        """
        return self._response.evaluate4(time, station, freq, direction, station0, tile0)
