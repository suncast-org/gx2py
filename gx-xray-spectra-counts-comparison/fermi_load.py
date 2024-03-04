import astropy.io.fits as fits
import astropy.time as atime
import astropy.units as u
import numpy as np

from dataclasses import dataclass
from datetime import datetime as dtime

NOMINAL_DT = 256 << u.ms
NOMINAL_BURST_DT = 64 << u.ms
MISSION_BASE_TIME = atime.Time(dtime.fromisoformat('2001-01-01T00:00'), format='datetime')
DETECTOR_AREA = 100 << u.cm**2

nearest = lambda a, v: np.argmin(np.abs(a - v))


@dataclass
class FermiDat:
    counts: u.Quantity
    dt: u.Quantity
    time_bins: dtime
    energy_bins: u.Quantity


@dataclass
class FermiTimeTaggedEventList:
    events: u.Quantity
    energy_bins: u.Quantity
    channels: np.ndarray
    times: np.ndarray

    @property
    def pretty_times(self) -> atime.Time:
        return MISSION_BASE_TIME + atime.TimeDelta(self.times, format='sec')

    @classmethod
    def from_fits(cls, fn: str):
        with fits.open(fn) as f:
            energy_bins = np.column_stack(
                (f['ebounds'].data['e_min'], f['ebounds'].data['e_max'])
            ) << u.keV
            channels = f['ebounds'].data['channel']

            # leave as numbers for speed in slicing and other manipulations
            evt_times = f['events'].data['time']
            evts = f['events'].data['pha']
        return cls(evts, energy_bins, channels, evt_times)

    def slice(self, start_time: str | atime.Time, end_time: str | atime.Time):
        s, e = atime.Time(start_time), atime.Time(end_time)

        # convert start, end to just # seconds to speed up index
        # searching (use C searches, avoid astropy internals)
        # important when TTE has like 20,000,000 events ;)
        s = (s - MISSION_BASE_TIME).to_value(u.s)
        e = (e - MISSION_BASE_TIME).to_value(u.s)
        si, ei = nearest(s, self.times), nearest(e, self.times)

        self.times = self.times[si:ei+1]
        self.events = self.events[si:ei+1]

    @u.quantity_input
    def bin_to_spectrogram(self, ts: atime.Time, te: atime.Time, dt: u.s) -> dict[str, np.ndarray]:
        start = (ts - MISSION_BASE_TIME).to_value(u.s)
        stop = (te - MISSION_BASE_TIME).to_value(u.s)
        bins = np.arange(start, stop, dt.to_value(u.s))

        # binning indices to apply to energy
        bin_indices = np.digitize(self.times, bins)

        # E, t spectrogram to return
        spectrogram = np.zeros((bins.size - 1, self.channels.max()))

        for time_idx in range(bins.size - 1):
            these = (bin_indices == time_idx)
            binned_timestep, _ = np.histogram(self.events[these], self.channels)
            spectrogram[time_idx] = binned_timestep

        return {
            'time_bins': atime.Time(MISSION_BASE_TIME + atime.TimeDelta(bins, format='sec')),
            'spectrogram': spectrogram.T.astype(int)
        }


def load_ctime(fits_file: str) -> FermiDat:
    with fits.open(fits_file) as gbm:
        energy_bins = np.column_stack((gbm['ebounds'].data['e_min'], gbm['ebounds'].data['e_max'])) << u.keV
        cts = gbm['spectrum'].data['counts'] << u.ct

        dt = gbm['spectrum'].data['exposure'] << u.s
        time_bins = (
            atime.Time(MISSION_BASE_TIME, format='datetime') +
            atime.TimeDelta(gbm['spectrum'].data['time'], format='sec')
        )
        time_bins = np.concatenate((time_bins, [time_bins[-1] + dt[-1]]))
        time_bins = atime.Time(time_bins)
        time_bins.precision = 9
        time_bins = time_bins.to_datetime()

        return FermiDat(counts=cts, dt=dt, time_bins=time_bins, energy_bins=energy_bins)


def slice_ctime(dat: FermiDat, start: str | dtime, end: str | dtime) -> FermiDat:
    ''' slice ctime data by UTC time '''
    tslice_start, tslice_end = start, end
    if isinstance(start, str):
        tslice_start = dtime.fromisoformat(start)
    if isinstance(end, str):
        tslice_end = dtime.fromisoformat(end)

    si, ei = nearest(dat.time_bins, tslice_start), nearest(dat.time_bins, tslice_end)

    sliced_bins = dat.time_bins[si:ei]
    sliced_dt = atime.TimeDelta(dat.dt[si:ei-1]).to(u.s)
    sliced_cts = dat.counts[si:ei-1]

    return FermiDat(sliced_cts, sliced_dt, sliced_bins, dat.energy_bins)


def bin_down_ctime(dat: FermiDat) -> FermiDat:
    ''' bin ctime dat down in time to the big bins
        returns: rebinned and livetime-adjusted data
    '''
    if not np.all(np.isfinite(dat.counts)):
        raise ValueError("can only rebin GBM times when there are no data gaps")

    time_rebinned = []
    delta = 0.75
    i = 0
    nom = (1 - (dat.dt / NOMINAL_DT)) < delta
    while i < dat.counts.shape[0]:
        if nom[i]:
            time_rebinned.append(dat.counts[i])
            i += 1
        else:
            time_rebinned.append(dat.counts[i:i+4].sum(axis=0))
            i += 4

    aligned_cts = np.array(time_rebinned) << u.ct

    new_time_bins = (
        atime.Time(dat.time_bins[0], format='datetime') +
        NOMINAL_DT * np.arange(aligned_cts.shape[0] + 1)
    )

    dt = NOMINAL_DT * np.ones(aligned_cts.shape[0])
    return FermiDat(
        counts=aligned_cts,
        dt=dt,
        time_bins=new_time_bins.to_datetime(),
        energy_bins=dat.energy_bins
    )

