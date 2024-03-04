import astropy.table as atable
import astropy.units as u
from astropy.io import fits
import numpy as np

def stix_srm_for_sunxspex(fn: str, which: str='unattenuated') -> dict[str, np.ndarray]:
    try:
        att_index_map = {'unattenuated': 1, 'attenuated': 4}
        srm_idx = att_index_map[which]
    except KeyError:
        raise ValueError(
            f'{which} is not a valid attenuator choice (choose from attenuated, unattenuated)'
        )

    with fits.open(fn) as f:
        srm_dat = atable.QTable.read(f[srm_idx], format='fits')
        ct_energy_dat = atable.QTable.read(f['ebounds'], format='fits')
        area = f[srm_idx].header['geoarea'] << u.cm**2

        srm = np.array(srm_dat['MATRIX']) << (u.ct / u.keV / u.ph)
        model_edges = np.column_stack((srm_dat['ENERG_LO'], srm_dat['ENERG_HI']))
        count_edges = np.column_stack((ct_energy_dat['E_MIN'], ct_energy_dat['E_MAX']))
        srm *= area * np.diff(count_edges, axis=1).flatten()
        srm = srm.to(u.cm**2 * u.ct / u.ph)

    model_edges = model_edges.to_value(u.keV)
    count_edges = count_edges.to_value(u.keV)
    return {
        'srm': srm.to_value(u.cm**2 * u.ct / u.ph),
        'photon_channel_bins': model_edges,
        'photon_channel_mids': model_edges.mean(axis=1).flatten(),
        'photon_channel_binning': np.diff(model_edges, axis=1).flatten(),
        'count_channel_bins': count_edges,
        'count_channel_mids': count_edges.mean(axis=1).flatten(),
        'count_channel_binning': np.diff(count_edges, axis=1).flatten(),
        'area': area
    }

