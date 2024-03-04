# Example manipulating GX_Simulator output: counts spectrum from photon spectrum
## Please download the data files from [here](https://drive.google.com/drive/folders/1TlkJ4RoDRLIhGxWXoLCCKpD8aMarn1Cq) and read the [accompanying readme.txt](https://drive.google.com/file/d/1uYQ5aEZFd3Z6R-rv-TB-JXWVb4l1GkLv) for more information. The readme contents are replicated in this folder's README.md for convenience.

### `readme.txt` from Google Drive
The example converts the GX_Simulator produced X-ray photon spectrum
into a counts spectrum by folding the spectrum through an instrument response.

Care must be taken to properly integrate the photon spectrum into the
instrument photon model bins; the GX_Simulator output is a function of energy,
and is not binned initially.

The .sav file contains the X-ray photon spectrum in ph / keV / cm2 / s.
The data files contain counts produced by the STIX SSWIDL ground software
and the raw Fermi/GBM counts/response files.

Because Fermi/GBM rotates to achieve full-sky coverage, its instrument response
changes as a function of time. Therefore care must be taken to choose
the correct response matrix (RM) from its .rsp2 file.
The valid time ranges are populated in the FITS header of each entry.
In this example the response matrix index is hard-coded but examining the
DATE-OBS and DATE-END FITS header keywords could allow automated picking of the RM.
Problems arise when fitting across RM regions. Don't do it.

