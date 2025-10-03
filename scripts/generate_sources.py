from utils import *
import json
import random
import datetime
import argparse
import warnings
import astropy.units as u
from utils import Source
warnings.filterwarnings("ignore", category=UserWarning, append=True)


parser = argparse.ArgumentParser(description="Generate random sky sources.")
parser.add_argument('N', type=int, help='Number of sources to generate')
parser.add_argument('--center', type=str, default="10:00:28.3736076266 2:13:12.0064477050", help='Central RA in hms and Dec in dms')
parser.add_argument('--radius', type=str, default="1.25deg", help='FOV for source distribution. Use astropy units (e.g. 1.5deg, 2rad, 30arcmin). If not units, degrees assumed.')
parser.add_argument('--max_sep', type=float, default=0.0, help='Maximum separation between sources in arcsec (0 means no limit)')
parser.add_argument('--min_sep', type=float, default=0.0, help='Minimum separation between sources in arcsec (0 means no limit)')
parser.add_argument('--intensity', type=float, default=1, help='Base intensity in Jy. Range between 0.5I and 1.5I')
parser.add_argument('--prefix', type=str, default=None, help='Prefix for output filename')
parser.add_argument('--ref_freq', type=float, default=1310, help='Reference frequency for spectral index in MHz')
parser.add_argument('--major_axis', type=float, default=5.0, help='Major axis in arcsec')
parser.add_argument('--fits', action='store_true', help='Output in FITS format (default is JSON)')
parser.add_argument('--author', type=str, default='Diaz-Gonzalez, D.J.', help='Author name to include in FITS header')
parser.add_argument('--omegabeam', type=float, default=None, help='Beam area in arcsec^2 for converting flux to Jy/beam. Default is SKA-Mid-AA* beam')

args = parser.parse_args()

if (args.omegabeam is not None) and (args.omegabeam <= 0):
    raise ValueError("Beam area must be a positive value in arcsec^2")

if args.omegabeam is None:
    omegabeam = 3.38551 * u.arcsec  * 2.95 * u.arcsec * np.pi / (4 * np.log(2))
    print(f"No beam area provided, assuming SKA-Mid-AA* beam at 1.4 GHz ({omegabeam.to(u.arcsec**2):.2f})")
else:
    omegabeam = args.omegabeam * u.arcsec**2


N_sources = args.N

limit = u.Quantity(args.radius)
if not limit.unit.is_equivalent(u.deg):
    raise ValueError("Radius must be an angular quantity with astropy units (e.g. 1.5deg, 2rad, 30arcmin)")
# If no units, assume degrees
if limit.unit == u.dimensionless_unscaled:
    limit = limit * u.deg
    
center = SkyCoord(args.center, unit=(u.hourangle, u.deg), frame='icrs')

prefix = args.prefix
if prefix is None:
    prefix = datetime.datetime.now().strftime("%Y%m%d_%H%M") + '_'
if not prefix.endswith('_') and len(prefix) > 0:
    prefix = f"{prefix}_"

sources = []
while len(sources) < N_sources:

    intensity = args.intensity  * random.uniform(0.5, 1.5) * u.Jy 
    spec_index = 0 # random.uniform(0, 4) - 2
    major_axis = args.major_axis * random.uniform(0.5, 1.5) * u.arcsec
    minor_axis = major_axis * random.uniform(0.5, 1.0)
    coord  = center.directional_offset_by(position_angle=
                                          random.uniform(0,360)*u.deg, separation= ((limit - 2 * major_axis) * np.sqrt(random.random())))
    print("Generating source at ", coord.to_string('hmsdms'), " with offset ", coord.separation(center).to(u.arcsec))

    area_source = (np.pi * major_axis * minor_axis) / (4 * np.log(2))
    intensity = intensity * (area_source/omegabeam)  # Convert to Jy/beam assuming Gaussian source and beam
    pa = random.uniform(0, 180) * u.deg
    source = Source(ra=coord.ra.deg, dec=coord.dec.deg, I=intensity, spec_index=spec_index, ref_freq=args.ref_freq * u.MHz, major_axis=major_axis, minor_axis=minor_axis, pa=pa, rot_meas = 0.0 * u.rad/u.m**2)
    if args.min_sep > 0 or args.max_sep > 0:
        dists = np.array([source.coord.separation(s.coord).arcsec for s in sources])
        if args.min_sep > 0:
            if np.any(dists < args.min_sep):
                print(f" - Rejected (too close to another source, min_sep={args.min_sep} arcsec)")
                continue
        if args.max_sep > 0:
            if np.any(dists > args.max_sep):
                print(f" - Rejected (too far from another source, max_sep={args.max_sep} arcsec)")
                continue
    sources.append(source)

if not args.fits:
    json_data = [source.to_json() for source in sources]
    with open(f'{prefix}{N_sources:03d}_sources.json', 'w') as f:
        json.dump(json_data, f, indent=4)
    print (f"Saved {N_sources} sources to {prefix}{N_sources:03d}_sources.json")
else:
    from astropy.io import fits
    from astropy.table import Table
    ra_list = np.array([source.ra.to(u.deg).value for source in sources])
    dec_list = np.array([source.dec.to(u.deg).value for source in sources])
    I_list = np.array([source.I.to(u.Jy).value for source in sources])
    Q_list = np.array([source.Q.to(u.Jy).value for source in sources])
    U_list = np.array([source.U.to(u.Jy).value for source in sources])
    V_list = np.array([source.V.to(u.Jy).value for source in sources])


    sp_list = np.array([source.spec_index for source in sources])
    ref_freq_list = np.array([source.ref_freq.to(u.MHz).value for source in sources])
    rm_list = np.array([source.rot_meas.to(u.rad / u.m**2).value for source in sources])
    major_axis_list = np.array([source.major_axis.to(u.arcsec).value for source in sources])
    minor_axis_list = np.array([source.minor_axis.to(u.arcsec).value for source in sources])
    pa_list = np.array([source.pa.to(u.deg).value for source in sources])
    id_list = np.array([f"SRC{i:03d}" for i in range(N_sources)])



    cols = [
        fits.Column(name='ID', format='10A', array=id_list),
        fits.Column(name='RA', format='D', unit='deg', array=ra_list),
        fits.Column(name='DEC', format='D', unit='deg', array=dec_list),
        fits.Column(name='STK_I', format='E', unit='Jy/beam', array=I_list),
        fits.Column(name='STK_Q', format='E', unit='Jy/beam', array=Q_list),
        fits.Column(name='STK_U', format='E', unit='Jy/beam', array=U_list),
        fits.Column(name='STK_V', format='E', unit='Jy/beam', array=V_list),
        fits.Column(name='SPECIDX', format='E', array=sp_list),
        fits.Column(name='REFFREQ', format='E', unit='MHz', array=ref_freq_list),
        fits.Column(name='RM', format='E', unit='rad/m^2', array=rm_list),
        fits.Column(name='MAJ', format='E', unit='arcsec', array=major_axis_list),
        fits.Column(name='MIN', format='E', unit='arcsec', array=minor_axis_list),
        fits.Column(name='PA', format='E', unit='deg', array=pa_list),
        fits.Column(name='RESOLVED', format='L', array=[0]*N_sources),  # Placeholder column
        fits.Column(name='ISL_RMS', format='E', unit='Jy/beam', array=[random.random()*1e-2]*N_sources),  # Placeholder column
    ]

    coldefs = fits.ColDefs(cols)
    table = fits.BinTableHDU.from_columns(coldefs)


    prihdu = fits.PrimaryHDU()
    prihdu.header['AUTHOR'] = args.author
    prihdu.header['CREATOR'] = 'SKA-Synth'
    prihdu.header['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d")

    hdu1 = fits.HDUList([prihdu, table])
    hdu1.writeto(f'{prefix}{N_sources:03d}_sources.fits', overwrite=True)
    print (f"Saved {N_sources} sources to {prefix}{N_sources:03d}_sources.fits")



