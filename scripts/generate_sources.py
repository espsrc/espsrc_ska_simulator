from utils import *
import json
import random
import datetime
import argparse
import warnings
warnings.filterwarnings("ignore", category=UserWarning, append=True)


parser = argparse.ArgumentParser(description="Generate random sky sources.")
parser.add_argument('N', type=int, help='Number of sources to generate')
parser.add_argument('--ra_0', type=float, default=32.42007, help='Central RA in degrees')
parser.add_argument('--dec_0', type=float, default=-10.1839, help='Central Dec in degrees')
parser.add_argument('--radius', type=float, default=1.25, help='Limit in degrees for source distribution')
parser.add_argument('--intensity', type=float, default=30, help='Base intensity in Jy. Range between 0.5I and 1.5I')
parser.add_argument('--prefix', type=str, default=None, help='Prefix for output filename')

args = parser.parse_args()

N_sources = args.N
limit = args.radius * u.deg

ra_0 = args.ra_0 * u.deg
dec_0 = args.dec_0 * u.deg

prefix = args.prefix
if prefix is None:
    prefix = datetime.datetime.now().strftime("%Y%m%d_%H%M") + '_'
if not prefix.endswith('_') and len(prefix) > 0:
    prefix = f"{prefix}_"


sources = [Source(ra_0  + random.uniform(-1, 1) * limit, dec_0 + random.uniform(-1,1) * limit, args.intensity * random.uniform(0.5,1.5), spec_index= random.uniform(0, 4) -2) for _ in range(N_sources)]
minIntensity = np.min([source.I.to(u.Jy).value for source in sources])
maxIntensity = np.max([source.I.to(u.Jy).value for source in sources])
json_data = [source.to_json() for source in sources]
with open(f'{prefix}{N_sources:03d}_sources.json', 'w') as f:
    json.dump(json_data, f, indent=4)

# skyModel = SkyModel.from_json(json_data)
# skyModel.show(block=True, vmin = minIntensity, vmax=maxIntensity, cmap='viridis', title='Test SkyModel', figsize=(10, 10), cfun=np.abs)

