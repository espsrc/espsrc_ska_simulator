#from karabo.simulation.sky_model import SkyModel
from utils import SkyModel
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=UserWarning, append=True)

# PHASE_CENTER_RA = 150.12
# PHASE_CENTER_DEC = 2.21
# sky_model = SkyModel.get_MIGHTEE_Sky()
sky_model = SkyModel.get_GLEAM_Sky()

print (f"Number of sources in the sky model: {len(sky_model.sources)}")
# center = sky_model.get_center()
# PHASE_CENTER_RA = center.ra.deg
# PHASE_CENTER_DEC = center.dec.deg
# print(f"Phase Center: {center}")
# # print(f"Phase Center: RA={PHASE_CENTER_RA}, Dec={PHASE_CENTER_DEC}")

# # print (sky_model.wcs.crval)
# sky_model.explore_sky([PHASE_CENTER_RA, PHASE_CENTER_DEC], block=True, cfun=np.log10,)


