from karabo.simulation.sky_model import SkyModel

skyModel = SkyModel.get_MIGHTEE_Sky()
phase_center_ra = 150.12
phase_center_dec = 2.21
skyModel.explore_sky([phase_center_ra, phase_center_dec], block=True)