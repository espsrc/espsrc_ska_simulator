# from karabo.imaging.image import Image
import sys
from ska_img import SKAImage
from utils import Source, show_exc
from radio_beam import Beam
from astrodendro import Dendrogram
from astrodendro.analysis import PPStatistic
import astropy.units as u
from astropy.io import fits
import numpy as np
import argparse
import matplotlib.pyplot as plt
import math, time
import datetime

def create_suffix(s: str) -> str:
    # Create a random 5 alfanum character suffix
    if s == "":
        import random, string
        s = ''.join(random.choices(string.ascii_lowercase + string.digits, k=5))
    if not s.startswith("_"):
        s = "_" + s
    return s



if __name__ == "__main__":
    t0 = time.time()
    try:
        parser = argparse.ArgumentParser(
            description="Test source detection on a simulated sky model"
        )


        parser.add_argument( "files", type=str, nargs="*", help="List of input files (not used in this test)", )
        parser.add_argument( "--components", type=int, default=500, help="Number of components for source modeling", )
        parser.add_argument( "--snr-threshold", type=float, default=1.0, help="SNR threshold for source detection", )
        parser.add_argument( "--rms", type=float, default=0.0, help="RMS noise level in Jy/beam (if 0, it will be estimated from the image)", )
        parser.add_argument( "--trunks", action="store_true", help="Also plot trunks", )
        parser.add_argument("--suffix", type=str, default="", help="Suffix for output files", )
        parser.add_argument("--fsources", type=str, default="", help="File with source list to overlay", )
        parser.add_argument('--author', type=str, default='Diaz-Gonzalez, D.J.', help='Author name to include in FITS header')
        parser.add_argument('--beam', type=str, default='', help='Beam parameters to use if not in FITS header (e.g., "10arcsec,8arcsec,30deg")')
        args = parser.parse_args()


        fig = plt.figure(figsize=(6 * len(args.files), 6))
        suffix = create_suffix(args.suffix)


        # image: SKAImage
        png_filename = f"plot_astro_test{suffix}.png"
        sources = []
        default_beam = None
        if args.beam != "":
            try:
                beam_params = args.beam.split(",")
                major = u.Quantity(beam_params[0])
                minor = u.Quantity(beam_params[1])
                pa = u.Quantity(beam_params[2])
                default_beam = Beam(major=major, minor=minor, pa=pa)
                print (f"Using user-defined beam: {default_beam.major}, {default_beam.minor}, {default_beam.pa}")
            except Exception as e:
                print (f"Error parsing beam parameters: {e}")

        for idx,path in enumerate(args.files):
            image = SKAImage(path=path)
            # Filter components
            mad = image.mad()
            flat_data = image.data.flatten()
            sorted_indices = flat_data.argsort()[::-1]  # Indices of sorted data in descending order
            threshold_index = sorted_indices[args.components - 1]
            threshold_value = flat_data[threshold_index]
            image.data[image.data < threshold_value] = 0.0


            image.name = ('_').join(path.split("/")[-2].split("_")[2:]).replace("_", " ")
            if default_beam is None:
                default_beam = Beam(image.pixelarea**0.5 * 5, image.pixelarea**0.5 * 5, 0*u.deg)

            image = image.convolve_to_beam(default_beam)
            image.data_to_beam(default_beam)
            image.tofits(f"test_convolved{suffix}.fits")
            # Reload image after convolution
            image = SKAImage(path=f"test_convolved{suffix}.fits")
            beam = image.get_beam(defaultBeam=default_beam)
            wcs = image.wcs2d
            pixel_area = wcs.proj_plane_pixel_area()
            pixels_in_beam = math.ceil((beam.sr).to(u.arcsec**2) / (pixel_area.to(u.arcsec**2)))
            print (f"\tCenter: {image.get_center().to_string('hmsdms')}")
            print (f"\tBeam: {beam.major.to(u.arcsec):.4f} x {beam.minor.to(u.arcsec):.4f} ({beam.pa:.1f})")
            print (f"\tShape: {image.data2d.shape}")
            print (f"\tPixels in beam: {pixels_in_beam:.1f}")
            image.data_to_beam()
            if args.rms > 0.0:
                rms = args.rms
                print (f"\tUsing user-defined RMS: {rms:.10f} Jy/beam")
            else:
                rms = image.mad(verbose=True)
                print (f"\tEstimated RMS: {rms:.10f} Jy/beam")
            d = Dendrogram.compute(image.data2d, min_value=threshold_value, min_delta=rms, min_npix=pixels_in_beam)
            # p = d.plotter()
            leaves = d.leaves
            ax = fig.add_subplot(1, len(args.files), idx +1, projection=image.wcs.celestial)
            # ax.imshow(image.data2d, origin='lower', cmap='viridis', interpolation='nearest',vmin=0, vmax=np.percentile(image.data2d, 99))
            image.draw(img=ax, title=f"{image.name} ({len(leaves)} sources)", exp=0.5, cmapstr='viridis', colorbar=True, show_beam=True, nan2zero=True)
            metadata = {}
            print (f"\tData unit converted to Jy/beam? :> {image.bunit}")
            metadata["data_unit"] =image.bunit
            metadata["spatial_scale"] = math.sqrt(pixel_area.to(u.arcsec**2).value) * u.arcsec
            metadata["beam_major"] = beam.major
            metadata["beam_minor"] = beam.minor
        

            for leaf in leaves:
                mask = leaf.get_mask()
                ax.contour(mask, colors='red', linewidths=0.5, alpha=0.7)
                stat = PPStatistic(leaf, metadata=metadata)
                coord = image.pix2coords(stat.x_cen, stat.y_cen, skyMode=True)
                ra = coord.ra.deg
                dec = coord.dec.deg
                data = image.data2d[mask]
                rms_src =np.std(data)
                n_pix = leaf.get_npix()
                
                src = Source(ra=ra*u.deg, dec=dec*u.deg, obj_id=f"SRC{leaf.idx:03d}", I=stat.flux, spec_index=0, ref_freq=image.restfreq(), major_axis=stat.major_sigma, minor_axis=stat.minor_sigma, pa=stat.position_angle, isl_rms=rms_src, resolved=n_pix > pixels_in_beam)
                sources.append(src)
                
            if args.trunks:
                for trunk in d.trunk:
                    mask = trunk.get_mask()
                    ax.contour(mask, colors='blue', linewidths=0.5, alpha=0.7)
                    

            if args.fsources != "":
                from astropy.table import Table
                tsources = Table.read(args.fsources)
                sources = Source.from_table_in_fits(tsources)
                for source in sources:
                    image.draw_source(ax, source, color="cyan", ellipse=True, alpha=0.5, label=None, verbose=False, lw=1)

            if len(sources) > 0:
                N_sources = len(sources)

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
                id_list = np.array([f"SRC{i:03d}" for i in range(len(sources))])

                resolved_list = np.array([source.resolved for source in sources])
                isl_rms_list = np.array([source.isl_rms.to(u.Jy).value for source in sources])



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
                    fits.Column(name='RESOLVED', format='L', array=resolved_list),  # Placeholder column
                    fits.Column(name='ISL_RMS', format='E', unit='Jy/beam', array=isl_rms_list),  # Placeholder column
                ]
                coldefs = fits.ColDefs(cols)
                table = fits.BinTableHDU.from_columns(coldefs)


                prihdu = fits.PrimaryHDU()
                prihdu.header['AUTHOR'] = args.author
                prihdu.header['CREATOR'] = 'SKA-Synth'
                prihdu.header['DATE'] = datetime.datetime.now().strftime("%Y-%m-%d")

                path_fits = f"test_astro_sources{suffix}.fits"
                hdu1 = fits.HDUList([prihdu, table])
                hdu1.writeto(path_fits, overwrite=True)
                print (f"Saved {N_sources} sources to {path_fits}")



            


            print (f"\tFound {len(leaves)} sources")
        plt.tight_layout()
        plt.savefig(f"plot_dendro{suffix}.png")
    except Exception as e:
        print (show_exc(e))

    print (f"Elapsed time: {time.time() - t0:.1f} s")