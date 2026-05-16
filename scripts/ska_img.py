#!/usr/bin/env python
import warnings
warnings.simplefilter('ignore')
import os
import glob
import math
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm, Normalize, PowerNorm
from matplotlib.ticker import FormatStrFormatter
from photutils.aperture import SkyEllipticalAperture


from astropy.io import fits as pyfits
from astropy.io import fits
from astropy.wcs import WCS
from astropy import stats, units as u, constants as cte
from astropy.coordinates import SkyCoord, Angle, Galactic, ICRS, FK5, FK4


from scipy.ndimage.measurements import center_of_mass
from astrodendro import Dendrogram
from astrodendro.analysis import PPStatistic
from radio_beam import Beam
import reproject as rp

from utils import show_exc, Source, fit_lines_Npts_per_pixel
# polyfit

VERBOSE = True


# def _theta_image_from_pa_sky(center_sc: SkyCoord, pa_sky_deg: float, w: WCS) -> float:
#     """
#     Devuelve theta (rad) en el plano de imagen (x->derecha, y->ARRIBA),
#     a partir del PA astronómico (grados E de N). Usa un paso angular pequeño
#     y WCS para medir el ángulo real en píxeles.
#     """

#     if hasattr(pa_sky_deg, 'unit'):
#         pa_sky_deg = pa_sky_deg.to(u.deg).value

#     x0, y0 = w.world_to_pixel(center_sc)

#     step = 1.0 * u.arcsec  
#     p1 = center_sc.directional_offset_by(pa_sky_deg * u.deg, step)


#     x1, y1 = w.world_to_pixel(p1)


#     dx = x1 - x0
#     dy_up = -(y1 - y0)  
#     theta = np.arctan2(dy_up, dx) 
#     return theta

def arcsec2pc(omega, D):
    if not hasattr(omega,'unit'):
        omega = omega * u.arcsec
    if not hasattr(D,'unit'):
        D = D * u.kpc
    return (np.tan(omega.to(u.rad)/2) * 2 * D).to(u.pc)

def get_file(re, paths):
    for path in paths:
        if re in os.path.basename(path):
            return(path)

def mylog(*args, end='\n', flush=False, verbose=VERBOSE):
    if verbose:
        print(*args, end=end, flush=flush)

def f_p(args):
    [idx,idy,x,y,value,frame] = args
    coord = SkyCoord(x * u.deg, y*u.deg, frame=frame)
    return ([idx, idy, coord.icrs.ra.to(u.deg).value, coord.icrs.dec.to(u.deg).value, value])

def distance(m,ref, factor=(1.,1.)):
    aux = np.zeros(m.shape)
    for idy, row in enumerate(m):
        for idx,col in enumerate(row):
            aux[idy,idx]= np.sqrt( (factor[0]*(idy - ref[0]))**2 + ((factor[1] * (idx - ref[1]))**2))

    return(aux)

def distance_shape(m_shape,ref, factor=(1.,1.)):
    aux = np.zeros(m_shape)
    for idy, row in enumerate(aux):
        for idx,col in enumerate(row):
            aux[idy,idx]= np.sqrt( (factor[0]*(idy - ref[0]))**2 + ((factor[1] * (idx - ref[1]))**2))

    return(aux)

def Lineal(x, m, a):
    return (x*m + a)

class SKAImage:
    data = None
    header = None

    @classmethod
    def spectral_index(cls, list_images, mask=None, map_mode=False):
        if not map_mode:
            if len(list_images) == 2:
                img1 = list_images[0]
                img2 = list_images[1]
                freq1 = img1.restfreq().to(u.Hz).value
                freq2 = img2.restfreq().to(u.Hz).value

                if mask is not None:
                    orig_img1 = img1.data
                    img1.data[~mask] = 0
                    orig_img2 = img2.data
                    img2.data[~mask] = 0


                if freq1 is None or freq2 is None:
                    raise Exception ("Cannot compute spectral index: rest frequency not found in one of the images")

                if img1.data.shape != img2.data.shape:
                    raise Exception ("Cannot compute spectral index: images have different shapes")

                flux_ratio = img1.janskys() / img2.janskys()
                spectral_index = np.log10(flux_ratio) / np.log10(freq1 / freq2)
                if mask is not None:
                    img1.data = orig_img1
                    img2.data = orig_img2
                return spectral_index
            else:
                # Use S(nu) = S0 (nu/nu0)^alpha
                # Apply mask if provided
                if mask is not None:
                    orig_data = []
                    for img in list_images:
                        orig_data.append(img.data)
                        img.data[~mask] = 0 

                list_log_freqs = [np.log10(img.restfreq().to(u.Hz).value) for img in list_images]
                list_log_fluxes = [np.log10(img.janskys().to(u.Jy).value) for img in list_images]

                alpha, intercept = np.polyfit(list_log_freqs, list_log_fluxes, 1)
                if mask is not None:
                    for i, img in enumerate(list_images):
                        img.data = orig_data[i]
                return alpha    
        else:
            sorted_images = sorted(list_images, key=lambda x: x.restfreq())
            list_freqs = [img.restfreq().to(u.Hz).value for img in sorted_images]
            list_matrix = []
            for img in sorted_images:
                # img.data_to_pixel()
                list_matrix.append(img.data2d)
                # img.data_to_beam()
            list_matrix = np.array(list_matrix)

            alpha_map = fit_lines_Npts_per_pixel(list_matrix, list_freqs, mode="loglog", min_positive=1e-30)[0]
            header = sorted_images[0].header
            header['BUNIT'] = 'dimensionless'
            img_alpha = SKAImage(data=alpha_map, header=header)
            return img_alpha


    def __init__(self, path = None, data=None, step = None, x0 = 0, y0 = 0, unit=u.deg, interf=False, frame=None, header=None, extent=None, bunit=None, cube=False):
        if path is not None:
            if os.path.exists(path):
                self.path = path
                self.fits = pyfits.open(path)
                header = self.fits[0].header
                if "OBJECT" in header.keys():
                    self.name = header['OBJECT']
                else:
                    self.name = 'REGION'

                data = self.fits[0].data
                self.header = header

                if "IMGTYPE" in header.keys():
                    self.interf = header["IMGTYPE"] == 'INTERFEROMETER'
                    self.interf = False

                if "CUNIT1" not in header.keys():
                    header["CUNIT1"] = "deg"
                if "CUNIT2" not in header.keys():
                    header["CUNIT2"] = "deg"
                if "CDELT1" not in header.keys():
                    if "CD1_1" in header.keys():
                        header["CDELT1"] = header["CD1_1"]
                if "CDELT2" not in header.keys():
                    if "CD2_2" in header.keys():
                        header["CDELT2"] = header["CD2_2"]

                if frame is None:
                    if "GLON" in self.header['CTYPE1'].upper():
                        self.frame = Galactic
                        self.interf = False
                    else:
                        self.frame = FK5
                        self.interf= True
                else:
                    self.frame = frame

                if (self.frame != Galactic):
                    self.coords_orig = SkyCoord(ra=(header['CRVAL1'] * u.Unit(header['CUNIT1'])).to(u.deg), dec=(header['CRVAL2'] * u.Unit(header['CUNIT2'])).to(u.deg), frame=self.frame)
                    self.coords = (self.coords_orig.ra, self.coords_orig.dec, self.coords_orig.frame)
                else:
                    self.coords_orig = SkyCoord(Galactic((header['CRVAL1'] * u.Unit(header['CUNIT1'])).to(u.deg), (header['CRVAL2'] * u.Unit(header['CUNIT2'])).to(u.deg)), frame=Galactic)
                    self.coords = (self.coords_orig.l, self.coords_orig.b, self.coords_orig.frame)

                if ("BUNIT" in header.keys()):
                    self.header['BUNIT'] = header['BUNIT'].lower().replace("jy", "Jy")



                # header['NAXIS'] = 2
                self.step = ((header['CDELT1'] * u.Unit(header['CUNIT1'])).to(u.deg), (header['CDELT2'] * u.Unit(header['CUNIT2'])).to(u.deg))
                self.pixelarea = np.abs(self.step[0] * self.step[1])
                self.size = (header['NAXIS1'], header['NAXIS2'])
                #data = np.nan_to_num(data, False, 1e-20)
                self.data = data
                self.header = header
                x0 = self.coords[0].value - (header['CRPIX1']-1)*self.step[0].value
                y0 = self.coords[1].value - (header['CRPIX2']-1)*self.step[1].value

                self.extent = (x0, x0 + self.step[0].value * self.size[0], y0, y0 + self.step[1].value * self.size[1])
                self.axis=[np.arange(x0, x0 + self.step[0].value * self.size[0], self.step[0].value),np.arange(y0, y0 + self.step[1].value * self.size[1], self.step[1].value)]
            else:
                mylog ("Sorry. The file %s does not exists." % path)
                return None
        else:
            try:
                self.data = data
                self.size = data.shape
                self.interf = True
                if header is None:
                    self.step = step
                    self.extent = (x0, x0 + self.step[0].value * self.size[0], y0, y0 + self.step[1].value * self.size[1])
                    self.axis=[np.arange(x0, x0 + self.step[0].value * self.size[0], self.step[0].value),np.arange(y0, y0 + self.step[1].value * self.size[1], self.step[1].value)]
                    self.header = self.new_header()
                else:
                    # if (not self.cube):
                    #     self.header = self.fix_header(header)
                    # else:
                    #     self.header = header
                    if 'NAXIS' not in header.keys() and 'WCSAXES' in header.keys():
                        header['NAXIS'] = header['WCSAXES']
                    if 'NAXIS1' not in header.keys():
                        header['NAXIS1'] = self.data.shape[-1]
                    if 'NAXIS2' not in header.keys():
                        header['NAXIS2'] = self.data.shape[-2]
      
                    self.header = header
                    self.interf = True
                    try:
                        wcs = WCS(self.header)
                        x0, y0 = self.pix2coords(0, 0, 0)
                        x1, y1 = self.pix2coords(1, 1, 0)
                        xN, yN = self.pix2coords(self.header['NAXIS1'], self.header['NAXIS2'], 0)
                        # x0, y0 = wcs.wcs_pix2world(0, 0, 0)
                        # x1, y1 = wcs.wcs_pix2world(1, 1, 0)
                        # xN, yN = wcs.wcs_pix2world(self.header['NAXIS1'], self.header['NAXIS2'], 0)
                    except Exception as e:
                        print (show_exc(e))
                        (x0,xN,y0,yN) = (1*header['CDELT1'], self.size[0]*header['CDELT1'],1*header['CDELT2'], self.size[1]*header['CDELT2'])
                    self.extent = (x0,xN,y0,yN)
                    self.step = ((header['CDELT1'] * u.Unit(header['CUNIT1'])).to(u.deg), (header['CDELT2'] * u.Unit(header['CUNIT2'])).to(u.deg))
                    correction =  header['CDELT2']/ (x0/15 * math.cos(math.radians(180)))
                    self.pixelarea = np.abs(self.step[0] * self.step[1])
                    self.axis=[np.arange(x0, x0 + correction * self.size[0], correction),np.arange(y0, y0 + self.step[1].value * self.size[1], self.step[1].value)]
                    self.coords_orig = SkyCoord(ra=(x0 * unit).to(u.deg), dec=(y0 * unit).to(u.deg))
                    self.coords = (self.coords_orig.ra, self.coords_orig.dec, self.coords_orig.frame)
                    if "BUNIT" in header.keys():
                        self.header['BUNIT'] = header['BUNIT']
                try:
                    self.name = self.header['OBJECT']
                except:
                    self.name = 'UNSET'

                self.header['IMGTYPE'] = 'COMBINATION'
            except Exception as e:
                mylog (show_exc(e))

        if (not hasattr(self, 'frame')):
            self.frame = FK5

    @classmethod
    def new_model(cls, center, fov, pixels, freq=1420*u.MHz, beam = None):
        self = cls.__new__(cls)
        if not isinstance(center, SkyCoord):
            raise ValueError("center must be a SkyCoord object")
        
        if not hasattr(fov, 'unit'):
            if isinstance(fov, (int, float)):
                fov = fov * u.deg
            elif isinstance(fov, str):
                try:
                    fov = Angle(fov)
                except:
                    raise ValueError("fov must be an astropy Quantity with angular units (e.g. 1.5deg, 2rad, 30arcmin)")
            else:       
                raise ValueError("fov must be an astropy Quantity with angular units")
        if not fov.unit.is_equivalent(u.deg):
            raise ValueError("fov must be an angular quantity with astropy units (e.g. 1.5deg, 2rad, 30arcmin)")
        if not freq.unit.is_equivalent(u.Hz):
            raise ValueError("freq must be a frequency quantity with astropy units (e.g. 1.4GHz, 150MHz, 3THz)")
        hdu = pyfits.PrimaryHDU()
        self.header = hdu.header
        self.header['SIMPLE'] = True
        self.header['BITPIX'] = -32
        self.header['NAXIS'] = 2
        self.header['NAXIS1'] = pixels
        self.header['NAXIS2'] = pixels
        self.header['EXTEND'] = True
        self.name ="Model"
        self.header['CRVAL1'] = center.ra.to(u.deg).value
        self.header['CRVAL2'] = center.dec.to(u.deg).value
        self.header['CUNIT1'] = 'deg'
        self.header['CUNIT2'] = 'deg'
        self.header['CDELT1'] = -1*(fov.to(u.deg).value)/pixels
        self.header['CDELT2'] = (fov.to(u.deg).value)/pixels
        self.header['CRPIX1'] = pixels/2 + 0.5
        self.header['CRPIX2'] = pixels/2 + 0.5
        self.header['BUNIT'] = 'Jy/beam'
        self.header['RESTFRQ'] = freq.to(u.Hz).value
        
        if center.frame.name.lower() == "galactic":
            self.frame = Galactic
            self.header['CTYPE1'] = 'GLON-SIN'
            self.header['CTYPE2'] = 'GLAT-SIN'
            self.header['RADESYS'] = 'GALACTIC'
        else:
            self.frame = eval(center.frame.name.upper())
            self.header['CTYPE1'] = 'RA---SIN'
            self.header['CTYPE2'] = 'DEC--SIN'
            self.header['RADESYS'] = center.frame.name.upper()
        self.size = (pixels, pixels)
        self.data = np.zeros((pixels, pixels))

        if beam:
            self.set_beam(beam)
        return self




    def reproject(self, img):
        obj = WCS(img.header)
        array, footprint = rp.reproject_interp((self.data, WCS(self.header)), obj, img.data.shape)
        a =SKAImage(data=array,header=img.header)
        a.name = self.name
        a.header['AUTHOR'] ='IMGCombine by d.diaz@irya.unam.mx'
        if hasattr(img, 'frame'):
            a.frame = img.frame
        else:
            a.frame = FK5
        return a


    def set_beam(self, beam):
        try:
            self.header['BMIN'] = beam.minor.to(u.deg).value
            self.header['BMAJ'] = beam.major.to(u.deg).value
            self.header['BPA'] = beam.pa.value
        except:
            pass

    def data_to_pixel(self, verbose=False):
        if self.header['BUNIT'] == 'Jy/beam':
            mylog("\tConverting from Jy/beam to Jy/px...", flush=True, end= " ", verbose=verbose)
            self.data = self.data / (self.omega_beam()/self.omega_pix()).value #Convert Jy/beam to Jy/px
            self.header['BUNIT'] = 'Jy/px'
            self.units = "Jy/px"

    def data_to_beam(self, beam=None, verbose=False):
        if self.header['BUNIT'] != 'Jy/beam':
            mylog ("\tConverting from Jy/px to Jy/beam...", flush=True, end=" ")
            self.data *= (self.omega_beam(beam)/self.omega_pix()).value #Convert Jy/px to Jy/beam
            self.header['BUNIT'] = 'Jy/beam'
            self.units = "Jy/beam"

    @property
    def bunit(self):
        try:
            return u.Unit(self.header['BUNIT'])
        except Exception as e:
            mylog(show_exc(e))
            return u.Jy / u.beam
        return u.Jy / u.beam

    def restfreq (self, unit=u.GHz, default=None):
        try:
            units_in_header = u.Unit(self.wcs.wcs.cunit[2])
            if not units_in_header.is_equivalent(u.Hz):
                units_in_header = u.Hz
        except:
            units_in_header = u.Hz
        try:    
            if 'WAVELENG' in self.header.keys():
                return (self.header["WAVELENG"] * units_in_header).to(unit, equivalencies=u.spectral())
            if 'RESTFRQ' in self.header.keys():
                return (self.header['RESTFRQ'] * units_in_header).to(unit, equivalencies=u.spectral())
            elif 'REFFREQ' in self.header.keys():
                return (self.header['REFFREQ'] * units_in_header).to(unit, equivalencies=u.spectral())
            elif 'RESTFREQ' in self.header.keys():
                return (self.header['RESTFREQ'] * units_in_header).to(unit, equivalencies=u.spectral())
            elif 'CRVAL3' in self.header.keys():
                return (self.header['CRVAL3'] * units_in_header).to(unit, equivalencies=u.spectral())
            else:
                return default
        except Exception as e:
            mylog(show_exc(e))
            return default
    
    def wavelength (self, unit=u.um):
        return (self.restfreq().to(unit, equivalencies=u.spectral()))


    def pix2coords(self,x,y=None, freq=0, z = 0, skyMode=False):
        try:
            if y is None:
                (x,y) = x

            try:
                wcs = self.wcs2d
                if wcs.naxis < 2:
                    raise Exception ("WCS has less than 2 axis")
                if wcs.naxis > 2:
                    x0, y0, _, _ = wcs.wcs_pix2world(x, y, 0, 0, 0)
                else:
                    x0, y0 = wcs.wcs_pix2world(x, y, 0)
            except:
                wcs = self.fix_header(self.header)
                if wcs.naxis < 2:
                    raise Exception ("WCS has less than 2 axis")
                if wcs.naxis > 2:
                    x0, y0, _, _  = wcs.wcs_pix2world(x, y, 0, 0, 0)
                else:
                    x0, y0 = wcs.wcs_pix2world(x, y, 0)

            if skyMode:
                raUnit = u.Unit(self.header['CUNIT1'])
                decUnit = u.Unit(self.header['CUNIT2'])
                coord = SkyCoord(x0 * raUnit, y0 * decUnit, frame=self.frame)
                return coord
            else:
                return (x0,y0)
        except Exception as e:
            mylog(show_exc(e))
            raise e

    def coords2pix(self,coord):
        coord=coord.transform_to(self.frame)

        try:
            x = coord.l
            y = coord.b
        except:
            x = coord.ra
            y = coord.dec
        try:
            wcs = self.wcs
            if wcs.naxis > 2:
                x0, y0, _, _ = wcs.wcs_world2pix(x,y,0,0,0)
            else:
                x0, y0 = wcs.wcs_world2pix(x, y, 0)
        except:
            wcs = self.wcs2d
            x0, y0 = wcs.wcs_world2pix(x, y, 0)
        return (x0,y0)

    def world_to_pixel(self, coord):
        return self.coords2pix(coord)

    def pixel_to_world(self, x, y=None):
        return self.pix2coords(x,y)
    
    def get_center(self, frame=None):
        if frame == None:
            frame = self.frame
        try:
            if isinstance(frame, str):
                frame = eval(frame)
        except Exception as e:
            mylog(show_exc(e))
            frame = self.frame

        try:
            y = self.data.shape[-2]
            x = self.data.shape[-1]
            (x,y) = self.pix2coords(int(x/2), int(y/2))
            x = x * u.Unit(self.header['CUNIT1'])
            y = y * u.Unit(self.header['CUNIT2'])

            if self.frame == Galactic:
                coords = SkyCoord(l=x, b=y, frame=Galactic)
            else:
                coords = SkyCoord(ra=x, dec=y, frame=FK5)
            if self.frame != frame:
                coords = coords.transform_to(frame)
            return coords
        except Exception as e:
            mylog (show_exc(e))
            return (0,0)

    def __str__(self):
        mylog (self.name)

    def get_origin(self):
        if self.interf:
            pass
        else:
            pass

    def extent_to_coords(self, frame=None):
        (x0,x1,y0,y1) = self.extent 
        x0 = x0 * u.Unit(self.header['CUNIT1'])
        x1 = x1 * u.Unit(self.header['CUNIT1'])
        y0 = y0 * u.Unit(self.header['CUNIT2'])
        y1 = y1 * u.Unit(self.header['CUNIT2'])
        if self.frame == FK5:
            coords=(SkyCoord(ra=x0,dec=y0, frame=self.frame),SkyCoord(ra=x0,dec=y1, frame=self.frame),SkyCoord(ra=x1,dec=y1, frame=self.frame),SkyCoord(ra=x1,dec=y0, frame=self.frame))
        else:
            coords=(SkyCoord(l=x0,b=y0, frame=self.frame),SkyCoord(l=x0,b=y1, frame=self.frame),SkyCoord(l=x1,b=y1, frame=self.frame),SkyCoord(l=x1,b=y0, frame=self.frame))

        if frame == None:
            frame = self.frame
        try:
            if isinstance(frame, str):
                frame = eval(frame)
        except Exception as e:
            mylog(show_exc(e))
            frame = self.frame

        if frame != self.frame:
            return [coord.transform_to(frame) for coord in coords ]
        else:
            return [coord for coord in coords ]

    def get_beam(self, defaultBeam=None):
        try:
            beam = Beam.from_fits_header(self.header)
            return (beam)
        except Exception as e:
            if defaultBeam is not None:
                self.set_beam(defaultBeam)
                print ("Warning: Beam information not found in header. Using default beam of %s" % defaultBeam)
                return defaultBeam
            else:
                print ("Warning: Beam information not found in header and no default beam provided.")
                return None

    @property
    def beam_str(self):
        beam = self.get_beam()
        if beam is not None:
            return (f"{beam.major.to(u.arcsec).value:3.2f}\"x{beam.minor.to(u.arcsec).value:3.2f}\" @{beam.pa.to(u.deg).value:3.1f}°")
        else:
            return ("NoBeam")
        
    def zoom(self, zoom, in2d=True):
        ''' Zoom into a region of the image. 
            zoom = [ra_center, dec_center, radius] or [SkyCoord, radius]
        '''

        if in2d:
            orig_data = np.copy(self.data2d)
        else:
            orig_data = np.copy(self.data)
        if (len(zoom) == 2):
            zoom_coords, radius_zoom = zoom
        else:
            (ra_center, dec_center, radius_zoom) = zoom
            zoom_coords = SkyCoord(ra_center, dec_center, frame=self.frame, unit=(u.hourangle, u.deg))

        # if self.frame != zoom_coords.frame:
        #     zoom_coords = zoom_coords.transform_to(self.frame)
        ra_center = zoom_coords.ra.to(u.deg).value
        dec_center = zoom_coords.dec.to(u.deg).value
        radius_zoom = Angle(radius_zoom).to(u.deg).value
        # Extract pixels corresponding to the zoom limits in pixels coordinates
        # x_center, y_center = self.coords2pix(zoom_coords)
        low_left = SkyCoord(ra_center + radius_zoom, dec_center - radius_zoom, frame=self.frame, unit=u.deg)
        up_right = SkyCoord(ra_center - radius_zoom, dec_center + radius_zoom, frame=self.frame, unit=u.deg)

        wcs = self.wcs2d
        # (x0, y0) = wcs.wcs_world2pix(low_left.ra.to(u.deg).value, low_left.dec.to(u.deg).value, 0)
        (x0, y0) = self.coords2pix(low_left)
        (x1, y1) = self.coords2pix(up_right)

        #Create wcs with the zoomed region
        zoomHeader = self.header.copy()
        zoomHeader['CRPIX1'] = 1
        zoomHeader['CRPIX2'] = 1
        zoomHeader['CRVAL1'] = low_left.ra.to(u.deg).value
        zoomHeader['CRVAL2'] = low_left.dec.to(u.deg).value
        zoomHeader['NAXIS1'] = int(abs((x1 - x0)))
        zoomHeader['NAXIS2'] = int(abs((y1 - y0)))
        zoomHeader['CDELT1'] = self.header['CDELT1']
        zoomHeader['CDELT2'] = self.header['CDELT2']
        zoomHeader['CUNIT1'] = self.header['CUNIT1']
        zoomHeader['CUNIT2'] = self.header['CUNIT2']
        zoomHeader['CTYPE1'] = self.header['CTYPE1']
        zoomHeader['CTYPE2'] = self.header['CTYPE2']
        zoomHeader = self.fix_header(zoomHeader)


        
        if (in2d):
            subimg = orig_data[int(min(y0,y1)):int(max(y0,y1)), int(min(x0,x1)):int(max(x0,x1))]
        else:
            subimg = orig_data[:,:,int(min(y0,y1)):int(max(y0,y1)), int(min(x0,x1)):int(max(x0,x1))]

        zoomImg = SKAImage(data=subimg, header=zoomHeader)
        return zoomImg

    def integrate_channels_image(self, center = None, a = None, b = None, pa = None, channels = [], manual_linefree_ranges = None, mask=None):
        """
        Calcula un plano 2D del continuo integrado excluyendo líneas espectrales
        
        Parameters:
        - cubo: numpy array 3D (frecuencias, y, x)
        - sigma_umbral: umbral para detectar líneas
        - percentil: percentil para calcular el continuo (50=mediana, recomendado)
        """
        cubo = self.data[0,:,:,:]
        n_freq, ny, nx = cubo.shape
        sigma_umbral = 3.0
        percentil = 50
        
        # 1. Calcular estadísticas robustas para identificar líneas
        print("Calculando estadísticas para detectar líneas...")
        
        # Media robusta a lo largo del eje espectral
        mediana_espectral = np.median(cubo, axis=0)
        import scipy.stats as stats
        mad_espectral = stats.median_abs_deviation(cubo, axis=0, scale='normal')
        
        
        # 2. Crear máscara para excluir líneas
        mascara_sin_lineas = np.ones_like(cubo, dtype=bool)
        
        for i in range(ny):
            for j in range(nx):
                espectro = cubo[:, i, j]
                # Detectar píxeles que son líneas (outliers)
                umbral_superior = mediana_espectral[i, j] + sigma_umbral * mad_espectral[i, j]
                umbral_inferior = mediana_espectral[i, j] - sigma_umbral * mad_espectral[i, j]
                
                # Marcar como líneas los valores fuera del umbral
                mascara_lineas = (espectro > umbral_superior) | (espectro < umbral_inferior)
                mascara_sin_lineas[:, i, j] = ~mascara_lineas
        
        # 3. Calcular continuo usando percentil sobre píxeles no-linea
        print("Calculando continuo integrado...")
        continuo_2d = np.zeros((ny, nx))

        import multiprocessing as mp
        pool = mp.Pool(mp.cpu_count()-1)

        
        for i in range(ny):
            for j in range(nx):
                espectro = cubo[:, i, j]
                # mascara_validos = mascara_sin_lineas[:, i, j]
                
                # if np.sum(mascara_validos) > 0:
                #     # Usar percentil sobre los valores que no son líneas
                #     continuo_2d[i, j] = np.percentile(espectro[mascara_validos], percentil)
                # else:
                #     # Fallback: usar mediana de todo el espectro
                #     continuo_2d[i, j] = np.median(espectro)
                continuo_2d[i, j] = np.median(espectro)
    
        
        continuum_img = SKAImage(data=continuo_2d, header=self.fix_header(self.header))
        residuos = cubo - continuo_2d[None, :, :]
        resid_img = SKAImage(data=residuos, header=self.header)
        resid_img.name = "Continuum-subtracted cube"
        return continuum_img, resid_img
    
    
    def draw(self, img=None, exp=0.3, cmapstr='gray', scale='power', colorbar=False, plot=False, title=None, show_beam = True, filename=None, nan2zero = False):
        
        from astropy.wcs import WCS
        if self.data is not None:
            if nan2zero:
                self.data2d[np.isnan(self.data2d)] = 0.
         
            # Plot using matplotlib with WCS projection
            if img is None:
                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(1, 1, 1, projection=self.wcs.celestial)
            else:
                ax = img
            norm = None
            if scale == 'power':
                norm = PowerNorm(gamma=exp, clip=False)
            elif scale == 'log':
                norm = LogNorm()
            elif scale == 'linear':
                norm = Normalize()
                norm.autoscale(self.data2d)

            if colorbar:
                            # ax.imshow(image.data2d, origin='lower', cmap='viridis', interpolation='nearest',vmin=0, vmax=np.percentile(image.data2d, 99))

                im = ax.imshow(self.data2d, origin='lower', norm=norm, cmap=cmapstr, interpolation='nearest')
                cbar = plt.colorbar(im, ax=ax, extend='both', format=FormatStrFormatter('%.2e'))
                cbar.set_label(f"[{self.header['BUNIT']}]")
            else:
                im = ax.imshow(self.data2d, origin='lower', norm=norm, cmap=cmapstr, interpolation='nearest')

            if show_beam:
                try:
                    self.draw_beam(ax)
                except Exception as e:
                    print (show_exc(e))

            ax.set_xlabel('RA')
            ax.set_ylabel('Dec')
            if title is None:
                title = self.name
            ax.set_title('{}'.format(title))
            if plot:
                plt.show()
            if filename is not None:
                plt.savefig(filename)
            return ax
        return None


    def center_of_mass (self, mask=None, coords=True, sigma=0.0, verbose=False):
        mad = self.mad(verbose=verbose)
        if (mad is None):
            mad = 0.
        aux = np.copy(self.data)
        aux[np.isnan(aux)] = 0.
        
        aux[(aux < sigma * mad)] = 0
        if mask is not None:
            aux[~(mask)] = 0.
        cm = center_of_mass(aux)
        if not coords:
            return(cm)
        else:
            (ra,dec) = self.pix2coords(cm[1],cm[0])
            coord_cm = SkyCoord(ra,dec, frame=self.frame, unit=(u.deg, u.deg))
            return (coord_cm)


    def draw_beam(self, img, color="green"):
        try:
            beam = self.get_beam()
            (cx,cy) = (beam.major.to(u.arcsec).value * 0.5 / (self.omega_pix().to(u.arcsec**2).value**0.5),beam.major.to(u.arcsec).value * 0.5/ (self.omega_pix().to(u.arcsec**2).value**0.5))

            # Convert major and minor to pixel units
            major_pix =np.abs(beam.major.to(u.arcsec).value / (self.omega_pix().to(u.arcsec**2).value**0.5))
            minor_pix =np.abs(beam.minor.to(u.arcsec).value / (self.omega_pix().to(u.arcsec**2).value**0.5))
            # (cx,cy) = self.pix2coords(int(max(cx, 0.1 * self.data.shape[1])),int(max(cy, 0.1 * self.data.shape[0])))
            (cx,cy) = (major_pix +1 , major_pix + 1)
            beam_pa_value = np.sign(self.header['CDELT1']) * beam.pa.value
            if beam.minor.value != beam.major.value:
                frame = patches.Ellipse((cx,cy), width=minor_pix, height=major_pix, color=color, angle=beam_pa_value)
            else:
                frame = patches.Circle((cx,cy), major_pix, color=color)

            img.add_patch(frame)

        except Exception as e:
            mylog (show_exc(e))

    def sigma(self, value=3., noise=None):
        data = np.copy(self.data)
        try:
            mad = stats.median_absolute_deviation(data, axis=None, ignore_nan=True) 
            data = np.nan_to_num(data, False, mad)
            if noise  == None:
                noise = min(0,mad)
            data[np.where(data < value*mad)] = noise
            return data
        except Exception as e:
            mylog (show_exc(e))
            return data
        
    def angular_resolution(self, unit=u.arcsec):
        try:
            beam = self.get_beam()
            if beam is not None:
                return (beam.major.to(unit))
            else:
                return None
        except Exception as e:
            mylog(show_exc(e))
            return None
        
    def spw_range(self, unit=u.GHz):
        try:
            freq_max = None
            freq_min = self.restfreq()
            bandwidth = self.header['CDELT3'] * self.header['NAXIS3'] * u.Hz
            if freq_min is not None:
                freq_max = freq_min + bandwidth
                return (freq_min.to(unit, equivalencies=u.spectral()), freq_max.to(unit, 
                                                                                equivalencies=u.spectral()))

            if 'CRVAL3' in self.header.keys() and 'CDELT3' in self.header.keys() and 'NAXIS3' in self.header.keys():
                freq_min = (self.header['CRVAL3'] - (self.header['CDELT3'] * self.header['CRPIX3'])) * u.Hz
                freq_max = (self.header['CRVAL3'] + (self.header['CDELT3'] * (self.header['NAXIS3'] - self.header['CRPIX3']))) * u.Hz

            return (freq_min.to(unit, equivalencies=u.spectral()), freq_max.to(unit, equivalencies=u.spectral()))
        except Exception as e:
            return [0 * unit, 0 * unit]

    def show_headers(self, exclude = ['HISTORY'], onlykeys=False):
        if self.header is not None:
            for item in self.header:
                if item not in exclude:
                    if onlykeys:
                        mylog ("\t%s" % (item))
                    else:
                        mylog ("\t%s %s" % (item, self.header[item]))

    def field_show(self, frame=None):
        if self.data is None:
            mylog (self.coords)
        else:
            mylog (self.coords[0].transform_to(frame))

    def field(self):
        if self.data is not None:
            lng = Angle(self.coords[0]).dms
            lat = Angle(self.coords[1]).dms
            lng_range =  self.step[0] * self.size[0]
            lat_range =  self.step[1] * self.size[1]
            lng_end = Angle(self.coords[0] + lng_range).dms
            lat_end = Angle(self.coords[1] + lat_range).dms
            return(np.abs(lng_range * lat_range))
        return (0 * u.deg**2)

    def fov(self, unit=u.arcsec):
        if self.data is not None:
            lng = Angle(self.coords[0]) - (self.step[0] * self.header['CRPIX1'])
            lat = Angle(self.coords[1]) - (self.step[1] * self.header['CRPIX2'])
            if self.frame == Galactic:
                params = {'l':lng,'b':lat}
            else:
                params = {'ra':lng,'dec':lat}
            params['frame'] = self.frame
            coords = SkyCoord(**params)
            dx = (self.header['NAXIS1'] - self.header['CRPIX1']) * self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])
            dy = (self.header['NAXIS2'] - self.header['CRPIX2']) * self.header['CDELT2'] * u.Unit(self.header['CUNIT2'])
            lng_end =  Angle(self.coords[0]) + dx 
            lat_end =  Angle(self.coords[1]) + dy
            if self.frame == Galactic:
                params = {'l':lng_end,'b':lat_end}
            else:
                params = {'ra':lng_end,'dec':lat_end}
            params['frame'] = self.frame
            end_coords = SkyCoord(**params)
            if self.frame == Galactic:
                coords = coords.transform_to(FK5)
                end_coords = end_coords.transform_to(FK5)
            range_x = np.abs(coords.ra - end_coords.ra).to(unit)
            range_y = np.abs(coords.dec - end_coords.dec).to(unit)
            return {'x':range_x, 'y':range_y, 'area':range_x * range_y}
        return {'x':0*unit, 'y':0*unit, 'area':0*unit**2}
    
    @property
    def header2d(self):
        header2d = self.header.copy()
        header2d['NAXIS'] = 2
        header2d['WCSAXES'] = 2
        header2d['NAXIS1'] = self.header['NAXIS1']
        header2d['NAXIS2'] = self.header['NAXIS2']
        keys_to_remove = [  'NAXIS3', 'NAXIS4', 
                            'PC1_1', 'PC2_1', 'PC3_1', 'PC4_1', 'PC1_2', 'PC2_2',
                            'PC3_2', 'PC4_2', 'PC1_3', 'PC2_3', 'PC3_3', 'PC4_3',
                            'PC1_4', 'PC2_4', 'PC3_4', 'PC4_4',                            'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03', 'PC01_04', 'PC02_04',
                                'PC03_04', 'PC04_04', 'PC03_01', 'PC03_02', 'PC03_03', 'PC03_04',
                                'PC04_01', 'PC04_02', 'PC04_03', 'PC04_04',
                            'CTYPE3', 'CRVAL3', 'CDELT3', 'CRPIX3', 'CUNIT3', 'NAXIS3', 'CROTA3',
                            'CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4', 'NAXIS4', 'CROTA4',
                            'PV2_1', 'PV2_2', 'SPECSYS',
                            'ALTRVAL', 'ALTRPIX', 'VELREF', 'HISTORY',  'COMMENT']
        for key in keys_to_remove:
            if key in header2d.keys():
                del header2d[key]
        return header2d



    def fix_header(self, header):
        keys_to_remove = [  'NAXIS', 'NAXIS3', 'NAXIS4', 
                            'PC1_1', 'PC2_1', 'PC3_1', 'PC4_1', 'PC1_2', 'PC2_2',
                            'PC3_2', 'PC4_2', 'PC1_3', 'PC2_3', 'PC3_3', 'PC4_3',
                            'PC1_4', 'PC2_4', 'PC3_4', 'PC4_4',                            'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03', 'PC01_04', 'PC02_04',
                                'PC03_04', 'PC04_04', 'PC03_01', 'PC03_02', 'PC03_03', 'PC03_04',
                                'PC04_01', 'PC04_02', 'PC04_03', 'PC04_04',
                            'CTYPE3', 'CRVAL3', 'CDELT3', 'CRPIX3', 'CUNIT3', 'NAXIS3', 'CROTA3',
                            'CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4', 'NAXIS4', 'CROTA4',
                            'PV2_1', 'PV2_2', 'SPECSYS',
                            'ALTRVAL', 'ALTRPIX', 'VELREF', 'HISTORY', 'COMMENT']

        if "CUNIT1" not in header.keys():
                    header["CUNIT1"] = "deg"
        if "CUNIT2" not in header.keys():
            header["CUNIT2"] = "deg"
        if "CDELT1" not in header.keys():
            if "CD1_1" in header.keys():
                header["CDELT1"] = header["CD1_1"]
        if "CDELT2" not in header.keys():
            if "CD2_2" in header.keys():
                header["CDELT2"] = header["CD2_2"]

        hdu = pyfits.PrimaryHDU()
        if 'BUNIT' not in header.keys():
            hdu.header['BUNIT'] = 'Jy/beam'
        hdu.header['AUTHOR'] ='IMGCombine by d.diaz@irya.unam.mx'
        hdu.header['SIMPLE'] = header['SIMPLE'] if 'SIMPLE' in header.keys() else True
        hdu.header['BITPIX'] = -32
        hdu.header['EXTEND'] = header['EXTEND'] if 'EXTEND' in header.keys() else True
        for key in header:
            try:
                if key not in keys_to_remove and key not in hdu.header:
                    hdu.header[key] = header[key]
            except:
                pass
        hdu.header['NAXIS'] = 2


        if header['NAXIS'] > 2:
            hdu.header['IMGTYPE'] = 'INTERFEROMETER'
        else:
            if 'ORIGIN' in hdu.header.keys() and 'cornish' in header['ORIGIN'].lower():
                hdu.header['IMGTYPE'] = 'INTERFEROMETER'
            else:
                hdu.header['IMGTYPE'] = 'SINGLE-DISH'
        return hdu.header

    def new_header(self, hdu=None, beam=None):
        try:
            keys_to_remove = [  'PC1_1', 'PC2_1', 'PC3_1', 'PC4_1', 'PC1_2', 'PC2_2',
                                'PC3_2', 'PC4_2', 'PC1_3', 'PC2_3', 'PC3_3', 'PC4_3',
                                'PC1_4', 'PC2_4', 'PC3_4', 'PC4_4', 'CTYPE3', 'CRVAL3',
                                'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03', 'PC01_04', 'PC02_04',
                                'PC03_04', 'PC04_04', 'PC03_01', 'PC03_02', 'PC03_03', 'PC03_04',
                                'PC04_01', 'PC04_02', 'PC04_03', 'PC04_04',
                                'CDELT3', 'CRPIX3', 'CUNIT3', 'CTYPE4', 'CRVAL4', 'CDELT4',
                                'CRPIX4', 'CUNIT4', 'PV2_1', 'PV2_2', 'SPECSYS',
                                'ALTRVAL', 'ALTRPIX', 'VELREF', 'HISTORY', 'COMMENT','NAXIS3','NAXIS4']
            if not hdu:
                hdu = pyfits.PrimaryHDU()
            hdu.header['SIMPLE'] = True
            hdu.header['BITPIX'] = -32
            hdu.header['NAXIS'] = 2
            hdu.header['NAXIS1'] = self.size[0]
            hdu.header['NAXIS2'] = self.size[1]
            hdu.header['EXTEND'] = True
            if self.header:
                hdu.header['CRVAL1'] = self.header['CRVAL1']
                hdu.header['CRVAL2'] = self.header['CRVAL2']
                hdu.header['CDELT1'] = self.header['CDELT1']
                hdu.header['CDELT2'] = self.header['CDELT2']
                hdu.header['CRPIX1'] = self.header['CRPIX1']
                hdu.header['CRPIX2'] = self.header['CRPIX2']
            else:
                hdu.header['CRVAL1'] = 0
                hdu.header['CRVAL2'] = 0
                try:
                    hdu.header['CDELT1'] = self.step[0].value
                    hdu.header['CDELT2'] = self.step[1].value
                except:
                    hdu.header['CDELT1'] = 1
                    hdu.header['CDELT2'] = 1
                hdu.header['CRPIX1'] = 0
                hdu.header['CRPIX2'] = 0
            try:
                hdu.header['CUNIT1'] = self.step[0].unit.name
                hdu.header['CUNIT2'] = self.step[1].unit.name
            except:
                hdu.header['CUNIT1'] = 'deg'
                hdu.header['CUNIT2'] = 'deg'
            if not beam:
                beam = Beam(10 * u.Unit(hdu.header['CUNIT1']),10 * u.Unit(hdu.header['CUNIT2']),0 * u.deg)
            hdu.header['BMAJ'] = beam.major.value
            hdu.header['BMIN'] = beam.minor.value
            hdu.header['BPA'] = beam.pa.value
            if not self.interf:
                hdu.header['CTYPE1'] = 'GLON-CAR'
                hdu.header['CTYPE2'] = 'GLAT-CAR'
            else:
                hdu.header['CTYPE1'] = 'RA---SIN'
                hdu.header['CTYPE2'] = 'DEC--SIN'
                hdu.header['RADESYS'] = 'FK5'
            hdu.header['BUNIT'] = 'Jy/beam'
            hdu.header['AUTHOR'] ='IMGCombine by d.diaz@irya.unam.mx'

            if self.header:
                for key in self.header:
                    if key not in keys_to_remove and key not in hdu.header:
                        try:
                            hdu.header[key] = self.header[key]
                        except Exception as e:
                            mylog (key, show_exc(e))

            if self.interf:
                hdu.header['IMGTYPE'] = 'INTERFEROMETER'
            else:
                hdu.header['IMGTYPE'] = 'SINGLE-DISH'

            return hdu.header
        except Exception as e:
            mylog (show_exc(e))
            return None

    def tofits(self, filename=None, overwrite=True):
        hdu = pyfits.PrimaryHDU()
        hdu.header = self.header
        hdu.data = self.data
        hdu.header['ORIGIN'] = 'SKA-Image'
        hdu.header['OBJECT'] = self.name

        version = 0
        if filename is None:
            filename = 'v%d_%s' % (version,os.path.basename(self.path))
            while os.path.exists(filename):
                version +=1
                filename = 'v%d_%s' % (version,os.path.basename(self.path))
        hdu.writeto(filename, overwrite=overwrite)
        self.path = filename
        self.fits = [hdu]

    def get_maximum(self):
        data = np.copy(self.data)
        if (np.isnan(data).any()):
            np.nan_to_num(data, False, -1e-20)
        id_max = np.argmax(data)
        row = int(id_max / data.shape[1])
        col = id_max % data.shape[1]
        obj = WCS(self.header)
        x,y = obj.wcs_pix2world(col, row, 0)
        x *=  u.Unit(self.header['CUNIT1'])
        y *=  u.Unit(self.header['CUNIT2'])
        coords = SkyCoord(ra=x, dec=y, frame=FK5)
        return (coords, data[row,col], [row,col])

    def janskys_beam(self):
        return self.janskys() / self.Nbeams()

    def Nbeams(self):
       return (self.Npts().value * self.omega_pix() / self.omega_beam()) * u.beam

    def janskys_sum(self):
        return np.nansum(self.data) * (u.Jy/u.beam)

    def janskys_pixel(self):
        omega_pix = self.header['CDELT1'] * u.Unit(self.header['CUNIT1']) * self.header['CDELT2'] * u.Unit(self.header['CUNIT2'])
        u_px = u.def_unit('px', omega_pix)
        npts = (self.data.size - np.isnan(self.data).sum()) * u_px
        return self.janskys() /npts


    def janskys(self, freq=None, alpha_spec=3.5, center = None, a=None, b=None, pa=None, channels=[], sigma_clip=None):
        ''' Compute total flux in the image (or in a region if center is given)
            freq: frequency to scale the flux (if different from the rest frequency)
            alpha_spec: spectral index to scale the flux (default=3.5)
            center: SkyCoord of the center of the region to compute the flux (if None, the whole image is used)
            a: semi-major axis of the region (in pixels or Quantity with angular units)
            b: semi-minor axis of the region (in pixels or Quantity with angular units)
            pa: position angle of the region (in degrees, measured from North to East)
        '''
        if freq is not None:
            freq_factor = (freq.to(u.Hz).value / self.restfreq().to(u.Hz).value)**(alpha_spec)
        else:
            freq_factor = 1.
        mask = None

        if center is not None:
            data = self.data2d.copy()

            if a is None:
                a = data.shape[1]/2.
            if b is None:
                b = a
            if pa is None:
                pa = 0.
            else:
                pa = pa + 90 * u.deg
            # pa = pa + 90 * u.deg
            (x0,y0) = self.coords2pix(center)
            y, x = np.indices(data.shape)
            xp = (x - x0) * np.cos(np.radians(pa)) + (y - y0) * np.sin(np.radians(pa))
            yp = -(x - x0) * np.sin(np.radians(pa)) + (y - y0) * np.cos(np.radians(pa))
            if (isinstance(a, u.Quantity)):
                a = (a.to(u.Unit(self.header['CUNIT1'])) / (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])))
            if (isinstance(b, u.Quantity)):
                b = (b.to(u.Unit(self.header['CUNIT2'])) / (self.header['CDELT2'] * u.Unit(self.header['CUNIT2'])))



            mask = ((xp/a)**2 + (yp/b)**2) <= 1.

            wcs = WCS(self.header)
            if wcs.naxis == 2:
                data[~(mask)] = np.nan
            else:
                data = self.data.copy()
                data[:,:,~(mask)] = np.nan
        else:
            data = self.data.copy()




        beam = self.get_beam()
        total = 0.
        wcs = WCS(self.header)
        if wcs.naxis == 2:
            if sigma_clip is not None:
                mad = self.mad()
                if mask is not None:
                    mask = mask & (data > sigma_clip * mad)
                else:
                    mask = (data > sigma_clip * mad)
                data[~(mask)] = np.nan

            npts = data[~(np.isnan(data))].size
            mean = np.nanmean(data)
            omega_beam = beam.sr.to(u.deg**2)
            omega_pix = (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])) ** 2
            full_area = omega_pix * npts
            n_beams = full_area / omega_beam
            total = mean*n_beams * u.Jy * freq_factor

            return total
        else:

            data_cube = data.copy()
            
            if channels == []:
                channels = range(self.data.shape[-3])
            for ch in channels:
                data = data_cube[0,ch,:,:].copy()
                if sigma_clip is not None:
                    mad = self.mad(channel=ch)
                    if mask is not None:
                        mask = mask & (data > sigma_clip * mad)
                    else:
                        mask = (data > sigma_clip * mad)
                    data[~(mask)] = np.nan

                npts = data[~(np.isnan(data))].size
                mean = np.nanmean(data)
                omega_beam = beam.sr.to(u.deg**2)
                omega_pix = (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])) ** 2
                full_area = omega_pix * npts
                n_beams = full_area / omega_beam
                total += mean*n_beams * u.Jy * freq_factor
            return total
        
    def neg_janskys(self):
        beam = self.get_beam()
        npts = np.count_nonzero(self.data < 0.)
        mean = np.nanmean(self.data[self.data < 0.])
        omega_beam = beam.sr.to(u.deg**2)
        omega_pix = (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])) ** 2
        full_area = omega_pix * npts
        n_beams = full_area / omega_beam
        return mean*n_beams * u.Jy

    def mean(self):
        return(np.nanmean(self.data))

    def Npts(self):
        omega_pix = (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])) ** 2
        u_px = u.def_unit('px', omega_pix)
        return (self.data.size - np.isnan(self.data).sum()) * u_px

    def detect_sources(self, threshold=5., step=1, ratio_beam=1, mask=None, verbose=False):
        try:
            wcs = self.wcs2d
            pixel_area = wcs.proj_plane_pixel_area()
            beam = self.get_beam()

            pixels_in_beam = math.ceil((beam.sr).to(u.arcsec**2) / (pixel_area.to(u.arcsec**2)))
            rms = self.mad(verbose=verbose)
            d = Dendrogram.compute(self.data2d, min_value=threshold * rms, min_delta=step * rms, min_npix=pixels_in_beam * ratio_beam, wcs=wcs)
            leaves = d.leaves
            sources = []
            metadata = {}
            metadata['data_unit'] = self.header['BUNIT'] if 'BUNIT' in self.header.keys() else 'Jy/beam'
            metadata['data_unit'] = u.Unit(metadata['data_unit'])
            metadata['spatial_scale'] = pixel_area**0.5
            metadata['beam_major'] = beam.major
            metadata['beam_minor'] = beam.minor
            metadata['beam_pa'] = beam.pa
            metadata['wcs'] = wcs


            for i, leaf in enumerate(leaves):
                source = PPStatistic(leaf, metadata=metadata)
                sources.append(DendroSource(leaf, source, wcs))
            return sources

        except Exception as e:
            mylog (show_exc(e), verbose=verbose)
            return []

    @property
    def velocity(self):
        if 'ALTRVAL' in self.header.keys():
            return (self.header['ALTRVAL'] * u.m/u.s)
        else:
            return 0. * u.m/u.s
        
    @property   
    def redshift(self):
        return (self.velocity.to(u.km/u.s) / cte.c.to(u.km/u.s)).value

    def adjacent_coords(self, coord, pixels_distance=1):
        (x,y) = self.coords2pix(coord)
        coords = []
        for dx in range(-pixels_distance, pixels_distance+1):
            for dy in range(-pixels_distance, pixels_distance+1):
                if (dx != 0 or dy != 0):
                    (ra,dec) = self.pix2coords(x+dx,y+dy)
                    coord_adj = SkyCoord(ra,dec, frame=self.frame, unit=(u.deg, u.deg))
                    coords.append(coord_adj)
        return coords
    
    def convolve_to_beam(self, beam):
        from astropy.convolution import convolve_fft
        from radio_beam import Beam
        import warnings
        try:
            warnings.simplefilter('ignore')
            current_beam = self.get_beam()
 
            if current_beam is None:
                current_beam = beam
                kernel_array = current_beam.as_kernel(self.omega_pix()**0.5)
                new_data = np.copy(self.data2d)
                new_data = convolve_fft(self.data2d, kernel_array, normalize_kernel=True, allow_huge=True)          
            elif beam > current_beam:
                kernel = beam.deconvolve(current_beam)
                kernel_array = kernel.as_kernel(self.omega_pix()**0.5)
                new_data = np.copy(self.data2d)
                new_data = convolve_fft(self.data2d, kernel_array, normalize_kernel=True, allow_huge=True)
            else:
                new_data = np.copy(self.data2d)
            new_image = SKAImage(data=new_data, header=self.header2d)
            new_image.name = self.name + ' Convolved'
            new_image.set_beam(beam)
            return new_image
        except Exception as e:
            mylog (show_exc(e))

        return None

    def get_unit(self):
        if 'BUNIT' in self.header.keys():
            return u.Unit(self.header['BUNIT'])
        else:
            return u.Jy / u.beam

        

    def mad(self, iters=100, limit=0.99, stack=False, cutoff=2., mask=None, channel=None, data=None, verbose=False, is_model=False):
        import warnings
        try:
            warnings.simplefilter('ignore')
            counter = 0
            if data is None:
                if channel is not None:
                    data = np.copy(self.data[0,channel,:,:])
                else:
                    if len(self.data.shape) < 4:
                        data = np.copy(self.data2d)
                    else:
                        data = self.data2d
            try:
                if mask is not None:
                    data[mask] = np.nan
            except Exception as e:
                mylog ("Warning in combine.mad !!!! We can not masked the image", verbose=verbose)
            if is_model:
                data = data[((data > 0.) & ~(np.isnan(data)))]
            next_mad = stats.mad_std(data, axis=None, ignore_nan=True)
            mad = 1e3
            mad_stack = []
            mad_mask = (data < cutoff * mad)
            while ((next_mad / mad < limit and counter < iters) or (counter == 0)) and (mad_mask.sum() > 0):
                mad = next_mad
                next_mad = stats.mad_std(data[mad_mask], axis=None, ignore_nan=True)
                counter += 1
                mad_stack.append(mad)
                if next_mad <= 0:
                    break
                mad_mask = (data < cutoff * mad)

            if stack: 
                return(mad_stack)# * (u.Jy/beam)
            else:
                return(mad_stack[-1])# * (u.Jy/beam)
        except Exception as e:
            mylog(f"Iteration: {counter} :> {show_exc(e)}", verbose=verbose)
            return self.data2d.std()# * (u.Jy/beam)

    def madN(self, N=1, mad=None):
        if mad is None:
            mad = self.mad()
        if N == 1:
            return mad
        data = np.copy(self.data)
        mask = (data < N * mad)
        data[~(mask)] = np.nan 
        return  (stats.mad_std(data, axis=None, ignore_nan=True))

    def mad_mask(self, iters=100, limit=0.99, stack=False, cutoff=2., mask=None):
        import warnings
        try:
            warnings.simplefilter('ignore')
            data = np.copy(self.data)
            try:
                if mask is not None:
                    data[mask] = np.nan
            except Exception as e:
                mylog ("Warning in combine.mad !!!! We can not masked the image")
            counter = 0
            next_mad = stats.mad_std(data, axis=None, ignore_nan=True)
            mad = 1e3
            mad_stack = []
            mask = (data > np.min(data))
            while (next_mad / mad < limit and counter < iters) or (counter == 0):
                mad = next_mad
                mask = (data < cutoff * mad)
                next_mad = stats.mad_std(data[np.where(data < cutoff * mad)], axis=None, ignore_nan=True)
                counter += 1
                mad_stack.append(mad)
                if next_mad <= 0:
                    break

            if stack: 
                return(mad_stack, mask)# * (u.Jy/beam)
            else:
                return(mad_stack[-1], mask)# * (u.Jy/beam)
        except Exception as e:
            mylog(show_exc(e), verbose=True)

    def max(self):
        return (np.nanmax(self.data))

    def min(self):
        return (np.nanmin(self.data))

    def dr(self, iters=100):
        return (self.max() / self.mad(iters))

    def dr_neg(self, iters=100):
        return np.abs(self.min() / self.mad(iters))

    def poso_neg(self, mad=None, iters=100, sigma=-3):
        if not mad:
            mad = self.mad(iters)
        pixels = np.where(self.data < (sigma * mad))
        n = len(pixels[0])
        return (n/self.Npts().value)

    def lo_scale(self):
        if 'LO_SCALE' in self.header.keys():
            return (self.header['LO_SCALE'])
        else:
            return 1.

    def hi_scale(self):
        if 'HI_SCALE' in self.header.keys():
            return(self.header['HI_SCALE'])
        else:
            return 1.

    def omega_pix(self):
        try:
            return (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])) ** 2
        except:
            return (1 * u.Unit('deg')) **2

    def omega_beam (self, beam=None):
        if beam is None:
            beam = self.get_beam()
        omega_beam = beam.sr.to(u.deg**2)
        return omega_beam

 

    def mask_pb(self, primary_beam_image, limit=0.2):
        mask = ((primary_beam_image.data < limit) | (np.isnan(primary_beam_image.data)))
        try:
            self.data[mask] = np.nan
        except Exception as e:
            mylog ("Warning in combine.mask_pb!!!! We can not masked the image")

  

    def __unicode__(self):
        return ("%s %s" % (self.path, self.name))
    
    def draw_source(self, ax, coord, a=None, b = None, pa = None, color="green", ellipse=False, alpha=0.3, label = None, verbose=False, lw=1):
        if isinstance(coord, Source):
            src = coord
            coord = src.coord
            if a is None:
                a = src.major_axis / 2.
            if b is None:
                b = src.minor_axis / 2.
            if pa is None:
                pa = src.pa
    
        if verbose:
            mylog(f"Draw source at {coord.to_string('hmsdms')} a={a}, b={b}, pa={pa}, color={color}, ellipse={ellipse}, alpha={alpha}, label={label}")

        try:
            if b is None:
                b = a
            if pa is None:
                pa = 0. * u.deg
            
            
            
            (cx,cy) = coord.ra.value, coord.dec.value
            ax.plot(cx, cy, marker='+', color=color, markersize=15, transform=ax.get_transform('world'))
            if label is not None:
                ax.text(cx, cy, label, color=color, fontsize=12, transform=ax.get_transform('world'))

            if ellipse and a is not None:
                try:
                    if b is None:
                        b = a
                    if pa is None:
                        pa = 0. * u.deg
                    reg = SkyEllipticalAperture(coord, a, b, pa)
                    pixels = reg.to_pixel(self.wcs2d)
                    frame = pixels.plot(ax=ax, lw=lw, color=color, alpha=alpha, fill=False)
                    # ax.add_patch(frame)



                except Exception as e:
                    print(show_exc(e))




        except Exception as e:
            mylog (show_exc(e))
    
    @property
    def wcs(self):
        try:
            wcs = WCS(self.header)
            return wcs
        except Exception as e:
            print (show_exc(e))
            return None

    @property
    def wcs2d(self):
        try:
            wcs = WCS(self.header)
            if wcs.naxis >= 2:
                return WCS(self.fix_header(self.header))
            else:
                return wcs
        except Exception as e:
            print (show_exc(e))
            return None

    @property
    def data2d(self):
        if self.data is not None:
            if (self.wcs.naxis > 2):
                return self.data[0,0]
            else:
                return self.data
        return None

    def mask_outside_ellipse_fits(self, 
        center,                 # (x0_pix, y0_pix) in pixels OR SkyCoord in the sky
        a, b,                   # major and minor axes: in pixels (float) or angular units (Quantity, e.g. 30*u.arcsec)
        pa,                     # ellipse angle
        pa_kind="image",        # "image" (from +x to +y, CCW) or "sky" (astronomical PA: E of N, degrees)
        fill_value=np.nan,      # value to set outside the ellipse
    ):
        ax = plt.subplot(1,1,1)
        center  = self.coords2pix(center) if isinstance(center, SkyCoord) else center
        a       = np.abs( a.to(u.deg).value / self.header['CDELT1'] if isinstance(a, u.Quantity) else a )
        b       = np.abs( b.to(u.deg).value / self.header['CDELT1'] if isinstance(b, u.Quantity) else b )

        ellipse = patches.Ellipse(center, width=a, height=b, angle=pa.value + 90)
        yy, xx = np.mgrid[:self.data.shape[-2], :self.data.shape[-1]]


        coords = np.vstack((xx.ravel(), yy.ravel())).T  # (N, 2) array of (x, y) pairs
        mask_flat = ellipse.contains_points(coords)
        mask = mask_flat.reshape(self.data.shape[-2], self.data.shape[-1])

        masked = np.where(mask, self.data2d, fill_value)
        img = SKAImage(data=masked, header=self.fix_header(self.header))
        return img

    def contain_coord(self, coord):
        (x,y) = self.coords2pix(coord)
        if x >= 0 and x < self.data.shape[-1] and y >= 0 and y < self.data.shape[-2]:
            return True
        return False

class DendroSource:
    def __init__(self, structure, stats, wcs=None):
        self.structure = structure
        self.stats = stats
        self.wcs = wcs
    
    @property
    def flux(self):
        return self.stats.flux
    
    @property
    def coord(self):
        return SkyCoord(ra=self.stats.x_cen * u.deg, dec=self.stats.y_cen * u.deg, frame='icrs')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('paths', nargs='*')
    parser.add_argument('--headers', action='store_true', default=False)
    parser.add_argument('--fields', action='store_true', default=False)
    parser.add_argument('--image', action='store_true', default=False)
    parser.add_argument('--exp', type=float, default=0.3)
    args = parser.parse_args()
    if len(args.paths) == 0:
        fits = glob.glob('*.fits')
        fits = sorted(fits)
    else:
        fits = args.paths
    print ("Loading %d files" % len(fits))

    images = []
    for idx, fpath in enumerate(fits):
        images.append(SKAImage(fpath))

    if args.image:
        number_of_subplots=len(fits)
        img=images[0]
        wcs = WCS(img.header)
        img.draw(plt.subplot(1,1,1, projection=wcs), exp=args.exp, cmapstr='afmhot', colorbar=True)
        plt.show()
    else:
        for img in images:
            if args.headers:
                print (f'FILE: {fpath}')
                img.show_headers(onlykeys=False)
            if args.fields:
                img.field_show()
