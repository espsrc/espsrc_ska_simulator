#!/usr/bin/env python
import warnings
warnings.simplefilter('ignore')
import astropy.io.fits as pyfits
from astropy.modeling import models, fitting
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import glob, os
from matplotlib.colors import LogNorm, Normalize, PowerNorm
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from scipy.ndimage.measurements import center_of_mass
from pylab import *
import math
import requests
from commons import show_exc
import aplpy


import argparse
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle, Galactic, ICRS, FK5, FK4
from matplotlib.ticker import FormatStrFormatter
from scipy.ndimage.interpolation import rotate
from scipy.optimize import curve_fit
#from scipy.stats import chisquare
#from scipy import stats
import matplotlib.patches as patches
import matplotlib.pylab as pylab
from multiprocessing import Pool
from astropy.wcs import WCS
import reproject as rp
from astropy import constants, units as u, table, stats, coordinates, wcs, log, coordinates as coord, convolution, modeling, visualization; 
from astropy.io import fits
from radio_beam import Beam
import radio_beam.utils as rbu
from astropy.convolution import convolve_fft
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from matplotlib.ticker import FormatStrFormatter
from uvcombine.uvcombine import feather_kernel, fftmerge #, feather_compare, pbcorr, flux_match, match_flux_units, file_in



VERBOSE = True


# from uvcombine.uvcombine import feather_kernel, fftmerge #, feather_compare, pbcorr, flux_match, match_flux_units, file_in

# try:
#     import cv2
# except Exception as e:
#     print(str(e))

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

class BGPS: #BOLOCAM
    freqs = [(1.12 * u.mm).to(u.GHz, equivalencies=u.spectral()), (2.1 * u.mm).to(u.GHz, equivalencies=u.spectral())]
    labels = ['1.1mm','2.1mm']
    bmaj = 9.166670000000E-3 * u.deg
    bmin = 9.166670000000E-3 * u.deg
    bpa = 0 * u.deg

class MGPS90:
    freqs = [90.19] * u.GHz #Revisar el paper y sacar bien las freq
    labels = ['3mm']

class ALMAIMF:
    freqs = [216.2, 93.173402] * u.GHz
    labels = ['B6','B3']

class CORNISH:
    freqs = [5] * u.GHz
    beam = Beam((1.5 * u.arcsec).to(u.deg),(1.2 * u.arcsec).to(u.deg), 0 * u.deg)
    repository = "https://cornish.leeds.ac.uk/public/data_tls"

    @staticmethod
    def get_beam_from_file(img):
        from io import StringIO
        header = img.header


        if ("OBJECT" in header.keys()):
            filepath = os.path.join(CORNISH.repository, '{}.{}'.format(header['OBJECT'].replace(' ','_'), 'fits'))
            mylog ('Searching for {}.{}'.format(header['OBJECT'].replace(' ','_'), 'fits'))
            import subprocess
            tofile = os.path.join(os.path.realpath(os.path.dirname(img.path)),'{}.{}'.format(header['OBJECT'].replace(' ','_'), 'fits'))
            if not os.path.exists(tofile):
                mylog ("\tDownloading...")
                subprocess.run(["curl", "--cipher", 'DEFAULT:!DH', "--insecure", filepath, '-o', tofile])
                mylog ("\tDONE")
            beam = Beam.from_fits_header(fits.getheader(tofile))
            CORNISH.beam = beam

        return CORNISH.beam


class Image:
    data = None
    header = None

    def __init__(self, path = None, data=None, step = None, x0 = 0, y0 = 0, unit=u.deg, interf=False, frame=None, header=None, extent=None, bunit=None):
        if path is not None:
            if os.path.exists(path):
                self.path = path
                self.fits = pyfits.open(path)
                header = self.fits[0].header
                header = self.fix_header(header)
                self.header = header
                try:
                    if 'BMIN' not in header.keys():
                        if 'ORIGIN' in header.keys() and 'cornish' in header['ORIGIN'].lower():
                            beam = CORNISH.get_beam_from_file(self)
                            header['BMAJ'] = beam.major.value
                            header['BMIN'] = beam.minor.value
                            header['BPA'] = beam.pa.value
                        if 'BGPSVERS' in self.header.keys():
                            self.header['BMAJ'] = BGPS.bmaj.value
                            self.header['BMIN'] = BGPS.bmin.value
                except Exception as e:
                    mylog (show_exc(e))
                header['NAXIS'] = 2
                if "OBJECT" in header.keys():
                    self.name = header['OBJECT']
                else:
                    self.name = 'REGION'
                if (len(self.fits[0].data.shape) > 2):
                    data = self.fits[0].data[0,0]
                    header['NAXIS'] = 2
                    #self.interf = True
                else:
                    data = self.fits[0].data
                    #self.interf = False

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



                header['NAXIS'] = 2
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
                    self.header = self.fix_header(header)
                    self.interf = True
                    try:
                        wcs = WCS(self.header)
                        x0, y0 = wcs.wcs_pix2world(0, 0, 0)
                        x1, y1 = wcs.wcs_pix2world(1, 1, 0)
                        xN, yN = wcs.wcs_pix2world(self.header['NAXIS1'], self.header['NAXIS2'], 0)
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

    def reproject(self, img):
        obj = WCS(img.header)
        array, footprint = rp.reproject_interp((self.data, WCS(self.header)), obj, img.data.shape)
        a = Image(data=array,header=img.header)
        a.name = self.name
        a.header['AUTHOR'] ='IMGCombine by d.diaz@irya.unam.mx'
        if hasattr(img, 'frame'):
            a.frame = img.frame
        else:
            a.frame = FK5
        return a

    def zoom(self):
        img = self
        (y,x) = img.data.shape
        rows = []
        cols = []

        for i in range(y):
            if (not (np.isnan(img.data[i,:]).all())):
                rows.append(i)
        for j in range(x):
            if (not (np.isnan(img.data[:,j]).all())):
                cols.append(j)

        img.header['CRPIX1'] -= min(cols)
        img.header['CRPIX2'] -= min(rows)
        img.data = img.data[min(rows):max(rows),min(cols):max(cols)]
        img.header['NAXIS1'] = img.data.shape[1]
        img.header['NAXIS2'] = img.data.shape[0]
        return (img)

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

    def data_to_beam(self, verbose=False):
        if self.header['BUNIT'] != 'Jy/beam':
            mylog ("\tConverting from Jy/px to Jy/beam...", flush=True, end=" ")
            self.data *= (self.omega_beam()/self.omega_pix()).value #Convert Jy/px to Jy/beam
            self.header['BUNIT'] = 'Jy/beam'
            self.units = "Jy/beam"


    def restfreq (self, unit=u.GHz):
        if 'WAVELENG' in self.header.keys():
            return (1.12 * u.mm).to(u.Hz, equivalencies=u.spectral())
        if 'RESTFRQ' in self.header.keys():
            return (self.header['RESTFRQ'] * u.Hz).to(unit)
        elif 'REFFREQ' in self.header.keys():
            return (self.header['REFFREQ'] * u.Hz).to(unit)
        elif 'ORIGIN' in self.header.keys() :
            if 'cornish' in self.header['ORIGIN']:
                return CORNISH.freqs[0].to(u.Hz).to(unit)
        elif 'BGPSVERS' in self.header.keys():
            return BGPS.freqs[0].to(u.Hz).to(unit)

    def pix2coords(self,x,y=None):
        if y is None:
            (x,y) = x
        try:
            wcs = WCS(self.header)
            x0, y0 = wcs.wcs_pix2world(x, y, 0)
        except:
            wcs = WCS(self.fix_header(self.header))
            x0, y0 = wcs.wcs_pix2world(x, y, 0)
        return (x0,y0)

    def coords2pix(self,coord):
        coord=coord.transform_to(self.frame)

        try:
            x = coord.l
            y = coord.b
        except:
            x = coord.ra
            y = coord.dec
        try:
            wcs = WCS(self.header)
            x0, y0 = wcs.wcs_world2pix(x, y, 0)
        except:
            wcs = WCS(self.fix_header(self.header))
            x0, y0 = wcs.wcs_world2pix(x, y, 0)
        return (x0,y0)

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
            (y,x) = self.data.shape
            (x,y) = self.pix2coords(x/2, y/2)
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

    def get_beam(self):
        try:
            try:
                beam = Beam.from_fits_header(self.header)
                return (beam)
            except:
                if 'BGPSVERS' in self.header.keys():
                     return(Beam(BGPS.bmaj,BGPS.bmin))
        except Exception as e:
            mylog(show_exc(e), verbose=True, flush=True)
            return None
        return None

    def draw(self, img=None, exp=1., rotated = True, plot=False, rectangle = None, alpha=1, contour=False, return_data=False, barcolor=True, reversecolor=True, cmapstr='viridis', scale='power', linestyles='solid', levels = [-5,5,10,20,40,80,160,320], title=None, draw_axes=True, clip=99.9, barcolor_units='Jy/beam', show_beam=True, subplot=(1,1,1), vmin=None, vmax=None, cmap=None, hidex=False, hidey=False, fontsize='xx-large', zoom=None, nancolor=None, contour_colors=None, scalebar=None, barcolor_discrete=False):
        from astropy.wcs import WCS
        if self.data is not None:
            if img is None:
                fig = plt.figure()
            else:
                fig = img

            tmp_path = "{}.fits".format(time.time())
            self.tofits(tmp_path)
            f1 = aplpy.FITSFigure(tmp_path, figure=fig, subplot=subplot, alpha=alpha)
            if cmap is not None:
                cmapstr = cmap
            if not contour:
                f1.show_colorscale(cmap=cmapstr, exponent=exp, stretch=scale, vmin=vmin, vmax=vmax)
            else:
                try:
                    colors = contour_colors
                    if levels is not None:
                        levels = levels * self.mad()
                    f1.show_colorscale(cmap=cmapstr, exponent=exp, stretch=scale, vmin=vmin, vmax=vmax)
                    f1.show_contour(tmp_path, colors=colors, levels=levels, stretch=scale)
                except Exception as e:
                    print (show_exc(e))

            if zoom is not None:
                center = self.get_center()
                fov = self.fov()
                f1.recenter(zoom[0].ra.to(u.deg), zoom[0].dec.to(u.deg), width=zoom[1] * fov['x'].to(u.deg).value, height=zoom[2] * fov['y'].to(u.deg).value)

            if hidex:
                f1.tick_labels.hide_x()
                f1.axis_labels.hide_x()
            if hidey:
                f1.tick_labels.hide_y()
                f1.axis_labels.hide_y()

            if barcolor:
                f1.add_colorbar()
                f1.colorbar.show()
                f1.colorbar.set_location('right')
                f1.colorbar.set_font(size=fontsize)
                if barcolor_units is not None:
                    f1.colorbar.set_axis_label_text(f"[{barcolor_units}]")
                else:
                    f1.colorbar.set_axis_label_text(f"[{self.header['BUNIT']}]")
                if barcolor_discrete:
                    f1.colorbar.set_style('discrete')
                f1.colorbar.set_axis_label_font(size=fontsize) #weight='bold'
            if scalebar is not None:
                f1.add_scalebar(scalebar.to(u.deg).value)
                f1.scalebar.set_font(size=fontsize)
                f1.scalebar.set_linestyle(linestyles)
                f1.scalebar.set_color('black')
                f1.scalebar.set_label(f"{scalebar[1]:.1f} {scalebar[2]}")
                f1.scalebar.set_linewidth(2)

            if title is None:
                title = self.name
            f1.set_title('{}'.format(title))
            if nancolor:
                f1.set_nan_color(nancolor)

            if show_beam:
                try: 
                    f1.add_beam()
                except:
                    pass

            if plot:
                plt.show()


            try:
                os.remove(tmp_path)
            except Exception as e:
                print (show_exc(e))
            return f1

    def center_of_mass (self, mask=None, coords=True):
        aux = np.copy(self.data)
        aux[np.isnan(aux)] = 0.
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
            (cx,cy) = self.pix2coords(int(max(cx, 0.1 * self.data.shape[1])),int(max(cy, 0.1 * self.data.shape[0])))
            if beam.minor.value != beam.major.value:
                try:
                    frame = patches.Ellipse((cx,cy), width=beam.minor.value, height=beam.major.value, color=color, angle=beam.pa.value, transform=img.get_transform('world'))
                except:
                    frame = patches.Ellipse((cx,cy), width=beam.minor.value, height=beam.major.value, color=color, angle=beam.pa.value)
            else:
                try:
                    frame = patches.Circle((cx,cy), beam.major.value, color=color, angle=beam.pa.value, transform=img.get_transform('world'))
                except:
                    frame = patches.Circle((cx,cy), beam.minor.value, color=color, angle=beam.pa.value)
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

    def to_fourier_space(self, min_beam_fraction=0.1, beam=None, pixscale=None, weights=None, deconvolve=True):
        if weights is not None:
            if not weights.shape == self.data.shape:
                raise ValueError("weights must be an array with the same shape as"
                                 " the high-res data.")
        if beam is None:
            beam = self.get_beam()

        nax2,nax1 = self.data.shape
        lowresfwhm = beam.major
        pixscale = (self.omega_pix()**0.5).to(u.deg)

        im_low = self.data.copy()/ (self.omega_beam() / self.omega_pix()).value #Convert Jy/beam to Jy/px
        #im_low = self.data
    
        if weights is None:
            kfft, ikfft = feather_kernel(nax2, nax1, lowresfwhm, pixscale)
            kfft = np.fft.fftshift(kfft)
            ikfft = np.fft.fftshift(ikfft)
        else:
            kfft = np.abs(np.fft.fftshift(np.fft.fft2(weights)))
            kfft /= np.max(kfft)
            ikfft = 1 - kfft

        yy,xx = np.indices([nax2, nax1])
        rr = ((xx-(nax1-1)/2.)**2 + (yy-(nax2-1)/2.)**2)**0.5
        angscales = nax1/rr * pixscale 

        fft_lo = np.fft.fftshift(np.fft.fft2(np.nan_to_num(im_low)))

        if deconvolve:
            fft_lo_deconvolved = fft_lo / kfft
        else:
            fft_lo_deconvolved = fft_lo 


        below_beamscale = kfft < min_beam_fraction
        #mask = (angscales >= SAS) & (angscales <= LAS) & (~below_beamscale)
        fft_lo_deconvolved[below_beamscale] = np.nan
        return (angscales.to(u.arcsec), np.abs(fft_lo_deconvolved))

    def magnitude_spectrum(self, sigma=0.):
        return 20 * np.log(np.abs(self.to_fourier_space(sigma)))

    def draw_magnitude_spectrum(self, img=None, sigma=0., cmapstr='afmhot', alpha=1., barcolor=False, exp=0.3):
        data = 20 * np.log(np.abs(self.to_fourier_space(sigma)))
        if img is None:
            fig, img = plt.subplots(1,1)#, figsize=(15,15))
        x_range = np.linspace(0,data.shape[0] * self.pixelarea**0.5, data.shape[0]).to(u.arcsec)
        x_range -= np.max(x_range) / 2.
        y_range = np.linspace(0,data.shape[1] * self.pixelarea**0.5, data.shape[1]).to(u.arcsec)
        y_range -= np.max(y_range) / 2.
        extent = (np.min(x_range).value, np.max(x_range).value, np.min(y_range).value, np.max(y_range).value)
        norm = PowerNorm(gamma=exp, clip=False)
        pos_net_clipped = img.imshow(data , cmap=cmapstr, alpha=alpha, norm=norm, extent=extent)
        if barcolor:
            cbar = plt.colorbar(pos_net_clipped, ax=img, extend='both')
            cbar.minorticks_on()
        img.set_title('2D FFT Power Spectrum')
        return (img)

    def to_image_space(self, freq_data, mask=None, filter_pass=None):
        if filter_pass is not None:
            freq_data *= freq_data * filter_pass
        if mask is None:
            f_ishift = np.fft.ifftshift(freq_data)
        else:
            freq_data[mask] = np.min(freq_data)
            f_ishift = np.fft.ifftshift(freq_data)
        img_back = np.fft.ifft2(f_ishift)
        img_back = np.real(img_back)
        # Ver los residuos de obtener la imagen con la parte imaginaria
        return img_back

    def extract_coords(self,x_min, x_max=None, y_min=None, y_max=None):
        if type(x_min) == tuple:
            (x_min, x_max, y_min, y_max) = x_min

        (idx_min, idx_max, idy_min, idy_max) = (np.argmin(np.abs(self.axis[0] - x_min)),np.argmin(np.abs(self.axis[0] - x_max)),np.argmin(np.abs(self.axis[1] - y_min)),np.argmin(np.abs(self.axis[1] - y_max)))
        area = self.data[min(idy_min, idy_max):max(idy_min,idy_max),min(idx_min, idx_max):max(idx_min,idx_max)]
        x0 = np.min(np.array([x_min,x_max]) * sign(self.step[0].value)) * sign(self.step[0].value)
        y0 = np.min(np.array([y_min,y_max]) * sign(self.step[1].value)) * sign(self.step[1].value)
        new_img = Image(data=area, step = self.step, x0=x0, y0=y0)
        return new_img

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

    def fix_header(self, header):
        keys_to_remove = [  'PC1_1', 'PC2_1', 'PC3_1', 'PC4_1', 'PC1_2', 'PC2_2',
                            'PC3_2', 'PC4_2', 'PC1_3', 'PC2_3', 'PC3_3', 'PC4_3',
                            'PC1_4', 'PC2_4', 'PC3_4', 'PC4_4', 
                            'CTYPE3', 'CRVAL3', 
                            'CDELT3', 'CRPIX3', 'CUNIT3', 'NAXIS3',
                            'CTYPE4', 'CRVAL4', 'CDELT4', 'CRPIX4', 'CUNIT4', 'NAXIS4',
                            'PV2_1', 'PV2_2', 'SPECSYS',
                            'ALTRVAL', 'ALTRPIX', 'VELREF', 'HISTORY', 'COMMENT']
        hdu = pyfits.PrimaryHDU()
        if 'BUNIT' not in header.keys():
            hdu.header['BUNIT'] = 'Jy/beam'
        hdu.header['AUTHOR'] ='IMGCombine by d.diaz@irya.unam.mx'
        hdu.header['SIMPLE'] = True
        hdu.header['BITPIX'] = -32
        hdu.header['NAXIS'] = 2
        hdu.header['EXTEND'] = True
        for key in header:
            if key not in keys_to_remove and key not in hdu.header:
                try:
                    hdu.header[key] = header[key]
                except Exception as e:
                    mylog (key, show_exc(e))

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
                hdu.header['CTYPE2'] = 'DEC---SIN'
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

    def tofits(self, filename=None):
        hdu = pyfits.PrimaryHDU()
        hdu.header = self.header
        hdu.data = self.data
        hdu.header['ORIGIN'] = 'IMGCombine'
        hdu.header['OBJECT'] = self.name

        version = 0
        if filename is None:
            filename = 'v%d_%s' % (version,os.path.basename(self.path))
            while os.path.exists(filename):
                version +=1
                filename = 'v%d_%s' % (version,os.path.basename(self.path))
        hdu.writeto(filename, overwrite=True)
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
        if not self.interf and False:
            coords = SkyCoord(l=x, b=y, frame=Galactic)
        else:
            coords = SkyCoord(ra=x, dec=y, frame=FK5)

        return (coords, data[row,col], [row,col])

    def plot(self, coord, img_frame = None, marker='x', sigma=3.):
        if not img_frame:
            fig, img_frame= plt.subplots(1,1)
            plot = True
        else:
            plot = False

        self.draw(img_frame, sigma)
        try:
            img_frame.scatter(coord.ra.value, coord.dec.value, marker=marker)
        except:
            img_frame.scatter(coord.l.value, coord.b.value, marker=marker)
        if plot:
            plt.show()

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


    def janskys(self, freq=None, alpha_spec=3.5):
        if freq is not None:
            freq_factor = (freq.value / self.restfreq().value)**(alpha_spec) 
        else:
            freq_factor = 1.

        beam = self.get_beam()
        npts = self.data.size - np.isnan(self.data).sum()
        mean = self.mean()
        omega_beam = beam.sr.to(u.deg**2)
        omega_pix = (self.header['CDELT1'] * u.Unit(self.header['CUNIT1'])) ** 2
        full_area = omega_pix * npts
        n_beams = full_area / omega_beam
        return mean*n_beams * u.Jy * freq_factor

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

    def mad(self, iters=100, limit=0.99, stack=False, cutoff=2., mask=None):
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
            while (next_mad / mad < limit and counter < iters) or (counter == 0):
                mad = next_mad
                next_mad = stats.mad_std(data[np.where(data < cutoff * mad)], axis=None, ignore_nan=True)
                counter += 1
                mad_stack.append(mad)
                if next_mad <= 0:
                    break

            if stack: 
                return(mad_stack)# * (u.Jy/beam)
            else:
                return(mad_stack[-1])# * (u.Jy/beam)
        except Exception as e:
            mylog(show_exc(e), verbose=True)

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

    def omega_beam (self):
        beam = self.get_beam()
        omega_beam = beam.sr.to(u.deg**2)
        return omega_beam

    def stats(self, title=None, iters=100, mask=None, fmt="standar", headers=True):
        try:
            if not title:
                title = self.name

            omega_pix = self.omega_pix()
            mad = self.mad(iters)
            if fmt == "latex":
                if headers:
                    mylog ("Region       & ", end=" ")
#                     mylog ("Sum          & ", end=" ")
                    mylog ("FluxDensity  [%s] & " % self.janskys().unit, end=" ")
                    mylog ("Jy / Beam    & ", end=" ")
                    mylog ("Mad          & ", end=" ")
                    mylog ("Max          & ", end=" ")
                    mylog ("Min          & ", end=" ")

                    mylog ("Lores Scale  \\\\")
                    #mylog ("Hires Scale  ")

                mylog ("%s & " % (title), end = ' ')
                mylog (" %.4lf &" % (self.janskys().value), end=' ')
                mylog (" %.4lf &" % (self.janskys_beam().value), end=' ')
                mylog (" %.4lf &" % mad, end=' ')
                mylog (" %.4lf &" % self.max(), end=' ')
                mylog (" %.4lf &" % self.min(), end=' ')
#                 mylog (" %.4lf, %.4lf dB &" % (self.dr(iters),20*np.log10(self.dr(iters))), end=' ')
#                 mylog (" %.4lf, %.4lf dB &" % (self.dr_neg(iters),20*np.log10(self.dr_neg(iters))), end=' ')
#                 mylog (" %.8lf &" % (self.poso_neg(mad=mad)), end=' ')
                mylog (" %.4lf \\\\" % (self.lo_scale()))
                #mylog (" %.8lf " % (self.hi_scale()))
            else:
                mylog ("%s " % (title))
                mylog ("\tSum         : %s | %s" % (self.janskys_sum(), self.janskys_pixel()))

                mylog ("\tFluxDensity : %s" % self.janskys())
                mylog ("\tJy / Beam   : %s " % (self.janskys_beam()))
                mylog ("\tField       : %s" % (self.field()))
                mylog (("\t%s" % self.get_beam()).replace('Beam:','Beam        :'))
                mylog ("\tBeam Area   : %s" % (self.omega_beam().to(u.arcsec**2)))
                mylog ("\tPixel Area  : %s" % (self.omega_pix().to(u.arcsec**2)))
                mylog ("\tPixel Size  : %s x %s" % ((abs(self.header['CDELT1']) * u.Unit(self.header['CUNIT1'])).to(u.arcsec),abs(self.header['CDELT2']) * u.Unit(self.header['CUNIT2']).to(u.arcsec)))
                mylog ("\tNumPixels   : %s" % self.Npts())
                mylog ("\tNumBeams    : %s" % (self.Nbeams()))
                mylog ("\tMad         : %.4lf" % mad)
                mylog ("\tMax         : %.4lf" % self.max())
                mylog ("\tMin         : %.4lf" % self.min())
                mylog ("\tDR          : %.4lf, %.4lf dB" % (self.dr(iters),20*np.log10(self.dr(iters))))
                mylog ("\tDR-Minus    : %.4lf, %.4lf dB" % (self.dr_neg(iters),20*np.log10(self.dr_neg(iters))))
                mylog ("\tPoso Neg    : %.8lf" % (self.poso_neg(mad=mad)))
                mylog ("\tLores Scale : %.8lf" % (self.lo_scale()))
                mylog ("\tHires Scale : %.8lf" % (self.hi_scale()))
            if mask:
                self.data = original
        except Exception as e:
            mylog(show_exc(e))

    def mask_pb(self, primary_beam_image, limit=0.2):
        mask = ((primary_beam_image.data < limit) | (np.isnan(primary_beam_image.data)))
        try:
            self.data[mask] = np.nan
        except Exception as e:
            mylog ("Warning in combine.mask_pb!!!! We can not masked the image")

    def angular_scales(self, png_path = None, plot=True):
        lowresfwhm = np.sqrt(self.get_beam().major.to(u.arcsec) * self.get_beam().minor.to(u.arcsec))
        mgps = self

        im_low = mgps.data

        # If weights are given, they must match the shape of the hires data
        weights = 1.

        (nax2,nax1) = mgps.data.shape
        pixscale = (self.pixelarea**0.5).to(u.deg).value

        kfft, ikfft = feather_kernel(nax2, nax1, lowresfwhm, pixscale)
        kfft = np.fft.fftshift(kfft)
        ikfft = np.fft.fftshift(ikfft)

        yy,xx = np.indices([nax2, nax1])
        rr = ((xx-(nax1-1)/2.)**2 + (yy-(nax2-1)/2.)**2)**0.5
        angscales = nax1/rr * pixscale*u.deg

        fft_lo = np.fft.fftshift(np.fft.fft2(np.nan_to_num(im_low * weights)))
        fft_lo_deconvolved = fft_lo / kfft

        below_beamscale = kfft < 0.1
        below_beamscale_plotting = kfft < 1e-3
        fft_lo_deconvolved[below_beamscale_plotting] = np.nan

        mask = (angscales >= lowresfwhm)  & (~below_beamscale)
        assert mask.sum() > 0

        if plot:
            import pylab as pl
            pl.clf()
            pl.title(self.name)
            pl.subplot(1,1,1)
            srt = np.argsort(angscales.to(u.arcsec).value[~below_beamscale_plotting])
            pl.loglog(angscales[mask].to(u.arcsec).value, np.abs(fft_lo_deconvolved[mask]), 'b.', alpha=1, label='Lo-res')
            ylim = pl.gca().get_ylim()
            pl.grid()
            if not png_path:
                plt.show()
            else:
                plt.savefig(os.path.join(png_path, '{0}_compare.png'.format(int(prefix))))

        return (angscales[mask])

    def __unicode__(self):
        return ("%s %s" % (self.path, self.name))

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

    images = []
    for idx, fpath in enumerate(fits):
        images.append(Image(fpath))

    if args.image:
        number_of_subplots=len(fits)
        
        if number_of_subplots > 1:
            ncols = 2
            nrows = math.ceil(number_of_subplots  / ncols)
            figsize = (10*ncols,10*nrows)
            wcs = WCS(images[0].header)
            for idx,img in enumerate(images):
                ax = plt.subplot(nrows, ncols, idx+1, projection=wcs)
                img.draw(ax, exp=args.exp, cmapstr='afmhot')
#             for i in range(number_of_subplots, len(axs)):
#                 fig.delaxes(axs[i])
            plt.show()
        else:
            img=images[0]
            wcs = WCS(img.header)
            img.draw(plt.subplot(1,1,1, projection=wcs), exp=args.exp, cmapstr='afmhot')
            plt.show()
    else:
        for img in images:
            if args.headers:
                print (f'FILE: {fpath}')
                img.show_headers(onlykeys=False)
            if args.fields:
                img.field_show()

        if len(images) > 1:
            combinations = []
            for idx,img in enumerate(images[:-1]):
                for idy,img2 in enumerate(images[idx+1:]):
                    combinations.append(Combination(img, img2))

