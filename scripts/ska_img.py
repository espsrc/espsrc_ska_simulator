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

from astropy.io import fits as pyfits
from astropy.io import fits
from astropy.wcs import WCS
from astropy import stats, units as u
from astropy.coordinates import SkyCoord, Angle, Galactic, ICRS, FK5, FK4

from scipy.ndimage.measurements import center_of_mass

from radio_beam import Beam
import reproject as rp

from utils import show_exc

VERBOSE = True


def _theta_image_from_pa_sky(center_sc: SkyCoord, pa_sky_deg: float, w: WCS) -> float:
    """
    Devuelve theta (rad) en el plano de imagen (x->derecha, y->ARRIBA),
    a partir del PA astronómico (grados E de N). Usa un paso angular pequeño
    y WCS para medir el ángulo real en píxeles.
    """

    if hasattr(pa_sky_deg, 'unit'):
        pa_sky_deg = pa_sky_deg.to(u.deg).value

    x0, y0 = w.world_to_pixel(center_sc)

    step = 1.0 * u.arcsec  
    p1 = center_sc.directional_offset_by(pa_sky_deg * u.deg, step)


    x1, y1 = w.world_to_pixel(p1)


    dx = x1 - x0
    dy_up = -(y1 - y0)  
    theta = np.arctan2(dy_up, dx) 
    return theta

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

    def __init__(self, path = None, data=None, step = None, x0 = 0, y0 = 0, unit=u.deg, interf=False, frame=None, header=None, extent=None, bunit=None):
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
        elif 'RESTFREQ' in self.header.keys():
            return (self.header['RESTFREQ'] * u.Hz).to(unit)
        elif 'CRVAL3' in self.header.keys():
            return (self.header['CRVAL3'] * u.Hz).to(unit)
        else:
            return None

    def pix2coords(self,x,y=None):
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
                wcs = self.wcs2d
                if wcs.naxis < 2:
                    raise Exception ("WCS has less than 2 axis")
                if wcs.naxis > 2:
                    x0, y0, _, _  = wcs.wcs_pix2world(x, y, 0, 0, 0)
                else:
                    x0, y0 = wcs.wcs_pix2world(x, y, 0)
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
            wcs = self.wcs2d
            x0, y0 = wcs.wcs_world2pix(x, y, 0)
        except:
            wcs = self.wcs2d
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

    def get_beam(self):
        try:
            beam = Beam.from_fits_header(self.header)
            return (beam)
        except Exception as e:
            mylog(show_exc(e), verbose=True, flush=True)
            return None
        return None

    def zoom(self, zoom):
        ''' Zoom into a region of the image. 
            zoom = [ra_center, dec_center, radius]
        '''
        orig_data = np.copy(self.data)
        (ra_center, dec_center, radius_zoom) = zoom
        zoom_coords = SkyCoord(ra_center, dec_center, frame=self.frame, unit=(u.hourangle, u.deg))
        if self.frame != zoom_coords.frame:
            zoom_coords = zoom_coords.transform_to(self.frame)
        ra_center = zoom_coords.ra.to(u.deg).value
        dec_center = zoom_coords.dec.to(u.deg).value
        radius_zoom = Angle(radius_zoom).to(u.deg).value
        # Extract pixels corresponding to the zoom limits in pixels coordinates
        # x_center, y_center = self.coords2pix(zoom_coords)
        low_left = SkyCoord(ra_center + radius_zoom, dec_center - radius_zoom, frame=self.frame, unit=u.deg)
        up_right = SkyCoord(ra_center - radius_zoom, dec_center + radius_zoom, frame=self.frame, unit=u.deg)
        (x0, y0) = self.coords2pix(low_left)
        (x1, y1) = self.coords2pix(up_right)

        #Create wcs with the zoomed region
        zoomHeader = self.header
        zoomHeader['CRPIX1'] = 1
        zoomHeader['CRPIX2'] = 1
        zoomHeader['CRVAL1'] = low_left.ra.to(u.deg).value
        zoomHeader['CRVAL2'] = low_left.dec.to(u.deg).value
        zoomHeader['NAXIS1'] = int(abs((x1 - x0)))
        zoomHeader['NAXIS2'] = int(abs((y1 - y0)))
        subimg = orig_data[int(min(y0,y1)):int(max(y0,y1)), int(min(x0,x1)):int(max(x0,x1))]

        zoomImg = SKAImage(data=subimg, header=zoomHeader)
        return zoomImg

    
    def draw(self, img=None, exp=0.3, cmapstr='gray', scale='power', colorbar=False, plot=False, title=None, show_beam = True):
        
        from astropy.wcs import WCS
        if self.data is not None:
         
            # Plot using matplotlib with WCS projection
            if img is None:
                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(1, 1, 1, projection=self.wcs2d)
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
                im = ax.imshow(self.data2d, origin='lower', norm=norm, cmap=cmapstr)
                cbar = plt.colorbar(im, ax=ax, extend='both', format=FormatStrFormatter('%.2e'))
                cbar.set_label(f"[{self.header['BUNIT']}]")
            else:
                im = ax.imshow(self.data2d, origin='lower', norm=norm, cmap=cmapstr)

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
            return ax
        return None


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
            beam_pa_value = np.sign(self.header['CDELT1']) * beam.pa.value
            if beam.minor.value != beam.major.value:
                try:
                    frame = patches.Ellipse((cx,cy), width=beam.minor.value, height=beam.major.value, color=color, angle=beam_pa_value, transform=img.get_transform('world'))
                except:
                    frame = patches.Ellipse((cx,cy), width=beam.minor.value, height=beam.major.value, color=color, angle=beam_pa_value)
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
        hdu.header['NAXIS'] = 2
        hdu.header['EXTEND'] = header['EXTEND'] if 'EXTEND' in header.keys() else True
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
        hdu.header['ORIGIN'] = 'SKA-Image'
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


    def janskys(self, freq=None, alpha_spec=3.5, center = None, a=None, b=None, pa=None):
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

        data = self.data2d.copy()
        if center is not None:
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
            data[~(mask)] = np.nan
        

        beam = self.get_beam()
        npts = data.size - np.isnan(data).sum()
        mean = np.nanmean(data)
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

 

    def mask_pb(self, primary_beam_image, limit=0.2):
        mask = ((primary_beam_image.data < limit) | (np.isnan(primary_beam_image.data)))
        try:
            self.data[mask] = np.nan
        except Exception as e:
            mylog ("Warning in combine.mask_pb!!!! We can not masked the image")

  

    def __unicode__(self):
        return ("%s %s" % (self.path, self.name))
    
    def draw_source(self, ax, coord, a, b = None, pa = None, color="green", ellipse=False, alpha=0.3, label = None, verbose=False):
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

            if ellipse:

                beam_pa_value = 180 -pa.value
                # beam_pa_value = (_theta_image_from_pa_sky(coord, pa.value, self.wcs2d) * u.rad).to(u.deg).value
                minor = min(a,b)
                major = max(a,b)
                if minor != major:
                    try:
                        frame = patches.Ellipse((cx,cy), width=minor.to(u.deg).value, height=major.to(u.deg).value, color=color, angle=beam_pa_value, transform=ax.get_transform('world'), fill=False, alpha=alpha)
                    except Exception as e:
                        print(show_exc(e))
                        frame = patches.Ellipse((cx,cy), width=minor.to(u.deg).value, height=major.to(u.deg).value, color=color, angle=beam_pa_value, fill=False, alpha=alpha)
                else:
                    try:
                        frame = patches.Circle((cx,cy), major.to(u.deg).value, color=color, angle=pa.value, transform=ax.get_transform('world'), fill=False, alpha=alpha)
                    except:
                        frame = patches.Circle((cx,cy), minor.value, color=color, angle=pa.value, fill=False, alpha=alpha)
                ax.add_patch(frame)


        except Exception as e:
            mylog (show_exc(e))
    
    @property
    def wcs(self):
        try:
            wcs = WCS(self.header)
            return wcs
        except:
            return None

    @property
    def wcs2d(self):
        try:
            wcs = WCS(self.header)
            if wcs.naxis >= 2:
                return WCS(self.fix_header(self.header))
            else:
                return None
        except:
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
