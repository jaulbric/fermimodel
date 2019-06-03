#!/usr/bin/env python
import numpy as np
import astropy.io.fits as pyfits
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame
from astropy.coordinates import SkyCoord
import warnings
import os
import datetime
from scipy.integrate import trapz
from tempfile import mkdtemp

d2r = np.pi/180.

class MaskTypeError(Exception):
    """Raise if the type of mask requested is not allowed"""
    pass

class NaxisError(Exception):
    """Raise if the energy axis of the image cannot be found."""
    pass

def angsep(xref, yref, center=(0.,0.)):
    """Calculates angular separation between two points on the sky. Parameters and returns are in degrees."""
    xref = xref*d2r
    yref = yref*d2r
    center = (center[0]*d2r, center[1]*d2r)
    return np.arccos((np.cos(yref)*np.cos(center[1])*np.cos(xref - center[0])) + (np.sin(yref)*np.sin(center[1])))/d2r

def pixsep(xref, yref, center=(0., 0.)):
    """Calculates pixel separation between two points in the image. Parameters and returns are in pixels."""
    return np.sqrt((xref - center[0])**2 + (yref - center[1])**2)

def genRadialMask(xcoords, ycoords, radius, radius2, angle, center, frame):
    """Generate a radial mask"""
    if frame != 'pixel':
        dist = angsep(xcoords, ycoords, center)
    else:
        dist = pixsep(xcoords, ycoords, center)
    
    if radius2 is None:
        radius2 = radius

    mask = np.where(dist <= radius, np.ones(xcoords.shape), np.zeros(xcoords.shape))
    return mask

def genSquareMask(xcoords, ycoords, extent):
    """Generate a square mask"""
    mask = np.where((xcoords > extent[0]) & (xcoords < extent[1]) & (ycoords > extent[2]) & (ycoords < extent[3]), np.ones(xcoords.shape), np.zeros(xcoords.shape))
    return mask

def integrateMapCube(data, xcoords, ycoords, energies):
    """Integrate the entire map cube and return the total flux in #/cm^2/s"""
    dflux = np.array([trapz(trapz(data[idx,:,:], np.cos(ycoords*d2r)*xcoords*d2r, axis=1), ycoords[:,0]*d2r, axis=0) for idx in range(len(energies))])
    flux = trapz(dflux, energies)
    return flux

def maskFits(fitsfile, out='maskedimage.fits', img_hdu=None, mask_type=None, radius=180., radius2=None, angle=0., center=(0., 0.), extent=[-180., 180., -90., 90.], frame='galactic', unit='degree', clobber=False, float_min=1.17549e-38):
    """Mask the fits image

    Parameters
    ----------
    fitsfile : str
        Path to the FITS file containing the the image which will be masked.
    out : str (Optional)
        Name of the output FITS file for which the mask wiil be applied
    img_hud : int or float (Optional)
        Name or integer for the FITS hdu containing the image. Default is 'Primary'.
    mask_type : str
        The geometry of the mask to be applied. Choices are 'radial' or 'square'.
    radius : float (Optional)
        Radius of the mask if mask_type is radial. Default is 180.
    radius2 : float (Optional)
        Second radius of the mask if the mask is not symmetric. Default is to use a symmetric mask.
    angle : float (Optional)
        Rotation angle of the ellipse. Default is 0.
    center : tuple (Optional)
        Center coordinates (C1, C2) of the radial mask. Default is (0., 0.).
    extent : list (Optional)
        [xmin, xmax, ymin, ymax] extent of the square mask. Default is [-180., 180., -90., 90.].
    frame : str (Optional)
        Coordinate frame to use with mask coordinates. Choices are 'galactic', 'icrs', 'fk5', 'pixel'. Default is 'galactic'.
    unit : str (Optional)
        Units of coordinates. Default is 'degree'.
    clobber : bool
        Flag to overwrite a file of the same name if it exists. Default is False.
    float_min : float
        Minimum float value to use for pixels in the image after masking. gtobssim doesn't like pixel values <= 0. Default is 1.17549e-38.

    Returns
    -------
    out : str or tuple
        If the input fits image is 2D the output is the full path to the masked fits image. If the input fits image is 3D the output is a tuple whose first entry is the full path to the masked fits image and whose second entry is the integrated flux of the masked fits image.

    """
    if frame == 'galactic':
        frame_str = '(GLAT, GLON)'
    elif frame == 'icrs':
        frame_str = '(RA, DEC)'
    elif frame == 'fk5':
        frame_str = '(RAJ2000, DECJ2000)'
    elif frame == 'pixel':
        frame_str = '(PIX1, PIX2)'
    else:
        raise IOError("Invalid frame {0}".format(frame))
    
    fitsfile = os.path.expandvars(fitsfile).replace("$(FERMI_DIR)", os.environ.get("FERMI_DIR")) if os.environ.get("FERMI_DIR") is not None else os.path.expandvars(fitsfile)
    hdu_list = pyfits.open(fitsfile)
    
    if img_hdu is None:
        # Assume the data image is in the primary hdu
        header = hdu_list['PRIMARY'].header
        img_hdu = 'PRIMARY'
    else:
        try:
            img_hdu = int(img_hdu)
        except ValueError:
            pass
        header = hdu_list[img_hdu].header

    wcs = WCS(header, naxis=2)
    # Probably shouldn't use these private attributes, but there doesn't seem to be another way to get the length of each axis from the WCS object.
    naxis1 = wcs._naxis1
    naxis2 = wcs._naxis2

    xgrid, ygrid = np.meshgrid(np.arange(naxis1), np.arange(naxis2))

    if frame != 'pixel':
        flat_xcoords, flat_ycoords = wcs.wcs_pix2world(xgrid.flatten(), ygrid.flatten(), 0)
        xcoords = flat_xcoords.reshape(naxis2, naxis1)
        ycoords = flat_ycoords.reshape(naxis2, naxis1)
    else:
        xcoords = xgrid
        xcoords = ygrid

    del xgrid
    del ygrid

    if mask_type == 'radial':
        if frame != 'pixel':
            c = SkyCoord(center[0], center[1], frame=frame, unit=unit)
            c_new = c.transform_to(wcs_to_celestial_frame(wcs))
            center_new = map(float, c_new.to_string().split(' '))
        else:
            center_new = center

        mask = genRadialMask(xcoords, ycoords, radius, radius2, angle, center_new, frame=frame)
        hdu_list[img_hdu].header['history'] = '{0} Applied radial mask to data.'.format(datetime.datetime.today().strftime('%d %B %Y'))
        hdu_list[img_hdu].header['history'] = 'radius={0}, radius2={1}, angle={2}, center={3} {4}'.format(radius, radius2, angle, center, frame_str)
    elif mask_type == 'square':
        if frame != 'pixel':
            c = SkyCoord([extent[0], extent[2]], [extent[1], extent[3]], frame=frame, unit=unit)
            c_new = c.transform_to(wcs_to_celestial_frame(wcs))
            l, t = map(float, c_new[0].to_string().split(' '))
            r, b = map(float, c_new[1].to_string().split(' '))
            extent_new = [l, r, b, t]
        else:
            extent_new = extent
        if 'GLON' in wcs.wcs.ctype[0]:
            mask = genSquareMask(np.where(xcoords < 180., xcoords, xcoords - 360.), ycoords, extent_new)
        else:
            mask = genSquarMask(xcoords, ycoords, extent_new)
        hdu_list[img_hdu].header['history'] = '{0} Applied square mask to data.'.format(datetime.datetime.today().strftime('%d %B %Y'))
        hdu_list[img_hdu].header['history'] = '[left, right, top, bottom]={0} {1}'.format(extent, frame_str)
    else:
        raise MaskTypeError("{0} is not a supported mask geometry.".format(mask_type))

    data_shape = hdu_list[img_hdu].data.shape
    xcoords_shape = xcoords.shape
    ycoords_shape = ycoords.shape

    naxis = header['NAXIS']
    if naxis > 2:
        tile_shape = [1,1,1]
        E_idx = 0
        for idx in range(1,4):
            naxis_E = None
            if header['CTYPE{0}'.format(idx)] not in wcs.wcs.ctype:
                E_idx = idx
                naxis_E = header['NAXIS{0}'.format(idx)]
                tile_shape[naxis - idx] = naxis_E
        if naxis_E is None:
            raise NaxisError("Could not find the energy axis of the image.") 
        
        try:
            # hdu_list[img_hdu].data = hdu_list[img_hdu].data*np.tile(mask, tuple(tile_shape))
            hdu_list[img_hdu].data *= np.tile(mask, tuple(tile_shape))
        except MemoryError:
            for idx in range(naxis_E):
                hdu_list[img_hdu].data[idx,:,:] *= mask
        del mask

        # Have to memmap data arrays because of intermediate arrays created by trapz
        tmp_dir = mkdtemp()
        try:
            energies_shape = hdu_list['ENERGIES'].data.energy.shape
        except AttributeError:
            energies_shape = hdu_list['ENERGIES'].data.Energy.shape

        m_data_path = os.path.join(tmp_dir, 'data_array.dat')
        m_energies_path = os.path.join(tmp_dir, 'energy_array.dat')
        m_xcoords_path = os.path.join(tmp_dir, 'xcoords_array.dat')
        m_ycoords_path = os.path.join(tmp_dir, 'ycoords_array.dat')

        m_data = np.memmap(m_data_path, dtype=np.float32, mode='w+', shape=data_shape)
        m_energies = np.memmap(m_energies_path, dtype=np.float32, mode='w+', shape=energies_shape)
        m_xcoords = np.memmap(m_xcoords_path, dtype=np.float32, mode='w+', shape=xcoords_shape)
        m_ycoords = np.memmap(m_ycoords_path, dtype=np.float32, mode='w+', shape=ycoords_shape)

        if wcs.wcs.cdelt[0] < 0.:
            m_data[:,:,:] = np.flip(hdu_list[img_hdu].data, axis=2)[:,:,:]
            m_xcoords[:,:] = np.flip(np.where(xcoords <= 180., xcoords, xcoords - 360.), axis=1)[:,:]
        else:
            m_data[:,:,:] = hdu_list[img_hdu].data[:,:,:]
            m_xcoords[:,:] = np.where(xcoords <= 180, xcoords, xcoords - 360.)[:,:]
        m_ycoords[:,:] = ycoords[:,:]
        try:
            m_energies[:] = hdu_list['ENERGIES'].data.energy[:]
        except AttributeError:
            m_energies[:] = hdu_list['ENERGIES'].data.Energy[:]
        del m_data
        del m_xcoords
        del m_ycoords
        del m_energies

    else:
        hdu_list[img_hdu].data *= mask
        del mask

    if not os.path.isabs(out):
        out = os.path.join(os.getcwd(), out)

    try:
        hdu_list[img_hdu].data[hdu_list[img_hdu].data < float_min] = float_min
    except MemoryError:
        for idx in range(naxis_E):
            hdu_list[img_hdu].data[idx,:,:][hdu_list[img_hdu].data[idx,:,:] < float_min] = float_min

    try:
        hdu_list.writeto(out)
    except IOError:
        if clobber:
            os.remove(out)
            hdu_list.writeto(out)
        else:
            raise
    hdu_list.close()

    if naxis > 2:
        m_data = np.memmap(m_data_path, dtype=np.float32, mode='r', shape=data_shape)
        m_xcoords = np.memmap(m_xcoords_path, dtype=np.float32, mode='r', shape=xcoords_shape)
        m_ycoords = np.memmap(m_ycoords_path, dtype=np.float32, mode='r', shape=ycoords_shape)
        m_energies = np.memmap(m_energies_path, dtype=np.float32, mode='r', shape=energies_shape)

        flux = integrateMapCube(m_data, m_xcoords, m_ycoords, m_energies)

        del m_data
        del m_xcoords
        del m_ycoords
        del m_energies
        del xcoords
        del ycoords
        return out, flux
    else:
        return out