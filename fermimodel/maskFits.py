#!/usr/bin/env python
import argparse
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
cfloat_min = 1.e-37

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

def MaskFits(fitsfile, out='maskedimage.fits', img_hdu=None, mask_type=None, radius=180., radius2=None, angle=0., center=(0., 0.), extent=[180., -180., 90., -90.], frame='galactic', unit='degree', clobber=False):
    """Mask the fits image"""
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
        energies_shape = hdu_list['ENERGIES'].data.energy.shape

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
        m_energies[:] = hdu_list['ENERGIES'].data.energy[:]
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
        # hdu_list[img_hdu].data = np.where(hdu_list[img_hdu].data > 0, hdu_list[img_hdu].data, cfloat_min*np.ones(hdu_list[img_hdu].data.shape))
        hdu_list[img_hdu].data[hdu_list[img_hdu].data < cfloat_min] = cfloat_min
    except MemoryError:
        for idx in range(naxis_E):
            hdu_list[img_hdu].data[idx,:,:][hdu_list[img_hdu].data[idx,:,:] < cfloat_min] = cfloat_min

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

    # Debug
    # import matplotlib.pyplot as plt
    # fig = plt.figure()
    # plt.imshow(mask, origin='lower')
    # plt.show()

def cli():
    """Command line interface"""
    helpString = "Mask fits image."
    parser = argparse.ArgumentParser(description=helpString)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('input', type=str, help="Fits file containing the image to which the program will apply a mask")
    parser.add_argument('-o', '--output', type=str, default='masked.fits', help="Optional filename to save the masked data. Default is masked.fits")
    parser.add_argument('-m', '--mask', type=str, choices=['radial', 'square'], default='radial', help='Shape of the mask.')
    parser.add_argument('-ih', '--image_hdu', help='HDU containing the image to be masked. If none is input it is assumed to live in the PRIMARY HDU.')
    parser.add_argument('-cl', '--clobber', action='store_true', help='Flag to overwrite a file of the same name. Default is false.')

    frame_group = parser.add_mutually_exclusive_group()
    frame_group.add_argument('-fk5', '--J2000', action='store_const', const='fk5', help='Flag sets coordinates of mask center to RAJ2000, DECJ2000')
    frame_group.add_argument('-icrs', '--celestial', action='store_const', const='icrs', help='Flag sets coordinates of mask center to RA, DEC')
    frame_group.add_argument('-gal', '--galactic', action='store_const', const='galactic', help='Flag sets coordinates of mask center to GLON, GLAT')
    frame_group.add_argument('-pix', '--pixel', action='store_const', const='pixel', help='Flag sets coordinates of mask center to PIXEL1, PIXEL2')
    parser.add_argument('-u', '--unit', type=str, default='degree', help='Units of mask coordinates. Default is degrees.')

    radial = parser.add_argument_group('radial', 'Parameters for radial mask')
    square = parser.add_argument_group('square', 'Parameters for square mask')

    radial.add_argument('-r', '--radius', type=float, default=180., help='Radius of the mask.')
    radial.add_argument('-r2', '--radius2', type=float, help='Second radius of ellipse. Not yet implemented')
    radial.add_argument('-a', '--angle', type=float, help='Angle of radius with respect to the horizontal axis. not yet implemented.')
    radial.add_argument('-xc', '--horizontal-center', type=float, help='Horizontal coordinate of mask center.')
    radial.add_argument('-yc', '--vertical-center', type=float, help='Vertical coordinate of mask center.')

    square.add_argument('-xmin', '--horizontal-min', type=float, help='Minimum horizontal coordinate value.')
    square.add_argument('-xmax', '--horizontal-max', type=float, help='Maximum horizontal coorindate value.')
    square.add_argument('-ymin', '--vertical-min', type=float, help='Minimum vertical coordinate value.')
    square.add_argument('-ymax', '--vertical-max', type=float, help='Maximum vertical coordinate value.')

    args = parser.parse_args()

    if args.J2000 is not None:
        frame = args.J2000
    elif args.celestial is not None:
        frame = args.celestial
    elif args.galactic is not None:
        frame = args.celestial
    elif args.pixel is not None:
        frame = args.pixel
    else:
        print "cannot mask image with frame type {0}".format(args.frame)
        exit()

    out = MaskFits(args.input, out=args.output, img_hdu=args.image_hdu, mask_type=args.mask, radius=args.radius, radius2=args.radius2, angle=args.angle, center=args.center, extent=args.extent, frame=frame, unit=args.unit, clobber=args.clobber)

    if isinstance(out, tuple):
        print "Output file saved at {0}".format(out[0])
        print "Calculated flux is {0} photons/cm^2/s".format(out[1])
    else:
        print "Output file saved at {0}".format(out)

    # # Debug
    # import matplotlib.pyplot as plt
    # out, flux = MaskFits(args.input, out=args.output, img_hdu=args.image_hdu, mask_type=args.mask, radius=args.radius, radius2=args.radius2, angle=args.angle, center=(args.horizontal_center, args.vertical_center), extent=[args.horizontal_min, args.horizontal_max, args.vertical_min, args.vertical_max], clobber=args.clobber)
    # print flux
    # gd = pyfits.open(out)
    # wcs = WCS(gd[0].header, naxis=2)
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection=wcs)
    # plt.imshow(gd[0].data[0])
    # plt.show()

if __name__ == '__main__': cli()