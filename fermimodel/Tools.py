import os
import numpy as np
import astropy.io.fits as pyfits
from astropy import units
from scipy.integrate import quad, trapz
from scipy.special import gamma, gammaincc
from xml.dom import minidom

# class GetFluxError(Exception):
#     """Raise when we cannot calculate the flux of the source from parameters"""
#     pass

# class WriteSpectrumError(Exception):
#     """Raised when spectrum type cannot be written to file."""
#     pass

# class HeaderCheckError(Exception):
#     """Raised when diffuse emission is missing essential headers."""
#     pass

# class ExtendedTemplateError(Exception):
#     """Raise when an extended template cannot be found."""
#     pass

def getFlux(spectype, emin, emax, **spectrumargs):
    """Calculate integrated flux between emin and emax using spectral parametes

    For power law spectral type the indefinite integral of the differential flux is

        $\frac{E K \left( \frac{E}{E_{0}} \right)^{- \gamma}}{1 - \gamma}$

    For log parabola spectral type the indefinite integral has no close form and must be numerically integrated

    For power law with subexponential cutoff the indefinite integral of the differential flux is

        $- \frac{K a^{\frac{\gamma - 1}{b}} E^{\gamma}_{0} e^{a E^{b}_{0}} \Gamma \left( \frac{1 - \gamma}{b}, a E^{b} \right)}{b}$

    Where $\Gamma \left( \frac{1 - \gamma}{b}, a E^{b} \right)$ is the (upper) incomplete gamma function.

    Parameters
    ----------
    spectype : str
        The spectral type in the global model (PowerLaw, LogParabola, PLSuperExpCutoff).
    emin : float
        Minimum energy for integrated flux
    emax : float
        Maximum energy for integrated flux
    pl_flux_density : float, optional
        The differential flux at pivot_energy for the PowerLaw fit, in photons/cm2/MeV/s.
    lp_flux_density : float, optional
        The differential flux at pivot_energy for the LogParabola fit, in photons/cm2/MeV/s.
    plec_flux_density : float, optional
        The differential flux at pivot_energy for the PLSuperExpCutoff fit, in photons/cm2/MeV/s.
    pivot_energy : float, optional
        The energy, in MeV, at which the error in the differential photon flux is minimal (i.e., the decorrelation energy for the power-law fit).
    pl_index : float, optional
        The photon index for the PowerLaw fit.
    lp_index : float, optional
        Photon Index at Pivot Energy for LogParabola fit.
    lp_beta : float, optional
        The curvature parameter (beta) for LogParabola spectrum types.
    plec_index : float, optional
        The low energy photon index for PLSuperExpCutoff fit.
    plec_expfactor : float, optional
        The exponential factor for PLSuperExpCutoff fit.
    plec_exp_index : float, optional
        The exponential index for PLSuperExpCutoff fit.

    Returns
    -------
    flux : float
        flux in units of #/m^2/s

    """
    flux = 0.
    if spectype == 'PowerLaw':
        plf= spectrumargs['pl_flux_density']
        p = spectrumargs['pivot_energy']
        pli = spectrumargs['pl_index']
        integrated_dflux = lambda energy: (energy*plf*(energy/p)**(-pli))/(1. - pli)
        flux = integrated_dflux(emax) - integrated_dflux(emin)
    elif spectype == 'LogParabola':
        lpf = spectrumargs['lp_flux_density']
        p = spectrumargs['pivot_energy']
        lpb = spectrumargs['lp_beta']
        lpi = spectrumargs['lp_index']
        dflux = lambda energy: lpf*((energy/p)**(- lpi - lpb*np.log(energy/p)))
        flux, fluxunc = quad(dflux, emin, emax, epsabs=0)
    elif spectype == 'PLSuperExpCutoff2':
        cof = spectrumargs['plec_flux_density']
        plecef = spectrumargs['plec_expfactor']
        pleci = spectrumargs['plec_index']
        p = spectrumargs['pivot_energy']
        plecei = spectrumargs['plec_exp_index']
        integrated_dflux = lambda energy: - (cof*(plecef**((pleci - 1.)/plecei))*(p**pleci)*np.exp(plecef*(p**plecei))*gamma((1.-pleci)/plecei)*gammaincc((1.-pleci)/plecei, plecef*(energy**plecei)))/plecei
        flux = integrated_dflux(emax) - integrated_dflux(emin)
        if np.isnan(flux):
            dflux = lambda energy: cof*((energy/p)**(-pleci))*np.exp(plecef*(p**plecei - energy**plecei))
            flux, fluxunc = quad(dflux, emin, emax, epsabs=0)
    elif spectype == 'BrokenPowerLaw':
        plf = spectrumargs['pl_flux_density']
        p = spectrumargs['pivot_energy']
        pli = spectrumargs['pl_index']
        ebreak = spectrumargs['ebreak']
        gamma2 = spectrumargs['gamma2']
        integrated_dflux1 = lambda energy: (energy*plf*(energy/p)**(-pli))/(1. - pli)
        integrated_dflux2 = lambda energy: (energy*plf*(energy/p)**(-gamma2))/(1. - gamma2)
        flux = integrated_dflux2(emax) - integrated_dflux2(ebreak) + integrated_dflux1(ebreak) - integrated_dflux1(emin)
    else:
        raise GetFluxError("FluxError: Cannot calculate flux for {0} spectral type.".format(spectype))

    if np.isnan(flux):
        raise GetFluxError("FluxError: Calculated flux is nan.")

    return flux*1.e4

def writeSpectrum(name, spectype, emin, emax, directory, float_min=1.e-37, **spectrumargs):
    """Write the spectrum to file.

    Parameters
    ----------
    spectype : str
        The spectral type in the global model (PowerLaw, LogParabola, PLSuperExpCutoff).
    name : str
        Name of the source. The spectrum file will be saved under this name.
    emin : float
        Minimum energy for integrated flux
    emax : float
        Maximum energy for integrated flux
    directory : str
        Directory in which to save spectrum file
    float_min : float
        Minimum value for differential flux.
    pl_flux_density : float
        The differential flux at pivot_energy for the PowerLaw fit, in photons/cm2/MeV/s.
    lp_flux_density : float
        The differential flux at pivot_energy for the LogParabola fit, in photons/cm2/MeV/s.
    plec_flux_density : float
        The differential flux at pivot_energy for the PLSuperExpCutoff fit, in photons/cm2/MeV/s.
    pivot_energy : float
        The energy, in MeV, at which the error in the differential photon flux is minimal (i.e., the decorrelation energy for the power-law fit).
    pl_index : float
        The photon index for the PowerLaw fit.
    lp_index : float
        Photon Index at Pivot Energy for LogParabola fit.
    lp_beta : float
        The curvature parameter (beta) for LogParabola spectrum types.
    plec_index : float
        The low energy photon index for PLSuperExpCutoff fit.
    plec_expfactor : float
        The exponential factor for PLSuperExpCutoff fit.
    plec_exp_index : float
        The exponential index for PLSuperExpCutoff fit.

    Returns
    -------
    filename : str
        Full path to spectrum file

    """
    # Errors in gtobssim, probably from flux normalization, cause gtobssim to crash if the differential flux is below a certain number. Set the minimum float to something C++ should always be able to handle. Comment this out if we ever fix this bug.
    name = name if name[0] != '_' else name[1:]
    filename = os.path.join(directory,name.replace(" ","")) + '.dat'

    with open(filename, "wb") as f: 
        if spectype == 'PowerLaw':
            plf = spectrumargs['pl_flux_density']
            p = spectrumargs['pivot_energy']
            pli = spectrumargs['pl_index']
            num = np.ceil(10.*(np.log10(emax) - np.log10(emin))).astype(int) + 1 # Number of energy bins, roughly 10 per decade
            E = np.logspace(np.log10(emin), np.log10(emax), num=num)
            dflux = plf*((E/p)**(-pli))
        elif spectype == 'LogParabola':
            lpf = spectrumargs['lp_flux_density']
            p = spectrumargs['pivot_energy']
            lpi = spectrumargs['lp_index']
            lpb = spectrumargs['lp_beta']
            num = np.ceil(10.*(np.log10(emax) - np.log10(emin))).astype(int) + 1 # Number of energy bins, roughly 10 per decade
            E = np.logspace(np.log10(emin), np.log10(emax), num=num)
            dflux = lpf*((E/p)**(- lpi - lpb*np.log(E/p)))
        elif spectype == 'PLSuperExpCutoff2':
            cof = spectrumargs['plec_flux_density']
            p = spectrumargs['pivot_energy']
            pleci = spectrumargs['plec_index']
            plecef = spectrumargs['plec_expfactor']
            plecei = spectrumargs['plec_exp_index']
            num = np.ceil(10.*(np.log10(emax) - np.log10(emin))).astype(int) + 1 # Number of energy bins, roughly 10 per decade
            E = np.logspace(np.log10(emin), np.log10(emax), num=num)
            dflux = cof*((E/p)**(-pleci))*np.exp(plecef*((p**plecei) - (E**plecei)))
        elif spectype == 'BrokenPowerLaw':
            plf = spectrumargs['pl_flux_density']
            p = spectrumargs['pivot_energy']
            pli = spectrumargs['pl_index']
            gamma2 = spectrumargs['gamma2']
            ebreak = spectrumargs['ebreak']
            num1 = np.ceil(10.*(np.log10(ebreak) - np.log10(emin))).astype(int) + 1 # Number of energy bins, roughly 10 per decade
            num2 = np.ceil(10.*(np.log10(emax) - np.log10(ebreak))).astype(int) + 1 # Number of energy bins, rooughly 10 per decade
            E1 = np.logspace(np.log10(emin), np.log10(ebreak), num=num1, endpoint=False)
            E2 = np.logspace(np.log10(ebreak), np.log10(emax), num=num2)
            dflux1 = plf*((E1/p)**(- pli))
            dflux2 = plf*((E2/p)**(- gamma2))
            dflux = np.append(dflux1, dflux2)
            E = np.append(E1, E2)
        else:
            raise WriteSpectrumError("WriteSpectrumError: {0} is not a spectral type that can be written to file.".format(spectype))

        dflux[dflux < float_min] = float_min
        lines = ['{0} {1}'.format(en, df) for en, df in zip(E, dflux)]
        
        f.write('\n'.join(lines))

    return filename

def checkHeader(fitsfile):
    """Checks the galactic diffuse model for the correct headers.

    As of 04/09/2019 the standard diffuse emission file gll_iem_v07.fits is missing some headers that causes gtobssim to crash. This function will raise an exception if the headers are missing.

    Parameters
    ----------
    fits : str
        Path to galactic diffuse file

    Returns
    -------
    True unless there are missing headers

    """
    replace_str = os.environ.get("FERMI_DIR")
    fitsfile = os.path.expandvars(fitsfile).replace("$(FERMI_DIR)", replace_str) if replace_str is not None else os.path.expandvars(fitsfile)
    hdr = pyfits.getheader(fitsfile, 0)
    required_headers = ['NAXIS1', 'CRPIX1', 'CRVAL1', 'CDELT1', 'NAXIS2', 'CRPIX2', 'CRVAL2', 'CDELT2', 'NAXIS3', 'CRPIX3', 'CRVAL2', 'CDELT3', 'CUNIT3', 'CTYPE3']
    missing_headers = set(required_headers).difference(hdr.keys())
    if len(missing_headers) > 0:
        raise HeaderCheckError("DiffuseHeaderError: {0} is missing headers ".format(fitsfile) + " ".join(missing_headers))
    else:
        return True

def getExtendedTemplate(fitsfile, extDir):
    """Searches Templates directory for extended template. Raises error if not found.

    Parameters
    ----------
    fitsfile : str
        Name of the template to search for
    extDir : str
        Full path to directory where templates are located

    Returns
    -------
    filename : str
        Full path to extended template

    """
    for root, dirs, files in os.walk(extDir):
        if fitsfile in files:
            return os.path.join(root, fitsfile)
        else:
            raise ExtendedTemplateError("TemplateError: Could not find {0} in Templates (location: {1})".format(fitsfile, root))

def getPos(ft1):
    """Get ROI information from events file header

    Parameters
    ----------
    ft1 : str
        Path to fits file containing event information

    Returns
    -------
    ra : float
        Right Ascension
    dec : float
        Declination
    rad : float
        Radius

    """
    file=pyfits.open(ft1)
    num=file[1].header['NDSKEYS']
    header=file[1].header
    right='POS(RA,DEC)'
    i=1
    keynum=0
    while i<=num:  #this step is necessary since it is not clear that the POS key word will have the same number always
        word='DSTYP%i' %i
        test=file[1].header[word]
        if(test==right):
            keynum=i
            i=num
        i+=1
    if(keynum==0):  #DSKEYS start numbering at 1, if this value hasn't been updated, KEYword doesn't exist
        print 'Error: No position keyword found in fits header (assuming position is RA and DEC.  Exiting...'
        exit()
    keyword='DSVAL%i' %keynum
    try:
        ra,dec,rad=header[keyword].strip('CIRCLE()').split(',') #gets rid of the circle and parenthesis part and splits around the comma
        float(ra)
    except:
        ra,dec,rad=header[keyword].strip('circle()').split(',')
    file.close()
    return float(ra),float(dec),float(rad)

def fileCheck(file):
    """Checks if a file exists"""
    if (not os.access(file,os.F_OK)):
        return 0
    return 1

#copied from Damien's macro
def parameter_element(free, name, maximum, minimum, scale, value):
    """Create an XML document parameter description element"""
    impl = minidom.getDOMImplementation()
    xmldoc_out = impl.createDocument(None,None,None)
    parameter = xmldoc_out.createElement('parameter')
    while float(value) >= float(maximum):
        maximum = float(maximum)*10.
    while float(value) <= float(minimum):
        minimum = float(minimum)/10.
    parameter.setAttribute('free', str(free))
    parameter.setAttribute('name', str(name))
    parameter.setAttribute('max', str(maximum))
    parameter.setAttribute('min', str(minimum))
    parameter.setAttribute('scale', str(scale))
    parameter.setAttribute('value', str(value))
    return parameter

def mybool(Input):
    return {'True':True,'False':False,'T':True,'F':False,'t':True,'f':False,'TRUE':True,'FALSE':False,"true":True,"false":False,"1":True,"0":False}.get(Input)

def angsep(xref1, yref1, xref2, yref2, input_units='degree', output_units='degree'):
    """Calculates angular separation between two points on the sky.

    Parameters
    ----------
    xref1 : float
        First horizontal coordinate
    yref1 : float
        First vertical coordinate
    xref2 : float
        Second horizontal coordinate
    yref2 : float
        Second vertical coordinate
    input_units : str (Optional)
        Units of input coordinates. Default is degree
    output_units : str (Optional)
        Units of output. Default is degree

    Returns
    -------
    sep : float
        Angular separation between (xref1, yref1) and (xref2, yref2)
    """
    xref1 *= units.Unit(input_units).to("rad")
    yref1 *= units.Unit(input_units).to("rad")
    xref2 *= units.Unit(input_units).to("rad")
    yref2 *= units.Unit(input_units).to("rad")
    return np.arccos((np.cos(yref1)*np.cos(yref2)*np.cos(xref1 - xref2)) + (np.sin(yref1)*np.sin(yref2)))*units.Unit("rad").to(output_units)