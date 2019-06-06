from xml.dom import minidom
import os

import Tools

from Exceptions import GetFluxError
from Exceptions import AddSourceError
from Exceptions import SpectrumError
from Exceptions import ExtendedTemplateError
# class AddSourceError(Exception):
#     """Raised when model cannot add a source."""
#     pass

# class SpectrumError(Exception):
#     """Raised when input spectrum type does not match allowed spectrum types"""
#     pass

def AddExtendedSource(name, spectype, spatialfunc, directory='', extDir='', ra=None, dec=None, glon=None, glat=None, major_axis=None, minor_axis=None, position_angle=None, efile='', emin=1.e2, emax=5.e5, frame='galactic', resolution='force-point', specfile=None, **spectrumargs):
    """Add an extended source to the model

    Parameters
    ----------
    name : str
        Name of the souce
    spectype : str
        Spectrum type. Options: PowerLaw, BrokenPowerLaw, LogParabola, PLSuperExpCutoff, PLSuperExpCutoff2, MapCube, FileSpectrum
    spatialfunc : str
        Spatial function describing source extension. Options: SpatialMap, RadialDisk, RadialGauss, Isotropic
    directory : str, optional
        Directory in which spectrum files will be stored
    extDir : str, optional
        Path to extended templates
    ra : float, optional
        Right Ascension of source in degrees
    dec : float, optional
        Decilination of source in degrees
    glon : float, optional
        Galactic longitude of source in degrees
    glat : float, optional
        Galactic latitude of source in degrees
    major_axis : float, optional
        The semi-major axis of the error ellipse at 95% confidence, in degrees.
    minor_axis : float, optional
        The semi-minor axis of the error ellipse at 95% confidence, in degrees.
    position_angle : float, optional
        The position angle of the 95%-confidence semi-major axis, from celestial North, positive toward increasing R.A. (eastward), in degrees.
    efile : str, optional
        Full path to extended template with which to model source. If only a name is given the extended templates directory is searched to retrieve the template.
    emin : float, optional
        Minimum energy for integrated flux. Default is 1.e2
    emax : float, optional
        Maximum energy for integrated flux. Default is 5.e5
    frame : str, optional
        Coordinate frame to use for source directions. Options: galactic, icrs. Default is galactic
    resolution : str, optional
        Resolution method if an extended template cannot be found. Options: force-point, skip, raise. Default is force-point
    specfile : str, optional
        Full path to spectrum file. If not specified the spectrum is calculated from spectral parameters.
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
    source : instance
        Document element containing source parameters
    extended : boolean
        Whether the source was modeled as an extended source or a point source

    """
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    source = xmldoc_out.createElement('source')
    source.setAttribute('name', name)
    spec = xmldoc_out.createElement('spectrum')
    spec.setAttribute('escale', 'MeV')

    # Spectral Type
    if spectype == 'PowerLaw':
        # Spatial Function
        if spatialfunc == 'RadialGauss':
            try:
                child1, child2 = GaussianSource(emin, emax, ra, dec, major_axis, minor_axis, position_angle, **spectrumargs)
                extended = True
            except SpectrumError as e:
                raise AddSourceError(e)
        elif spatialfunc == 'Isotropic':
            try:
                child1, child2 = Isotropic(emin, emax, **spectrumargs)
                extended = True
            except SpectrumError as e:
                raise AddSourceError(e)
        else:
            try:
                if os.path.isfile(efile):
                    fitsfile = efile
                else:
                    fitsfile = Tools.getExtendedTemplate(efile, extDir)
                child1, child2 = MapSource(emin, emax, fitsfile, **spectrumargs)
                extended = True
            except SpectrumError as e:
                raise AddSourceError(e)
            except ExtendedTemplateError as e:
                if resolution == 'force-point':
                    print "Could not find a template for {0}. Modeling source as a PowerLaw point source.".format(name)
                    try:
                       flux = spectrumargs['flux']
                    except KeyError:
                        try:
                            flux = Tools.getFlux(spectype, emin, emax, **spectrumargs)
                        except GetFluxError as e:
                            raise AddSourceError(e)
                    source.setAttribute('flux', "{}".format(flux))
                    try:
                        child1, child2 = GammaPointSource(spectype, emin, emax, frame=frame, ra=ra, dec=dec, glon=glon, glat=glat, **spectrumargs)
                        extended = False
                    except SpectrumError as e:
                        raise AddSourceError(e)
                elif resolution == 'skip':
                    raise AddSourceError(e)
                elif resolution == 'raise':
                    raise e
                else:
                    raise AddSourceError(e)
    elif spectype in ['LogParabola', 'PLSuperExpCutoff', 'PLSuperExpCutoff2']:
        try:
            if os.path.isfile(efile):
                fitsfile = efile
            else:
                fitsfile = Tools.getExtendedTemplate(efile, extDir)
            child1, child2 = FileSpectrumMap(name, spectype, emin, emax, fitsfile, directory, specfile=specfile, **spectrumargs)
            extended = True
        except SpectrumError as e:
            raise AddSourceError(e)
        except ExtendedTemplateError as e:
            if resolution == 'force-point':
                print "Could not find a template for {0}. Modeling source as a FileSpectrum point source.".format(name)
                try:
                    child1, child2 = FileSpectrum(name, spectype, emin, emax, directory, frame=frame, ra=ra, dec=dec, glon=glon, glat=glat, specfile=specfile, **spectrumargs)
                    extended = False
                except SpectrumError as e:
                    raise AddSourceError(e)
            elif resolution == 'skip':
                raise AddSourceError(e)
            elif resolutioin == 'raise':
                raise e
            else:
                raise AddSourceError(e)
    elif spectype == "MapCube":
        try:
            if os.path.isfile(efile):
                fitsfile = efile
            else:
                fitsfile = Tools.getExtendedTemplate(efile, extDir)
            child1, child2 = MapCube(spectype, emin, emax, fitsfile, directory, **spectrumargs)
            extended = True
        except SpectrumError as e:
            raise AddSourceError(e)
        except ExtendedTemplateError as e:
            raise AddSourceError(e)
    elif spectype == "BrokenPowerLaw":
        try:
            if os.path.isfile(efile):
                fitsfile = efile
            else:
                fitsfile = Tools.getExtendedTemplate(efile, extDir)
            child1, child2 = FileSpectrumMap(name, spectype, emin, emax, fitsfile, directory, specfile=specfile, **spectrumargs)
            extended = True
        except SpectrumError as e:
            raise AddSourceError(e)
        except ExtendedTemplateError as e:
            if resolution == 'force-point':
                print "Could not find a template for {0}. Modeling source as a BrokenPowerLaw point source.".format(name)
                try:
                    child1, child2 = GammaPointSource(spectype, emin, emax, frame=frame, ra=ra, dec=dec, glon=glon, glat=glat, **spectrumargs)
                    extended = False
                except SpectrumError as e:
                        raise AddSourceError(e)
            elif resolution == 'skip':
                raise AddSourceError(e)
            elif resolution == 'raise':
                raise e
            else:
                raise AddSourceError(e)
    elif spectype == "FileSpectrum":
        try:
            if os.path.isfile(efile):
                fitsfile = efile
            else:
                fitsfile = Tools.getExtendedTemplate(efile, extDir)
            child1, child2 = FileSpectrumMap(name, spectype, emin, emax, fitsfile, directory, specfile=specfile, **spectrumargs)
            extended = True
        except SpectrumError as e:
            raise AddSourceError(e)
        except ExtendedTemplateError as e:
            if resolution == 'force-point':
                print "Could not find a template for {0}. Modeling source as a FileSpectrum point source.".format(name)
                try:
                    child1, child2 = FileSpectrum(name, spectype, emin, emax, directory, frame=frame, ra=ra, dec=dec, glon=glon, glat=glat, specfile=specfile, **spectrumargs)
                    extended = False
                except SpectrumError as e:
                    raise AddSourceError(e)
            elif resolution == 'skip':
                raise AddSourceError(e)
            elif resolutioin == 'raise':
                raise e
            else:
                raise AddSourceError(e)
    else:
        raise AddSourceError("Cannot add extended source {0} with spectrum type {1} to model.".format(name, spectype))

    spec.appendChild(child1)
    spec.appendChild(child2)
    source.appendChild(spec)

    return source, extended

def AddPointSource(name, spectype, emin, emax, directory, ra=None, dec=None, glon=None, glat=None, frame='galactic', specfile=None, **spectrumargs):
    """Add a point source to the model

    Parameters
    ----------
    name : str
        Name of the souce
    spectype : str
        The spectral type in the global model (PowerLaw, LogParabola, PLSuperExpCutoff).
    emin : float
        Minimum energy for integrated flux
    emax : float
        Maximum energy for integrated flux
    directory : str
        Directory in which to save model files
    ra : float, optional
        Right Ascension of source in degrees
    dec : float, optional
        Decilination of source in degrees
    glon : float, optional
        Galactic longitude of source in degrees
    glat : float, optional
        Galactic latitude of source in degrees
    frame : str, optional
        Coordinate frame to use for source directions. Options: galactic, icrs.
    specfile : str, optional
        Full path to spectrum file. If not specified the spectrum is calculated from spectral parameters.
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
    source : instance
        Document element containing source parameters

    """
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    source = xmldoc_out.createElement('source')
    source.setAttribute('name', name)
    spec = xmldoc_out.createElement('spectrum')
    spec.setAttribute('escale', 'MeV')

    if (spectype == 'PowerLaw') or (spectype == 'Monochromatic') or (spectype == 'BrokenPowerLaw'):
        try:
            flux = spectrumargs['flux']
        except KeyError:
            try:
                flux = Tools.getFlux(spectype, emin, emax, **spectrumargs)
            except GetFluxError as e:
                raise AddSourceError(e)
        source.setAttribute('flux', "{}".format(flux))
        try:
            child1, child2 = GammaPointSource(spectype, emin, emax, frame=frame, ra=ra, dec=dec, glon=glon, glat=glat, **spectrumargs)
        except SpectrumError as e:
            raise AddSourceError(e)
    elif spectype in ['LogParabola', 'PLSuperExpCutoff2', 'PLSuperExpCutoff', 'FileSpectrum']:
        try:
            child1, child2 = FileSpectrum(name, spectype, emin, emax, directory, frame=frame, ra=ra, dec=dec, glon=glon, glat=glat, specfile=specfile, **spectrumargs)
        except SpectrumError as e:
            raise AddSourceError(e)
    else:
        raise AddSourceError("Cannot add point source {0} with spectrum type {1} to model.".format(name, spectype))

    spec.appendChild(child1)
    spec.appendChild(child2)
    source.appendChild(spec)

    return source

def MapSource(emin, emax, fitsfile, **spectrumargs):
    """Observation Simulation Map Source flux definition"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)

    try:
        flux = spectrumargs['flux']
    except KeyError:
        try:
            flux = Tools.getFlux('PowerLaw', emin, emax, **spectrumargs)
        except GetFluxError as e:
            raise SpectrumError(e)

    try:
        pli = spectrumargs['pl_index']
    except KeyError:
        raise SpectrumError("MapSourceError: Cannot create MapSource without power law index (pl_index).")

    spectrumClass = xmldoc_out.createElement('SpectrumClass')
    spectrumClass.setAttribute('name', 'MapSource')
    spectrumClass.setAttribute('params',"{0},{1},{2},{3},{4}".format(flux, pli, fitsfile, emin, emax))
    
    use_spectrum = xmldoc_out.createElement('use_spectrum')
    use_spectrum.setAttribute('frame', "galaxy")

    return spectrumClass, use_spectrum

def FileSpectrumMap(name, spectype, emin, emax, fitsfile, directory, specfile=None, **spectrumargs):
    """Observation Simulation FileSpectrumMap flux definition

    This flux definition takes all FileSpectrum parameters as well as the MapSource parameters
    """
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)

    if specfile is None:
        try:
            filename = Tools.writeSpectrum(name, spectype, emin, emax, directory, **spectrumargs)
        except WriteSpectrumError as e:
            raise SpectrumError(e)
    else:
        filename = specfile
        
    try:
        flux = spectrumargs['flux']
    except KeyError:
        try:
            flux = Tools.getFlux(spectype, emin, emax, **spectrumargs)
        except GetFluxError as e:
            raise SpectrumError(e)


    spectrumClass = xmldoc_out.createElement('SpectrumClass')
    spectrumClass.setAttribute('name', 'FileSpectrumMap')
    spectrumClass.setAttribute('params', "flux={0},fitsFile={1},specFile={2},emin={3},emax={4}".format(flux, fitsfile, filename, emin, emax))
    
    use_spectrum = xmldoc_out.createElement('use_spectrum')
    use_spectrum.setAttribute('frame', "galaxy")

    return spectrumClass, use_spectrum

def MapCube(spectype, emin, emax, efile, directory, **spectrumargs):
    """Observation Simulation MapCube flux definition"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    try:
        flux = spectrumargs['flux']
    except KeyError:
        try:
            flux = Tools.getFlux(spectype, emin, emax, **spectrumargs)
        except GetFluxError as e:
            raise SpectrumError(e)

    spectrumClass = xmldoc_out.createElement('SpectrumClass')
    spectrumClass.setAttribute('name', "MapCube")
    spectrumClass.setAttribute('params', "{0},{1}".format(flux, efile))

    use_spectrum = xmldoc_out.createElement('use_spectrum')
    use_spectrum.setAttribute('frame', "galaxy")

    return spectrumClass, use_spectrum

def GaussianSource(emin, emax, ra, dec, major_axis, minor_axis, position_angle, **spectrumargs):
    """Observation Simulation GaussianSource flux definition"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)

    try:
        flux = spectrumargs['flux']
    except KeyError:
        try:
            flux = Tools.getFlux("PowerLaw", emin, emax, **spectrumargs)
        except GetFluxError as e:
            raise SpectrumError(e)

    try:
        pli = spectrumargs['pl_index']
    except KeyError:
        raise SpectrumError("GaussianSourceError: GaussianSource source must include power law index (pl_index)")

    spectrumClass = xmldoc_out.createElement('SpectrumClass')
    spectrumClass.setAttribute('name', "GaussianSource")
    spectrumClass.setAttribute('params', "{0},{1},{2},{3},{4},{5},{6},{7},{8}".format(flux, pli, ra, dec, major_axis, minor_axis, position_angle, emin, emax))

    use_spectrum = xmldoc_out.createElement("use_spectrum")
    use_spectrum.setAttribute('frame', "galaxy")

    return spectrumClass, use_spectrum

def Isotropic(emin, emax, **spectrumargs):
    """Observation Simulation Isotropic flux definition"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    try:
        flux = spectrumargs['flux']
    except KeyError:
        try:
            flux = Tools.getFlux("PowerLaw", emin, emax, **spectrumargs)
        except GetFluxError as e:
            raise SpectrumError(e)

    spectrumClass = xmldoc_out.createElement("SpectrumClass")
    spectrumClass.setAttribute("name", "Isotropic")
    spectrumClass.setAttribute("params", "{0},{1},{2},{3}".format(flux, pli, emin, emax))

    use_spectrum = xmldoc_out.createElement("use_spectrum")
    use_spectrum.setAttribute("frame", "galaxy")

    return spectrumClass, use_spectrum

def GammaPointSource(spectype, emin, emax, frame='galactic', ra=None, dec=None, glon=None, glat=None, **spectrumargs):
    """Gamma-ray point source flux definition"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)

    particle = xmldoc_out.createElement("particle")
    particle.setAttribute("name", "gamma")

    if spectype == "Monochromatic":
        try:
            e = spectrumargs['energy']
        except KeyError:
            raise SpectrumError("GammaPointSourceError: Energy must be supplied in order to add a monochromatic source.")
            
        energy = xmldoc_out.createElement("energy")
        energy.setAttribute("e", "{}".format(e))
        particle.appendChild(energy)

    elif spectype == "PowerLaw":
        power_law = xmldoc_out.createElement("power_law")
        try:
            pli = spectrumargs['pl_index']
        except KeyError:
            raise SpectrumError("GammaPointSourceError: PowerLaw source must include power law index (pl_index)")

        power_law.setAttribute("emin", "{}".format(emin))
        power_law.setAttribute("emax", "{}".format(emax))
        power_law.setAttribute("gamma", "{}".format(pli))
        particle.appendChild(power_law)

    elif spectype == "BrokenPowerLaw":
        power_law = xmldoc_out.createElement("power_law")
        try:
            pli = spectrumargs['pl_index']
            gamma2 = spectrumargs['gamma2']
            ebreak = spectrumargs['ebreak']
        except KeyError:
            raise SpectrumError("GammaPointSourceError: BrokenPowerLaw source must include low energy power law index (pl_index), high energy power law index (gamma2), and break energy (ebreak).")

        power_law.setAttribute("emin", "{}".format(emin))
        power_law.setAttribute("emax", "{}".format(emax))
        power_law.setAttribute("gamma", "{}".format(pli))
        power_law.setAttribute("ebreak", "{}".format(ebreak))
        power_law.setAttribute("gamma2", "{}".format(gamma2))
        particle.appendChild(power_law)
    else:
        raise SpectrumError("GammaPointSourceError: Cannot create source with spectrum type {0}".format(spectype))

    if frame == 'galactic':
        direction = xmldoc_out.createElement("galactic_dir")
        direction.setAttribute("l", "{}".format(glon))
        direction.setAttribute("b", "{}".format(glat))
    elif frame in ['icrs', 'fk5']:
        direction = xmldoc_out.createElement("celestial_dir")
        direction.setAttribute("ra", "{}".format(ra))
        direction.setAttribute("dec", "{}".format(dec))
    else:
        raise SpectrumError("GammaPointSourceError: Cannot create a direction with the frame {}".format(frame))

    return particle, direction

def FileSpectrum(name, spectype, emin, emax, directory, frame='galactic', ra=None, dec=None, glon=None, glat=None, specfile=None, **spectrumargs):
    """Observation Simulation FileSpectrum flux definition"""
    xmldoc_out = minidom.getDOMImplementation().createDocument(None, None, None)
    
    try:
        flux = spectrumargs['flux']
    except KeyError:
        try:
            flux = Tools.getFlux(spectype, emin, emax, **spectrumargs)
        except GetFluxError as e:
            raise SpectrumError(e)

    if specfile is None:
        try:
            filename = Tools.writeSpectrum(name, spectype, emin, emax, directory, **spectrumargs)
        except WriteSpectrumError as e:
            raise SpectrumError(e)
    else:
        filename = specfile

    spectrumClass = xmldoc_out.createElement("SpectrumClass")
    spectrumClass.setAttribute("name", "FileSpectrum")
    spectrumClass.setAttribute("params", "flux={0},specFile={1}".format(flux, filename))

    if frame == 'galactic':
        direction = xmldoc_out.createElement("galactic_dir")
        direction.setAttribute("l", "{}".format(glon))
        direction.setAttribute("b", "{}".format(glat))
    elif frame in ['icrs', 'fk5']:
        direction = xmldoc_out.createElement("celestial_dir")
        direction.setAttribute("ra", "{}".format(ra))
        direction.setAttribute("dec", "{}".format(dec))
    else:
        raise SpectrumError("FileSpectrumError: Cannot create a direction with the frame {}".format(frame))

    return spectrumClass, direction