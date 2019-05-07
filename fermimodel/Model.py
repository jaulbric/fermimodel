import os
from astropy.coordinates import SkyCoord
from xml.dom import minidom

import Tools
import AddCatalogSources
import SimulationSources
import BuildRegion

class model:
    def __init__(self, name='mymodel', eventsfile=None, catalog=None, out=None, roi=None, frame='fk5', unit='degree', allsky=False, model_type='likelihood'):
        """Model class containing methods and attributes for building models for use with Fermitools

        Parameters
        ----------
        name : str (Optional)
            Name of the model. Default is 'mymodel'.
        eventsfile : str (Optional)
            Path to file containing events data. This is only used for region of interest information
        catalog : str (Optional)
            Path to source catalog.
        out : str (Optional)
            Name of the output model files.
        roi : tuple (Optional)
            Tuple defining ROI (horizontal center, vertical center, radius).
        frame : str (Optional)
            Coordinate frame to use for region of interest argument, e.g. galactic, icrs. Any astropy frame is allowed
        unit : str (Optional)
            Units for region of interest argument. Any astropy unit is allowed.
        allsky : bool (Optional)
            Flag to set region of interest to the entire sky. By default this sets the region of interest center to the galactic center (l=0, b=0).

        """
        self.name = name
        self.srcs = catalog

        if out is None:
            self.out = os.path.join(os.getcwd(), self.name + '.xml')
        else:
            self.out = out
        
        if roi is not None:
            self.setROI(roi=roi, frame=frame, unit=unit, allsky=allsky)
        elif eventsfile is not None:
            self.roi = Tools.getPos(eventsfile)
        else:
            self.roi = None
        
        if not model_type in ['likelihood', 'simulation']:
            raise IOError("Model type must be 'likelihood' or 'simulation'.")
        else:
            self.model_type = model_type

    def setROI(self, roi=None, frame='icrs', unit='degree', allsky=False):
        if allsky:
            c = SkyCoord(0., 0., frame='galactic', unit='degree')
            self.roi = (c.fk5.ra.degree, c.fk5.dec.degree, 180.)
        elif roi is not None:
            c = SkyCoord(roi[0], roi[1], frame=frame, unit=unit)
            self.roi = (c.fk5.ra.degree, c.fk5.dec.degree, roi[2])
        else:
            raise IOError("Could not define the region of interest.")

    def Print(self):
        """Print self"""
        excluded_keys = ['model', 'Sources']
        for key, val in self.__dict__.items():
            if key in excluded_keys:
                continue
            elif key == 'srcs':
                print 'Source list file: ', val
            elif key == 'out':
                print 'Output file name: ', val
            elif key == 'srcs':
                print 'Catalog : ', 
            elif key == 'model_type':
                print 'Model Type :', val
            else:
                print "{0} : {1}".format(key, val)

    def loadCatalog(self, catalog=None, GDfile=None, GDname='GalacticDiffuse', ISOfile=None, ISOname='IsotropicDiffuse', normsOnly=False, extDir='', radLim=-1, maxRad=None, ExtraRad=10., sigFree=5., varFree=True, psForce=False, E2CAT=False, makeRegion=False, GIndexFree=False, wd=None, oldNames=False, emin=1.e2, emax=5.e5, frame='fk5', extSrcRes='force-point', apply_mask=False, GDflux=8.4631675, ISOflux=0.):
        """Include sources in the catalog to the model. Optionaly include a galactic and isotropic diffuse models.

        Parameters
        ----------
        catalog : str
            Path to source catalog
        GDfile : str (Optional)
            Path to galactic diffuse emission file
        GDname : str (Optional)
            Name of galactic diffuse emission model
        ISOfile : str (Optional)
            Path to isotropic diffuse emission file
        ISOname : str (Optional)
            Name of isotropic diffuse emission file
        normsOnly : bool (Optional)
            Flag to only set normalization parameters free
        extDir : str (Optional)
            Directory with extended source templates
        radLim : float (Optional)
            Radius in degrees from center of ROI beyond which source parameters are fixed
        maxRad : float (Optional)
            Absolute maximum radius beyond which sources are fixed. This may be necessary when doing binned analysis and a variable source beyond radLim would be set free but this source is beyond the boundaries of the square region used for the binned likelihood
        ExtraRad : float (Optional)
            Radius beyond ROI radius out to which sources will be included with fixed parameters. Default of 10 degrees is good for analysis starting around 100 MeV, but for higher energy fits this can be decreased.
        sigFree : float (Optional)
            Significance below which source parameters are fixed, even if within radLim.
        varFree : float (Optional)
            Variability index above which source parameters are free. If beyond radLim and/or below sigFree only the normalization parameters is set free. Currently not implemented for building from xml catalog.
        psForce : bool (Optional)
            Flag to force exentended sources to be point sources
        E2CAT : bool (Optional)
            Flag to force use of catalog names for extended sources (only matters is using catalog FITS file).
        makeRegion : bool (Optional)
            Flag to also generate ds9 region file.
        GIndexFree : bool (Optional)
            The galactic diffuse is given a power-law spectral shape but by default the index is frozen. Setting this flag to True allows that to be free for additional freedom in diffuse fit.
        wd : str (Optional)
            Path to directory in which output files will be saved. If an absolute path for the output file names is given it will override this argument.
        oldNames : bool (Optional)
            Sets use of old naming convention. Underscore before name and no spaces. Default is False for likelihood models. This is automatically set to true for simulation models.
        emin : float (Optional)
            Minimum energy in MeV for source simulation. This should match simulation criteria. Default is 100 MeV.
        emax : float (Optional)
            Maximum energy in MeV for source simulation. This should match simulation criteria. Default is 5.e5 MeV.
        GDflux : float (Optional)
            Integrated flux to use for the galactic diffuse emission model in photons/m^2/s. Default is 8.4631675.
        ISOflux : float (Optional)
            Integrated flux to use for the isotropic diffuse emission model in photons/m^2/s. If 0.0 flux is calculated on the fly by the simulator. Default is 0.0.

        Returns
        -------

        """
        if (catalog is None) and (self.srcs is None):
            raise IOError("No catalog input")
        elif catalog is None:
            catalog = self.srcs
        

        self.var = varFree
        self.psF = psForce
        self.E2C = E2CAT
        self.nO = normsOnly
        extDir = (extDir if extDir != '' else '$(FERMI_DIR)/data/pyBurstAnalysisGUI/templates')
        self.extD = extDir
        self.ER = ExtraRad
        self.sig = sigFree
        self.reg = makeRegion
        self.GIF = GIndexFree
        self.wd = os.path.abspath(wd)
        self.extSrcRes = extSrcRes
        if makeRegion:
            rhold = os.path.splitext(os.path.basename(self.name))[0]
            self.regFile = os.path.join(self.wd, 'ROI_' + rhold + '.reg')
        print 'Creating file and adding sources from Catalog {0}'.format(catalog)

        if self.model_type == 'simulation':
            self.model, self.Sources = AddCatalogSources.simulation(srcs=catalog, roi=self.roi, var=self.var, psF=self.psF, E2C=self.E2C, nO=self.nO, extDir=extDir, ER=self.ER, sig=self.sig, reg=self.reg, GIF=self.GIF, wd=self.wd, extSrcRes=self.extSrcRes, GDfile=GDfile, GDname=GDname, ISOfile=ISOfile, ISOname=ISOname, oldNames=oldNames, emin=emin, emax=emax, frame=frame, apply_mask=apply_mask, GDflux=GDflux, ISOflux=ISOflux)
            self.srcs = catalog
        elif self.model_type == 'likelihood':
            self.model, self.Sources = AddCatalogSources.likelihood(srcs=catalog, roi=self.roi, var=self.var, psF=self.psF, E2C=self.E2C, nO=self.nO, extDir=extDir, ER=self.ER, sig=self.sig, reg=self.reg, GIF=self.GIF, wd=self.wd, GDfile=GDfile, GDname=GDname, ISOfile=ISOfile, ISOname=ISOname, oldNames=oldNames, frame=frame)
            self.srcs = catalog
        else:
            raise IOError("Model type must be either simulation or likelihood. To {0}".format(self.model_type))

    def writeXML(self, directory=None, out=None):
        """Write the model to XML

        Parameters
        ----------
        directory : str, optional
            Path to directory in which the XML file will be saved. Defaults is the work directory of srcList
        out : str, optional
            Name of the XML file. Default is output name of srcList

        Returns
        -------
        filename : str

        """
        if directory is None:
            directory = self.wd
        if out is None:
            if os.path.isabs(self.out):
                out = os.path.basename(self.out)
            else:
                out = self.out

        filename = os.path.join(directory,out)

        if self.model_type == 'simulation':
            try:
                source_library = self.model.firstChild
                source_library.writexml(open(filename, "wb"), indent="", addindent="\t", newl="\n")
                # print self.model.toprettyxml()
                # for node in self.model.childNodes:
                    # node.writexml(open(filename, "wb"), indent="", addindent="  ", newl="\n")
                return filename
            except AttributeError as e:
                if not hasattr(self, "model"):
                    print "There is no model to write to file."
                    return None
                else:
                    raise e
        elif self.model_type == 'likelihood':
            try:
                self.model.writexml(open(filename, "wb"), indent="", addindent="\t", newl="\n")
                return filename
            except AttributeError as e:
                if not hasattr(self, "model"):
                    print "There is no model to write to file."
                    return None
                else:
                    raise e

    def writeSrcList(self, directory=None, out=None):
        """Write source list to file

        Parameters
        ----------
        directory : str
            Path to directory in which the XML file will be saved.
        out : str
            Name of the XML file

        Returns
        -------

        """
        if directory is None:
            directory = self.wd
        if out is None:
            if os.path.isabs(self.out):
                out = os.path.splitext(os.path.basename(self.out))[0] + ".txt"
            else:
                out = os.path.splitext(self.out)[0] + ".txt"

        filename = os.path.join(directory,out)

        with open(filename, "wb") as f:
            try:
                f.write('\n'.join(self.Sources.keys()))
                return filename
            except AttributeError as e:
                if not hasattr(self, "Sources"):
                    print "There are no sources in the model."
                    return None
                else:
                    raise e

    def addSource(self, name, spectrum_type, ra=None, dec=None, glon=None, glat=None, major_axis=None, minor_axis=None, position_angle=None, spatial_function=None, extended_template=None, emin=1.e2, emax=5.e5, frame='galactic', specFile=None, **spectrumargs):
        """Add a single source to the model.

        Parameters
        ----------
        name : str
            Name of the source to add.
        spectrum_type : str
            Spectral type of the source. Options: Monochromatic, PowerLaw, BrokenPowerLaw, LogParabola, PLSuperExpCutoff2, MapCube, FileSpectrum
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
        spatial_function : str, optional
            Spatial function describing source extension. Options: None, SpatialMap, RadialDisk, RadialGauss, Isotropic
        extended_template : str, optional
            Full path to extended template.
        emin : float, optional
            Minimum energy for integrated flux. Default is 1.e2 MeV
        emax : float, optional
            Maximum energy for integrated flux. Default is 5.e5 MeV
        frame : str, optional
            Coordinate frame to use for source directions. Options: galactic, icrs
        specFile : str, optional
            Name (and path) to ASCII file containing spectrum. The first column contains the energies, and second column contains the differential flux at those energies. Energy scale should be in MeV.
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
        flux : float, optional
            Integrated flux in photons/cm^2/s. Default is to calculate the flux from spectral parameters.

        Returns
        -------

        """
        allowed_spectypes = ["Monochromatic", "PowerLaw", "BrokenPowerLaw", "LogParabola", "PLSuperExpCutoff2", "MapCube", "FileSpectrum"]
        if spectrum_type not in allowed_spectypes:
            raise AddSourceError("Cannot add a source with spectrum type {0}. Options are: {1}".format(spectrum_type, ", ".join(allowed_spectypes)))

        if spatial_function is None:
            source = SimulationSources.AddPointSource(name, spectrum_type, emin, emax, self.wd, ra=ra, dec=dec, glon=glon, glat=glat, frame=frame, specfile=specFile, **spectrumargs)
            try:
                self.model.appendChild(source)
                print "Added point source {0} to model.".format(name)
            except AttributeError as e:
                if not hasattr(self, "model"):
                    self.model = minidom.getDOMImplementation().createDocument(None,'source_library',None)
                    self.model.documentElement.setAttribute('title', 'source_library')
                    self.model.appendChild(source)
                    print "Added point source {0} to model.".format(name)
                else:
                    raise e

        else:
            source, modeled_extended = AddExtendedSource(name, spectrum_type, spatial_function, directory=self.wd, extDir=self.extD, ra=ra, dec=dec, glon=glon, glat=glat, major_axis=major_axis, minor_axis=minor_axis, position_angle=position_angle, efile=extended_template, emin=emin, emax=emax, frame=frame, resolution='skip', specfile=specFile, **spectrumargs)
            try:
                self.model.appendChild(source)
                print "Added extended source {0} to model.".format(name)
            except AttributeError as e:
                if not hasattr(self, "model"):
                    self.model = minidom.getDOMImplementation().createDocument(None,'source_library',None)
                    self.model.documentElement.setAttribute('title', 'source_library')
                    self.model.appendChild(source)
                    print "Added point source {0} to model.".format(name)
                else:
                    raise e

        print "Run writeXML and writeSrcList to update model files."

    def removeSource(self, *name):
        """Remove source from the model.

        Parameters
        ----------
        name : str
            Variable number of source names to remove.

        Returns
        -------

        """
        for n in name:
            found = False
            try:
                for src in self.model.getElementsByTagName("source"):
                    if src.getAttribute("name") == n:
                        self.model.documentElement.removeChild(src)
                        self.Sources.pop(n)
                        found = True
                        print "Removed {0} from model.".format(n)

                if not found:
                    print "Source {0} is not in model.".format(n)
            except AttributeError as e:
                if not hasattr(self, "model"):
                    print "Model has not been built, so there are not sources to remove."
                elif not hasattr(self, "Sources"):
                    print "No source list."
                else:
                    raise e

        print "Run writeXML and writeSrcList to update model files."

    def BuildRegion(self, regFile=None, frame='fk5'):
        """Build the DS9 region file for the model

        Parameters
        ----------
        regFile : str
            Name of the output file

        Returns
        -------
        regFile : str
            Name of the output file

        """
        if regFile is None:
            try:
                regFile = self.regFile
            except AttributeError:
                if hasattr(self, "wd"):
                    regFile = os.path.join(self.wd, self.name + '.reg')
                else:
                    regFile = self.name + '.reg'

        try:
            out = BuildRegion.BuildRegion(regFile, self.Sources, model_type=self.model_type, frame=frame)
        except AttributeError as e:
            if not hasattr(self, "Sources"):
                print "No sources in the model to include in the region file."
            else:
                raise e

        return out
