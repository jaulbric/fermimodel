# FermiModel

Documentation is in progess.

This package is for use with the [Fermitools](https://github.com/fermi-lat/Fermitools-conda/wiki/Quickstart-Guide). The main function of the package is to produce model files for use in analysis with gtlike and simulations with gtobssim. This is essentially an extension of [make4FGLxml.py](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/make4FGLxml.py), written by Tyrel Johnson.

## Installation

Download the package

```bash
# Get the version of the latest release
version=`curl -s https://api.github.com/repos/jaulbric/fermimodel/releases/latest | grep "tag_name" | cut -d ":" -f 2 | tr -d " "-\"-,`
# download and untar the package
wget https://github.com/jaulbric/fermimodel/archive/$version.tar.gz && mkdir /fermimodel && tar xvfz $version.tar.gz -C /fermimodel --strip-components 1 && rm $version.tar.gz
```

Move to the package root directory and install

```bash
cd /path/to/FermiModel
python setup.py install
```
This will also generate an executable `fermimodel` and place it in the user's bin.

## Usage

### As Python Module

The package can be imported in python as usual

```python
import fermimodel
```
Models can then be built by creating an instance of the model class

```python
mymodel = fermimodel.model()
```

> fermimodel.**model(name='mymodel', eventsfile=None, catalog=None, out=None, roi=None, frame='fk5', unit='degree', allsky=False, model_type='likelihood')**
>
> Model class containing methods and attributes for building models for use with Fermitools
>
> Parameters
>
>     name : str (Optional)
>         Name of the model. Default is 'mymodel'.
>     eventsfile : str (Optional)
>         Path to file containing events data. This is only used for region of interest information
>     catalog : str (Optional)
>         Path to source catalog.
>     out : str (Optional)
>         Name of the output model files.
>     roi : tuple (Optional)
>         Tuple defining ROI (horizontal center, vertical center, radius).
>     frame : str (Optional)
>         Coordinate frame to use for region of interest argument, e.g. galactic, icrs. Any astropy frame is allowed
>     unit : str (Optional)
>         Units for region of interest argument. Any astropy unit is allowed.
>     allsky : bool (Optional)
>         Flag to set region of interest to the entire sky. By default this sets the region of interest center to the galactic center (l=0, b=0).

Then load a catalog of sources or add sources individually (currently only for simulation models).

```python
mymodel.loadCatalog()
```

>fermimodel.model.**loadCatalog(catalog=None, GDfile=None, GDname='GalacticDiffuse', ISOfile=None, ISOname='IsotropicDiffuse', ISOpath="$(FERMI_DIR)/refdata/fermi/galdiffuse/isotropic_allsky.fits", normsOnly=False, extDir='', radLim=-1, maxRad=None, ExtraRad=10., sigFree=5., varFree=True, psForce=False, E2CAT=False, makeRegion=False, GIndexFree=False, wd=None, oldNames=False, emin=1.e2, emax=5.e5, frame='fk5', extSrcRes='force-point', apply_mask=False, GDflux=8.4631675, ISOflux=0.)**
>
>Include sources in the catalog to the model. Optionaly include a galactic and isotropic diffuse models.
>
> Parameters
>
>     catalog : str
>         Path to source catalog
>     GDfile : str (Optional)
>         Path to galactic diffuse emission file
>     GDname : str (Optional)
>         Name of galactic diffuse emission model
>     ISOfile : str (Optional)
>         Path to isotropic diffuse emission spectrum file
>     ISOname : str (Optional)
>         Name of isotropic diffuse emission file
>     ISOpath : str (Optional)
>         Path to isotropic diffuse spatial file
>     normsOnly : bool (Optional)
>         Flag to only set normalization parameters free
>     extDir : str (Optional)
>         Directory with extended source templates
>     radLim : float (Optional)
>         Radius in degrees from center of ROI beyond which source parameters are fixed
>     maxRad : float (Optional)
>         Absolute maximum radius beyond which sources are fixed. This may be necessary when doing binned analysis and a variable source beyond radLim would be set free but this source is beyond the boundaries of the square region used for the binned likelihood
>     ExtraRad : float (Optional)
>         Radius beyond ROI radius out to which sources will be included with fixed parameters. Default of 10 degrees is good for analysis starting around 100 MeV, but for higher energy fits this can be decreased.
>     sigFree : float (Optional)
>         Significance below which source parameters are fixed, even if within radLim.
>     varFree : float (Optional)
>         Variability index above which source parameters are free. If beyond radLim and/or below sigFree only the normalization parameters is set free. Currently not implemented for building from xml catalog.
>     psForce : bool (Optional)
>         Flag to force exentended sources to be point sources
>     E2CAT : bool (Optional)
>         Flag to force use of catalog names for extended sources (only matters is using catalog FITS file).
>     makeRegion : bool (Optional)
>         Flag to also generate ds9 region file.
>     GIndexFree : bool (Optional)
>         The galactic diffuse is given a power-law spectral shape but by default the index is frozen. Setting this flag to True allows that to be free for additional freedom in diffuse fit.
>     wd : str (Optional)
>         Path to directory in which output files will be saved. If an absolute path for the output file names is given it will override this argument.
>     oldNames : bool (Optional)
>         Sets use of old naming convention. Underscore before name and no spaces. Default is False for likelihood models. This is automatically set to true for simulation models.
>     emin : float (Optional)
>         Minimum energy in MeV for source simulation. This should match simulation criteria. Default is 100 MeV.
>     emax : float (Optional)
>         Maximum energy in MeV for source simulation. This should match simulation criteria. Default is 5.e5 MeV.
>     GDflux : float (Optional)
>         Integrated flux to use for the galactic diffuse emission model in photons/m^2/s. Default is 8.4631675.
>     ISOflux : float (Optional)
>         Integrated flux to use for the isotropic diffuse emission model in photons/m^2/s. If 0.0 flux is calculated on the fly by the simulator. Default is 0.0.

>   Returns

```python
mymodel.addSource()
```

> fermimodel.model.**addSource(name, spectrum_type, ra=None, dec=None, glon=None, glat=None, major_axis=None, minor_axis=None, position_angle=None, spatial_function=None, extended_template=None, emin=1.e2, emax=5.e5, frame='galactic', specFile=None, \*\*spectrumargs)**
>
> Add a single source to the model.
>
> Parameters
>
>     name : str
>         Name of the source to add.
>     spectrum_type : str
>         Spectral type of the source. Options: Monochromatic, PowerLaw, BrokenPowerLaw, LogParabola, PLSuperExpCutoff2, MapCube, FileSpectrum
>     ra : float, optional
>         Right Ascension of source in degrees
>     dec : float, optional
>         Decilination of source in degrees
>     glon : float, optional
>         Galactic longitude of source in degrees
>     glat : float, optional
>         Galactic latitude of source in degrees
>     major_axis : float, optional
>         The semi-major axis of the error ellipse at 95% confidence, in degrees.
>     minor_axis : float, optional
>         The semi-minor axis of the error ellipse at 95% confidence, in degrees.
>     position_angle : float, optional
>         The position angle of the 95%-confidence semi-major axis, from celestial North, positive toward increasing R.A. (eastward), in degrees.
>     spatial_function : str, optional
>         Spatial function describing source extension. Options: None, SpatialMap, RadialDisk, RadialGauss, Isotropic
>     extended_template : str, optional
>         Full path to extended template.
>     emin : float, optional
>         Minimum energy for integrated flux. Default is 1.e2 MeV
>     emax : float, optional
>         Maximum energy for integrated flux. Default is 5.e5 MeV
>     frame : str, optional
>         Coordinate frame to use for source directions. Options: galactic, icrs
>     specFile : str, optional
>         Name (and path) to ASCII file containing spectrum. The first column contains the energies, and second column contains the differential flux at those energies. Energy scale should be in MeV.

> \*\*spectrumargs
>
>     pl_flux_density : float, optional
>         The differential flux at pivot_energy for the PowerLaw fit, in photons/cm2/MeV/s.
>     lp_flux_density : float, optional
>         The differential flux at pivot_energy for the LogParabola fit, in photons/cm2/MeV/s.
>     plec_flux_density : float, optional
>         The differential flux at pivot_energy for the PLSuperExpCutoff fit, in photons/cm2/MeV/s.
>     pivot_energy : float, optional
>         The energy, in MeV, at which the error in the differential photon flux is minimal (i.e., the decorrelation energy for the power-law fit).
>     pl_index : float, optional
>         The photon index for the PowerLaw fit.
>     lp_index : float, optional
>         Photon Index at Pivot Energy for LogParabola fit.
>     lp_beta : float, optional
>         The curvature parameter (beta) for LogParabola spectrum types.
>     plec_index : float, optional
>         The low energy photon index for PLSuperExpCutoff fit.
>     plec_expfactor : float, optional
>         The exponential factor for PLSuperExpCutoff fit.
>     plec_exp_index : float, optional
>         The exponential index for PLSuperExpCutoff fit.
>     flux : float, optional
>         Integrated flux in photons/cm^2/s. Default is to calculate the flux from spectral parameters.

> Returns

Once the model has been built an XML file can be produced

```python
xmlfile_path = mymodel.writeXML()
```

> fermimodel.model.**writeXML(directory=None, out=None)**
>
> Write the model to XML
>
> Parameters
>
>     directory : str, optional
>         Path to directory in which the XML file will be saved. Defaults is the work directory of the model class
>     out : str, optional
>         Name of the XML file. Default is output name of model class
>
> Returns
>
>     filename : str
>         Path to the output XML file

If one is generating a model for use with gtobssim a source list file can also be produced

```python
srcList_path = mymodel.writeSrcList()
```

> fermimodel.model.**writeSrcList(directory=None, out=None)**
>
> Write source list to file
>
> Parameters
>
>     directory : str
>         Path to directory in which the file will be saved.
>     out : str
>         Name of the file
>

> Returns
>
>      filename : str
>         Full path to file
