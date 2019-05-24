# FermiModel

This package is for use with the [Fermitools](https://github.com/fermi-lat/Fermitools-conda/wiki/Quickstart-Guide). The main function of the package is to produce model files for use in analysis with gtlike and simulations with gtobssim. This is essentially an extension of [make4FGLxml.py](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/make4FGLxml.py), written by Tyrel Johnson and the command line syntax is nearly identical to that of [make4FGLxml.py](https://fermi.gsfc.nasa.gov/ssc/data/analysis/user/make4FGLxml.py).

## Installation

Download the package (This will save the package to a directory named `fermimodel` in the user's home directory)

#### BaSH
```bash
# Get the version of the latest release
version=`curl -s https://api.github.com/repos/jaulbric/fermimodel/releases/latest | grep "tag_name" | cut -d ":" -f 2 | tr -d " "-\"-,`
# download and untar the package
wget https://github.com/jaulbric/fermimodel/archive/$version.tar.gz \
&& mkdir ~/fermimodel \
&& tar xvfz $version.tar.gz -C ~/fermimodel --strip-components 1 \
&& rm $version.tar.gz
```

#### csh/tcsh

```bash
# Get the version of the latest release
set version=`curl -s https://api.github.com/repos/jaulbric/fermimodel/releases/latest | grep "tag_name" | cut -d ":" -f 2 | tr -d " "-\"-,`
# download and untar the package
wget https://github.com/jaulbric/fermimodel/archive/$version.tar.gz \
&& mkdir ~/fermimodel \
&& tar xvfz $version.tar.gz -C ~/fermimodel --strip-components 1 \
&& rm $version.tar.gz
```

Then move to the package root directory and install

```bash
cd ~/fermimodel
python setup.py install
```

## Usage

### As Python Module

The following instructions provide a simple example of how to generate a model for use with gtobssim. The roi is centered at the location of the Andromeda galaxy (M31) and will have a radius such that a 10 degree by 10 degree image can be contained within the roi. The Andromeda galaxy will be modeled as an elliptical gaussian with a power law emission spectrum. More information for each command can be found in the [Definitions](#definitions) section. The package can be imported in python as usual.

```python
import fermimodel
```
Models can then be built by creating an instance of the [model class](#model)

```python
mymodel = fermimodel.model(name='mymodel', out='mymodel.xml', frame='galactic', unit='degree', allsky=False, model_type='simulation')
```

The same model instance can be used to generate multiple models. The region of interest can be set either when the model instance is created or by directly calling the function [setROI](#modelsetROI).

```python
mymodel.setROI(roi=(10.8229, 41.2415, 7.071), frame='fk5', unit='degree', allsky=False)
```

Trying to run mission long simulations takes a very long time because the galactic diffuse emission flux is quite large (relatively). It is therefore useful to [mask](#maskFits) the galactic diffuse model so that pixels outside the region of interest have a very low flux.

```python
GDfile, GDflux = fermimodel.maskFits('$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v07.fits',
                                   out='gll_iem_v07_masked.fits', mask_type='radial',
                                   radius=7.071 + 10.,
                                   center=(10.8229, 41.2415),
                                   frame='fk5',
                                   unit='degree',
                                   clobber=True)
```

To add sources from the 4FGL catalog the [loadCatalog](#modelloadcatalog) function is used. For models that will be used by gtobssim it is necessary to have the parameter ```oldNames = True```. This will force all names to begin with an underscore, which is necessary because gtobssim doesn't like source names beginning with '4FGL'. If `model_type` was set to `'simulation'` when the model instance was initialized this will automatically set `oldNames = True`.

```python
mymodel.loadCatalog(catalog = 'gll_psc_v19.fit',
                    GDfile = GDfile,
                    GDname = 'gll_iem_v07',
                    ISOfile = '$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V2_v01.txt',
                    ISOname = 'iso_P8R3_SOURCE_V2_v01',
                    ISOPath = '$(FERMI_DIR)/refdata/fermi/galdiffuse/isotropic_allsky.fits',
                    extDir = '$(FERMI_DIR)/data/pyBurstAnalysisGUI/templates',
                    ExtraRad = 10.,
                    makeRegion = True,
                    wd='.',
                    emin = 1.e2,
                    emax = 5.e5,
                    frame = 'galactic',
                    extSrcRes = 'force-point',
                    apply_mask = False,
                    GDflux=GDflux)
```

Additional Sources can be add individually with [addSource](#modeladdsource). As an example we will add the Andromeda galaxy as a spatially extended source with a power law spectrum.

```python
mymodel.addSource('M31',
                  'PowerLaw',
                  ra=10.8229,
                  dec=41.2415,
                  major_axis=0.46,
                  minor_axis=0.11,
                  position_angle=62.,
                  spatial_function='RadialGauss',
                  emin=1.e2,
                  emax=5.e5,
                  pl_flux_density=4.74295524506e-13,
                  pivot_energy=887.615335377,
                  pl_index=2.52349635581)
```

Sources can also be removed from the model by calling the model class [removeSource](#modelremovesource) method. Since we have included a new definition for the Andromeda galaxy we should remove the old source definition (Note that source names in the model are different then how they appear in the 4FGL catalog. 4FGL catalog sources have names begining with an underscore, spaces replaced by underscores, '.' replaced with 'd', '+' replaced with 'p', and '-' replaced by 'm').

```python
mymodel.removeSource("_4FGL_J0043d2p4114")
```

gtobssim requires a path to an xml file containing the source definitions. We should then export the model to XML using the model class [writeXML](#modelwritexml) method. The type of XML output is determined by the `model_type` parameter when creating an instance of the model class.

```python
xmlfile_path = mymodel.writeXML()
```

gtobssim also requires a list of source names. This is useful because we may want to simulate a subset of the sources in the XML file. A list of sources in the model can be generated with the model class [writeSrcList](#modelwritesrclist) method.

```python
srcList_path = mymodel.writeSrcList()
```

### Via Command Line

The package installs an executable into the user's bin, `fermimodel`. Usage of the executable can be viewed using the `-h` or `--help` flags.

```bash
fermimodel -h
```

```txt
usage: fermimodel [-h] [-v] {buildmodel,maskfits} ...

Command line interface to the fermimodel tools.

positional arguments:
  {buildmodel,maskfits}
    buildmodel          Creates an xml model from the 4FGL catalog (FITS or
                        xml version) for a specific ROI. Coordinates of the
                        ROI center are taken from an input event file or
                        defined on the command line. For likelihood models
                        sources with free parameters within the original
                        extraction radius are chosen based on nearness to
                        center, significance, and variability. For simulation
                        models all sources within radLim + ExtraRad are
                        included in the model. Optionally, a file can be input
                        on the command line by using the '@' prefix.
                        Parameters inside this file should be one per line and
                        are read in order.
    maskfits            Mask fits image. The data HDU of the FITS file should
                        be a 2D image or a 3D cube with the 3rd axis being
                        energy (0 axis numpy indexing).

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

There are two commands for the executable: `buildmodel` and `maskfits`. For more information on these commands use the `-h` or `--help` flags.

```bash
fermimodel buildmodel -h
```

```txt
usage: fermimodel buildmodel [-h] [-n NAME] [-t {likelihood,simulation}]
                             [-o OUTPUTXML] [-G GALFILE] [-g GALNAME]
                             [-I ISOFILE] [-i ISONAME] [-e EXTDIR]
                             [-roi C1, C2, RADIUS] [-xc HORIZONTAL_CENTER]
                             [-yc VERTICAL_CENTER] [-icrs | -gal | -fk5]
                             [-u UNIT] [-r RADLIM] [-ER EXTRARAD] [-p] [-E2C]
                             [-m] [-wd WRITEDIR] [-ON] [-AS] [-N] [-R MAXRAD]
                             [-s SIGFREE] [-GIF] [--emin EMIN] [--emax EMAX]
                             [--isopath ISOPATH]
                             [--extSrcRes {force-point,skip,raise}] [-am]
                             [-gf GALFLUX] [-if ISOFLUX]
                             catalog [ev]

positional arguments:
  catalog               Catalog file to use, can be FITS or xml.
  ev                    Event file with ROI information in header.

optional arguments:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  Name of the model. This will be used for output
                        filenames if no other name is specified. Default is
                        'mymodel'.
  -t {likelihood,simulation}, --model-type {likelihood,simulation}
                        Type of model to generate. For likelihood analysis
                        this should be 'likelihood'. If the model is intended
                        for use with gtobssim this should be 'simulation'.
  -o OUTPUTXML, --outputxml OUTPUTXML
                        Name of output xml file, is set to overwrite files of
                        same name.
  -G GALFILE, --galfile GALFILE
                        Path to Galactic diffuse model. Default is $(FERMI_DIR
                        )/refdata/fermi/galdiffuse/gll_iem_v07.fits.
  -g GALNAME, --galname GALNAME
                        Name of Galactic diffuse component in output model.
                        Default is gll_iem_v07.
  -I ISOFILE, --isofile ISOFILE
                        Path to isotropic diffuse template for output model,
                        will default to P8R3 SOURCE class model.
  -i ISONAME, --isoname ISONAME
                        Name of isotropic diffuse component in output model,
                        default is for P8R3 SOURCE class.
  -e EXTDIR, --extDir EXTDIR
                        Path to directory with LAT extended source templates,
                        will default to STs default.
  -roi (C1, C2, RADIUS), --Region (C1, C2, RADIUS)
                        Region of Interest. This will override values taken
                        from event file if one is supplied.
  -xc HORIZONTAL_CENTER, --horizontal-center HORIZONTAL_CENTER
                        Horizontal coordinate of ROI center.
  -yc VERTICAL_CENTER, --vertical-center VERTICAL_CENTER
                        Vertical coordinate of ROI center.
  -icrs, --celestial    Flag sets coordinates of roi to RA, DEC.
  -gal, --galactic      Flag sets coordinates of roi to GLON, GLAT.
  -fk5, --J2000         Flag sets coordinates of roi to RAJ2000, DECJ2000.
  -u UNIT, --unit UNIT  Units for coordinates of roi using astropy units
                        convention. Default is 'degree'.
  -r RADLIM, --radLim RADLIM
                        Radius, in degrees, from ROI center beyond which all
                        source parameters should be fixed, will default to
                        selection radius. If --obssim flag is set this is used
                        as the radius of the ROI.
  -ER EXTRARAD, --ExtraRad EXTRARAD
                        Radius beyond event file ROI out to which sources will
                        be included in the model with all parameters fixed,
                        default is 10, good for analyses starting around a few
                        hundred MeV, can be decreased for high energy only
                        fits. If --obssim flag is set ExtraRad is added to
                        radLim when including sources.
  -p, --psForce         Flag to cast extended sources as point sources.
                        Default is False.
  -E2C, --E2CAT         Flag to use catalog names for extended sources.
                        Default is False.
  -m, --makeRegion      Flag to create ds9 region file as well as the xml
                        model. Default is False.
  -wd WRITEDIR, --writeDir WRITEDIR
                        Directory in which to write output files. Default is
                        the current directory.
  -ON, --oldNames       Flag to use the make2FGLxml style naming convention,
                        underscore before name and no spaces. Default is False
  -AS, --allsky         Generate model with all sources from the catalog.
                        Overides region of interest parameters.

Likelihood:
  These parameters only affect models intended for use in likelihood
  analysis.

  -N, --normsonly       Flag to only let the normalizations of parameters be
                        free, default is False.
  -R MAXRAD, --maxRad MAXRAD
                        Absolute maximum radius, in degrees, from ROI center
                        beyond which all source parameters should be fixed,
                        even variable sources will not be freed beyond this
                        radius, defaults to radLim value.
  -s SIGFREE, --sigFree SIGFREE
                        Average significance below which all source parameters
                        are fixed, defaults to 5. Note, if using the 3FGL
                        catalog xml file as input, this is actually a cut on
                        TS, so adjust accordingly.
  -GIF, --GIndexFree    Flag to use a power-law modification to the Galactic
                        diffuse model, spectrum and have the index be free.
                        Default is False

Simulation:
  These parameters only affect models intended for use with gtobssim.

  --emin EMIN           Minimum energy for integrated flux.
  --emax EMAX           Maximum energy for integrated flux.
  --isopath ISOPATH     Path to isotropic spatial file (a FITS image
                        containing pixels all set to 1.). Default is $(FERMI_D
                        IR)/refdata/fermi/galdiffuse/isotropic_allsky.fits
  --extSrcRes {force-point,skip,raise}
                        Resolution method for adding an extended source when
                        no extended template is found.
  -am, --apply-mask     Apply a region of interest mask to the diffuse models.
  -gf GALFLUX, --galflux GALFLUX
                        Integrated flux (photons/cm^2/s) to use for galactic
                        diffuse emmision model. Default is 0.00084631675.
  -if ISOFLUX, --isoflux ISOFLUX
                        Integrated flux (photons/cm^2/s) to use for isotropic
                        diffuse emission. Default is 0. (If 0. gtobssim will
                        integrate the spectrum)
```

The maskfits functionallity is also accessable via the command line by invoking the fermimodel executable with the command `maskfits`.

```bash
fermimodel maskfits -h
```

```txt
usage: fermimodel maskfits [-h] [-o OUTPUT] [-m {radial,square}]
                           [-ih IMAGE_HDU] [-cl] [-fk5 | -icrs | -gal | -pix]
                           [-u UNIT] [-roi C1, C2, RADIUS] [-r RADIUS]
                           [-r2 RADIUS2] [-a ANGLE] [-xc HORIZONTAL_CENTER]
                           [-yc VERTICAL_CENTER] [-xmin HORIZONTAL_MIN]
                           [-xmax HORIZONTAL_MAX] [-ymin VERTICAL_MIN]
                           [-ymax VERTICAL_MAX]
                           input

positional arguments:
  input                 Fits file containing the image to which a mask will be
                        applied.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Optional filename to save the masked data. Default is
                        masked.fits
  -m {radial,square}, --mask {radial,square}
                        Shape of the mask.
  -ih IMAGE_HDU, --image_hdu IMAGE_HDU
                        HDU containing the image to be masked. If none is
                        input it is assumed to live in the PRIMARY HDU.
  -cl, --clobber        Flag to overwrite a file of the same name. Default is
                        false.
  -fk5, --J2000         Flag sets coordinates of mask center to RAJ2000,
                        DECJ2000
  -icrs, --celestial    Flag sets coordinates of mask center to RA, DEC
  -gal, --galactic      Flag sets coordinates of mask center to GLON, GLAT
  -pix, --pixel         Flag sets coordinates of mask center to PIXEL1, PIXEL2
  -u UNIT, --unit UNIT  Units of mask coordinates. Default is degrees.

radial:
  Parameters for radial mask

  -roi (C1, C2, RADIUS), --Region (C1, C2, RADIUS)
                        Region of Interest
  -r RADIUS, --radius RADIUS
                        Radius of the mask.
  -r2 RADIUS2, --radius2 RADIUS2
                        Second radius of ellipse. Not yet implemented
  -a ANGLE, --angle ANGLE
                        Angle of radius with respect to the horizontal axis.
                        not yet implemented.
  -xc HORIZONTAL_CENTER, --horizontal-center HORIZONTAL_CENTER
                        Horizontal coordinate of mask center.
  -yc VERTICAL_CENTER, --vertical-center VERTICAL_CENTER
                        Vertical coordinate of mask center.

square:
  Parameters for square mask

  -xmin HORIZONTAL_MIN, --horizontal-min HORIZONTAL_MIN
                        Minimum horizontal coordinate value.
  -xmax HORIZONTAL_MAX, --horizontal-max HORIZONTAL_MAX
                        Maximum horizontal coorindate value.
  -ymin VERTICAL_MIN, --vertical-min VERTICAL_MIN
                        Minimum vertical coordinate value.
  -ymax VERTICAL_MAX, --vertical-max VERTICAL_MAX
                        Maximum vertical coordinate value.
```

## Definitions

### model

>fermimodel.**model(name='mymodel', eventsfile=None, catalog=None, out=None, roi=None, frame='fk5', unit='degree', allsky=False, model_type='likelihood')**
>
> Model class containing methods and attributes for building models for use with Fermitools
>
> #### Parameters
>
> <dl>
>   <dt>name : str (Optional)</dt>
>   <dd>Name of the model.</dd>
>   <dd>Default is 'mymodel'.</dd>
>   <dt>eventsfile : str (Optional)</dt>
>   <dd>Path to file containing events data. This is only used for region of interest information.</dd>
>   <dt>catalog : str (Optional)</dt>
>   <dd>Path to source catalog.</dd>
>   <dt>out : str (Optional)</dt>
>   <dd>Name of the output model files.</dd>
>   <dt>roi : tuple (Optional)</dt>
>   <dd>Tuple defining ROI (horizontal center, vertical center, radius).</dd>
>   <dt>frame : str (Optional)</dt>
>   <dd>Coordinate frame to use for region of interest argument, e.g. galactic, icrs. Any astropy frame is allowed.</dd>
>   <dd>Default is 'fk5'.</dd>
>   <dt>unit : str (Optional)</dt>
>   <dd>Units for region of interest argument. Any astropy unit is allowed.</dd>
>   <dd>Default is 'degree'</dd>
>   <dt>allsky : bool (Optional)</dt>
>   <dd>Flag to set region of interest to the entire sky. By default this sets the region of interest center to the galactic center (l=0, b=0).</dd>
>   <dd>Default is False.</dd>
>   <dt>model_type : str (Optional)</dt>
>   <dd>Model type. Choices are 'simulation' or 'likelihood'. Input parameters are identical for both model types but the output model will vary depending on whether the model will be used by gtlike or gtobssim.
    <dd>Default is 'likelihood'.</dd>
> </dl>

### model.setROI

> fermimodel.model.**setROI(roi=None, frame='fk5', unit='degree', allsky=False)**
>
> Set the Region of Interest for the model.
>
> #### Parameters
>
> <dl>
>   <dt>roi : tuple</dt>
>   <dd>(C1, C2, RADIUS)</dd>
>   <dt>frame : string</dt>
>   <dd>Coordinate frame for roi center coordinates C1 and C2</dd>
>   <dt>unit : str</dt>
>   <dd>Units of roi center coordinates, e.g. 'degree'</dd>
>   <dt>allsky : bool</dt>
>   <dd>Flag to set roi to GLAT = 0, GLON = 0, RADIUS = 180.</dd>  
>
> #### Returns

### maskFits

> fermimodel.**maskFits(fitsfile, out='maskedimage.fits', img_hdu=None, mask_type=None, radius=180., radius2=None, angle=0., center=(0., 0.), extent=[-180., 180., -90., 90.], frame='galactic', unit='degree', clobber=False, float_min=1.17549e-38)**
>
> Mask the fits image
>
> #### Parameters
>
> <dl>
>   <dt>fitsfile : str</dt>
>   <dd>Path to the FITS file containing the the image which will be masked.</dd>
>   <dt>out : str (Optional)</dt>
>   <dd>Name of the output FITS file for which the mask wiil be applied</dd>
>   <dt>img_hud : int or float (Optional)</dt>
>   <dd>Name or integer for the FITS hdu containing the image.</dd>
>   <dd>Default is 'Primary'.</dd>
>   <dt>mask_type : str</dt>
>   <dd>The geometry of the mask to be applied. Choices are 'radial' or 'square'.</dd>
>   <dt>radius : float (Optional)</dt>
>   <dd>Radius of the mask if mask_type is radial.</dd>
>   <dd>Default is 180.</dd>
>   <dt>radius2 : float (Optional)</dt>
>   <dd>Second radius of the mask if the mask is not symmetric.</dd>
>   <dd>Default is to use a symmetric mask.</dd>
>   <dt>angle : float (Optional)</dt>
>   <dd>Rotation angle of the ellipse.</dd>
>   <dd>Default is 0.</dd>
>   <dt>center : tuple (Optional)</dt>
>   <dd>Center coordinates (C1, C2) of the radial mask.</dd>
>   <dd>Default is (0., 0.).</dd>
>   <dt>extent : list (Optional)</dt>
>   <dd>[xmin, xmax, ymin, ymax] extent of the square mask.</dd>
>   <dd>Default is [-180., 180., -90., 90.].</dd>
>   <dt>frame : str (Optional)</dt>
>   <dd>Coordinate frame to use with mask coordinates. Choices are 'galactic', 'icrs', 'fk5', 'pixel'.</dd>
>   <dd>Default is 'galactic'.</dd>
>   <dt>unit : str (Optional)</dt>
>   <dd>Units of coordinates.</dd>
>   <dd>Default is 'degree'.</dd>
>   <dt>clobber : bool</dt>
>   <dd>Flag to overwrite a file of the same name if it exists.</dd>
>   <dd>Default is False.</dd>
>   <dt>float_min : float</dt>
>   <dd>Minimum float value to use for pixels in the image after masking. gtobssim doesn't like pixel values <= 0.</dd>
>   <dd>Default is 1.17549e-38.</dd>
> </dl>
>
> #### Returns
>
> <dl>
>   <dt>out : str or tuple</dt>
>   <dd>If the input fits image is 2D the output is the full path to the masked fits image. If the input fits image is 3D the output is a tuple whose first entry is the full path to the masked fits image and whose second entry is the integrated flux of the masked fits image.</dd>
> </dl>

### model.loadCatalog

>fermimodel.model.**loadCatalog(catalog=None, GDfile=None, GDname='GalacticDiffuse', ISOfile=None, ISOname='IsotropicDiffuse', ISOpath="$(FERMI_DIR)/refdata/fermi/galdiffuse/isotropic_allsky.fits", normsOnly=False, extDir='', radLim=-1, maxRad=None, ExtraRad=10., sigFree=5., varFree=True, psForce=False, E2CAT=False, makeRegion=False, GIndexFree=False, wd=None, oldNames=False, emin=1.e2, emax=5.e5, frame='fk5', extSrcRes='force-point', apply_mask=False, GDflux=0.00084631675, ISOflux=0.)**
>
>Include sources in the catalog to the model. Optionaly include a galactic and isotropic diffuse models.
>
> #### Parameters
>
> <dl>
>   <dt>catalog : str</dt>
>   <dd>Path to source catalog</dd>
>   <dt>GDfile : str (Optional)</dt>
>   <dd>Path to galactic diffuse emission file</dd>
>   <dt>GDname : str (Optional)</dt>
>   <dd>Name of galactic diffuse emission model</dd>
>   <dt>ISOfile : str (Optional)</dt>
>   <dd>Path to isotropic diffuse emission spectrum file</dd>
>   <dt>ISOname : str (Optional)</dt>
>   <dd>Name of isotropic diffuse emission file</dd>
>   <dt>ISOpath : str (Optional)</dt>
>   <dd>Path to isotropic diffuse spatial file</dd>
>   <dt>normsOnly : bool (Optional)</dt>
>   <dd>Flag to only set normalization parameters free</dd>
>   <dt>extDir : str (Optional)</dt>
>   <dd>Directory with extended source templates</dd>
>   <dt>radLim : float (Optional)</dt>
>   <dd>Radius in degrees from center of ROI beyond which source parameters are fixed</dd>
>   <dt>maxRad : float (Optional)</dt>
>   <dd>Absolute maximum radius beyond which sources are fixed. This may be necessary when doing binned analysis and a variable source beyond radLim would be set free but this source is beyond the boundaries of the square region used for the binned likelihood</dd>
>   <dt>ExtraRad : float (Optional)</dt>
>   <dd>Radius beyond ROI radius out to which sources will be included with fixed parameters. Default of 10 degrees is good for analysis starting around 100 MeV, but for higher energy fits this can be decreased.</dd>
>   <dt>sigFree : float (Optional)</dt>
>   <dd>Significance below which source parameters are fixed, even if within radLim.</dd>
>   <dt>varFree : float (Optional)</dt>
>   <dd>Variability index above which source parameters are free. If beyond radLim and/or below sigFree only the normalization parameters is set free. Currently not implemented for building from xml catalog.</dd>
>   <dt>psForce : bool (Optional)</dt>
>   <dd>Flag to force exentended sources to be point sources</dd>
>   <dt>E2CAT : bool (Optional)</dt>
>   <dd>Flag to force use of catalog names for extended sources (only matters is using catalog FITS file).</dd>
>   <dt>makeRegion : bool (Optional)</dt>
>   <dd>Flag to also generate ds9 region file.</dd>
>   <dt>GIndexFree : bool (Optional)</dt>
>   <dd>The galactic diffuse is given a power-law spectral shape but by default the index is frozen. Setting this flag to True allows that to be free for additional freedom in diffuse fit.</dd>
>   <dt>wd : str (Optional)</dt>
>   <dd>Path to directory in which output files will be saved. If an absolute path for the output file names is given it will override this argument.</dd>
>   <dt>oldNames : bool (Optional)</dt>
>   <dd>Sets use of old naming convention. Underscore before name and no spaces. Default is False for likelihood models. This is automatically set to true for simulation models.</dd>
>   <dt>emin : float (Optional)</dt>
>   <dd>Minimum energy in MeV for source simulation. This should match simulation criteria.</dd>
>   <dd>Default is 100 MeV.</dd>
>   <dt>emax : float (Optional)</dt>
>   <dd>Maximum energy in MeV for source simulation. This should match simulation criteria.</dd>
>   <dd>Default is 5.e5 MeV.</dd>
>   <dt>GDflux : float (Optional)</dt>
>   <dd>Integrated flux to use for the galactic diffuse emission model in photons/cm^2/s.</dd>
>   <dd>Default is 0.00084631675.</dd>
>   <dt>ISOflux : float (Optional)</dt>
>   <dd>Integrated flux to use for the isotropic diffuse emission model in photons/cm^2/s. If 0.0 flux is calculated on the fly by the simulator.</dd>
>   <dd>Default is 0.0.</dd>
> </dl>
>
> #### Returns

### model.addSource

> fermimodel.model.**addSource(name, spectrum_type, ra=None, dec=None, glon=None, glat=None, major_axis=None, minor_axis=None, position_angle=None, spatial_function=None, extended_template=None, emin=1.e2, emax=5.e5, frame='galactic', specFile=None, \*\*spectrumargs)**
>
> Add a single source to the model.
>
> #### Parameters
>
>  <dl>
>    <dt>name : str</dt>
>    <dd>Name of the source to add.</dd>
>    <dt>spectrum_type : str</dt>
>    <dd>Spectral type of the source. Options: Monochromatic, PowerLaw, BrokenPowerLaw, LogParabola, PLSuperExpCutoff, MapCube, FileSpectrum</dd>
>    <dt>ra : float (Optional)</dt>
>    <dd>Right Ascension of source in degrees</dd>
>    <dt>dec : float (Optional)</dt>
>    <dd>Decilination of source in degrees</dd>
>    <dt>glon : float (Optional)</dt>
>    <dd>Galactic longitude of source in degrees</dd>
>    <dt>glat : float (Optional)</dt>
>    <dd>Galactic latitude of source in degrees</dd>
>    <dt>major_axis : float (Optional)</dt>
>    <dd>The semi-major axis of the error ellipse at 95% confidence, in degrees.</dd>
>    <dt>minor_axis : float (Optional)</dt>
>    <dd>The semi-minor axis of the error ellipse at 95% confidence, in degrees.</dd>
>    <dt>position_angle : float (Optional)</dt>
>    <dd>The position angle of the 95%-confidence semi-major axis, from celestial North, positive toward increasing R.A. (eastward), in degrees.</dd>
>    <dt>spatial_function : str (Optional)</dt>
>    <dd>Spatial function describing source extension. Options: None, SpatialMap, RadialDisk, RadialGauss, Isotropic</dd>
>    <dt>extended_template : str (Optional)</dt>
>    <dd>Full path to extended template.</dd>
>    <dt>emin : float (Optional)</dt>
>    <dd>Minimum energy for integrated flux.</dd>
>    <dd>Default is 1.e2 MeV</dd>
>    <dt>emax : float (Optional)</dt>
>    <dd>Maximum energy for integrated flux.</dd>
>    <dd>Default is 5.e5 MeV</dd>
>    <dt>frame : str (Optional)</dt>
>    <dd>Coordinate frame to use for source directions.</dd>
>    <dt>specFile : str Optional</dt>
>    <dd>Name (and path) to ASCII file containing spectrum. The first column contains the energies, and second column contains the differential flux at those energies. Energy scale should be in MeV.</dd>
> </dl>
>
> #### \*\*spectrumargs
>
> <dl>
>   <dt>pl_flux_density : float (Optional)</dt>
>   <dd>The differential flux at pivot_energy for the PowerLaw fit, in photons/cm2/MeV/s.</dd>
>   <dt>lp_flux_density : float (Optional)</dt>
>   <dd>The differential flux at pivot_energy for the LogParabola fit, in photons/cm2/MeV/s.</dd>
>   <dt>plec_flux_density : float (Optional)</dt>
>   <dd>The differential flux at pivot_energy for the PLSuperExpCutoff fit, in photons/cm^2/MeV/s.</dd>
>   <dt>pivot_energy : float (Optional)</dt>
>   <dd>The energy, in MeV, at which the error in the differential photon flux is minimal (i.e., the decorrelation energy for the power-law fit).</dd>
>   <dt>pl_index : float (Optional)</dt>
>   <dd>The photon index for the PowerLaw fit.</dd>
>   <dt>lp_index : float (Optional)</dt>
>   <dd>Photon Index at Pivot Energy for LogParabola fit.</dd>
>   <dt>lp_beta : float (Optional)</dt>
>   <dd>The curvature parameter (beta) for LogParabola spectrum types.</dd>
>   <dt>plec_index : float (Optional)</dt>
>   <dd>The low energy photon index for PLSuperExpCutoff fit.</dd>
>   <dt>plec_expfactor : float (Optional)</dt>
>   <dd>The exponential factor for PLSuperExpCutoff fit.</dd>
>   <dt>plec_exp_index : float (Optional)</dt>
>   <dd>The exponential index for PLSuperExpCutoff fit.</dd>
>   <dt>flux : float (Optional)</dt>
>   <dd>Integrated flux in photons/cm^2/s.</dd>
>   <dd>Default is to calculate the flux from spectral parameters.</dd>
> </dl>
>
> #### Returns

### model.removeSource

> fermimodel.model.**.removeSource(self, \*name)**
>
> #### Parameters
>
> <dl>
>   <dt>name : str</dt>
>   <dd>Variable number of source names to remove.</dd>
> </dl>
>
> #### Returns

### model.writeXML

> fermimodel.model.**writeXML(directory=None, out=None)**
>
> Write the model to XML.
>
> #### Parameters
>
> <dl>
>   <dt>directory : str (Optional)</dt>
>   <dd>Path to directory in which the XML file will be saved.</dd>
>   <dd>Defaults is the work directory of the model class.</dd>
>   <dt>out : str (Optional)</dt>
>   <dd>Name of the XML file.</dd>
>   <dd>Default is output name of model class.</dd>
> </dl>
>
> #### Returns
>
> <dl>
>   <dt>filename : str</dt>
>   <dd>Full path to output XML file.</dd>
> </dl>

### model.writeSrcList

> fermimodel.model.**writeSrcList(directory=None, out=None)**
>
> Write source list to file.
>
> #### Parameters
>
> <dl>
>   <dt>directory : str (Optional)</dt>
>   <dd>Path to directory in which the file will be saved.</dd>
>   <dd>Default is the work directory of the model class.</dd>
>   <dt>out : str (Optional)</dt>
>   <dd>Name of the file</dd>
>   <dd>Default is the name of the model.</dd>
> </dl>
>
> #### Returns
>
> <dl>
>    <dt>filename : str</dt>
>    <dd>Full path to output file.</dd>
> </dl>