import argparse
import os

from fermimodel import Model

version = "02r04"

def cli():
    helpString="Creates an xml model from the 3FGL catalog (FITS or xml version) for a specific ROI,\
            coordinates of the ROI center are taken from an input event file,\
            the radius for including sources is 10 degrees beyond the extraction radius used in the event file,\
            sources with free parameters within the original extraction radius are chosen based on nearness to center, significance, and variability."
    parser=argparse.ArgumentParser(description=helpString)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s {0}".format(version))

    likelihood_group = parser.add_argument_group("Likelihood", "These parameters only affect models intended for use in likelihood analysis.")
    simulation_group = parser.add_argument_group("Simulation", "These parameters only affect models intended for use with gtobssim.")

    frame_group = parser.add_mutually_exclusive_group()

    parser.add_argument("catalog", type=str, help="Catalog file to use, can be FITS or xml.")
    parser.add_argument("ev", type=str, nargs="?", help="Event file with ROI information in header.")
    parser.add_argument("-n", "--name", type=str, default="mymodel", help="Name of the model. This will be used for output filenames if no other name is specified. Default is 'mymodel'.")
    parser.add_argument("-t", "--model-type", type=str, choices=['likelihood', 'simulation'], default='likelihood', help="Type of model to generate. For likelihood analysis this should be 'likelihood'. If the model is intended for use with gtobssim this should be 'simulation'.")
    parser.add_argument("-o", "--outputxml", type=str, default='mymodel.xml', help="Name of output xml file, is set to overwrite files of same name.")
    parser.add_argument("-G", "--galfile", type=str, default='$(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v07.fits', help="Path to Galactic diffuse model. Default is $(FERMI_DIR)/refdata/fermi/galdiffuse/gll_iem_v07.fits.")
    parser.add_argument("-g", "--galname", type=str, default='gll_iem_v07', help="Name of Galactic diffuse component in output model. Default is gll_iem_v07.")
    parser.add_argument("-I", "--isofile", type=str, default='$(FERMI_DIR)/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V2_v1.txt', help="Path to isotropic diffuse template for output model, will default to P8R3 SOURCE class model.")
    parser.add_argument("-i", "--isoname",type=str, default='iso_P8R3_SOURCE_V2_v1', help="Name of isotropic diffuse component in output model, default is for P8R3 SOURCE class.")
    likelihood_group.add_argument("-N", "--normsonly", action="store_true", help="Flag to only let the normalizations of parameters be free, default is False.")
    parser.add_argument("-e", "--extDir", type=str, default='', help="Path to directory with LAT extended source templates, will default to STs default.")#need to figure out what that is
    parser.add_argument("-roi", "--Region", type=float, nargs=3, metavar=("C1", "C2", "RADIUS"), help="Region of Interest. This will override values taken from event file if one is supplied.")
    
    frame_group.add_argument("-icrs", "--celestial", action="store_const", const="icrs", help="Flag sets coordinates of roi to RA, DEC.")
    frame_group.add_argument("-gal", "--galactic", action="store_const", const="galactic", help="Flag sets coordinates of roi to GLON, GLAT.")
    frame_group.add_argument("-fk5", "--J2000", action="store_const", const="fk5", help="Flag sets coordinates of roi to RAJ2000, DECJ2000.")

    parser.add_argument("-u", "--unit", type=str, default='deg', help="Units for coordinates of roi using astropy units convention. Default is 'degree'.")
    parser.add_argument("-r", "--radLim", type=float, default=-1., help="Radius, in degrees, from ROI center beyond which all source parameters should be fixed, will default to selection radius. If --obssim flag is set this is used as the radius of the ROI.")
    likelihood_group.add_argument("-R", "--maxRad", type=float, default=None, help="Absolute maximum radius, in degrees, from ROI center beyond which all source parameters should be fixed, even variable sources will not be freed beyond this radius, defaults to radLim value.")
    parser.add_argument("-ER", "--ExtraRad", type=float, default=10., help="Radius beyond event file ROI out to which sources will be included in the model with all parameters fixed, default is 10, good for analyses starting around a few hundred MeV, can be decreased for high energy only fits. If --obssim flag is set ExtraRad is added to radLim when including sources.")
    likelihood_group.add_argument("-s", "--sigFree", type=float, default=5., help="Average significance below which all source parameters are fixed, defaults to 5.  Note, if using the 3FGL catalog xml file as input, this is actually a cut on TS, so adjust accordingly.")
    parser.add_argument("-p","--psForce", action="store_true", help="Flag to cast extended sources as point sources. Default is False.")
    parser.add_argument("-E2C", "--E2CAT", action="store_true", help="Flag to use catalog names for extended sources. Default is False.")
    # parser.add_argument("-m","--makeRegion",type=mybool,default=True,help="Flag to create ds9 region file as well as the xml model, default is True.",choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    parser.add_argument("-m", "--makeRegion", action="store_true", help="Flag to create ds9 region file as well as the xml model. Default is False.")
    likelihood_group.add_argument("-GIF", "--GIndexFree", action="store_true", help="Flag to use a power-law modification to the Galactic diffuse model, spectrum and have the index be free. Default is False")
    #parser.add_argument("-ED","--edisp",type=mybool,default=False,help="Flag to turn on energy dispersion for free point and extended sources, never for diffuse backgrounds, default is False.",nargs="?",const=True,choices=['True','False','T','F','t','f','TRUE','FALSE','true','false',1,0])
    parser.add_argument("-wd", "--writeDir", type=str, default=os.getcwd(), help="Directory in which to write output files. Default is the current directory.")
    parser.add_argument("-ON", "--oldNames", action="store_true", help="Flag to use the make2FGLxml style naming convention, underscore before name and no spaces. Default is False")
    parser.add_argument("-AS", "--allsky", action="store_true", help="Generate model with all sources from the catalog. Overides region of interest parameters.")

    # simulation parameters
    simulation_group.add_argument("--emin", type=float, default=1.e2, help="Minimum energy for integrated flux.")
    simulation_group.add_argument("--emax", type=float, default=5.e5, help="Maximum energy for integrated flux.")
    simulation_group.add_argument("--extSrcRes", type=str, choices=["force-point", "skip", "raise"], default="force-point", help="Resolution method for adding an extended source when no extended template is found.")
    simulation_group.add_argument("-am", "--apply_mask", action="store_true", help="Apply a region of interest mask to the diffuse models.")
    simulation_group.add_argument("-gf", "--galflux", type=float, default=0.00084631675, help="Integrated flux (photons/cm^2/s) to use for galactic diffuse emmision model. Default is 0.00084631675.")
    simulation_group.add_argument("-if", "--isoflux", type=float, default=0., help="Integrated flux (photons/cm^2/s) to use for isotropic diffuse emission. Default is 0. (If 0. gtobssim will integrate the spectrum)")
    
    args=parser.parse_args()
    args.varFree=False#remove this later

    print "This is make4FGLxml version {0}".format(version)
    print "The default diffuse model files and names are for pass 8 and 4FGL and assume you have v11r5p3 of the Fermi Science Tools or higher."

    if args.model_type == 'simulation':
        args.oldNames = True

    if args.celestial is not None:
        frame = args.celestial
    elif args.galactic is not None:
        frame = args.galactic
    elif args.J2000 is not None:
        frame = args.J2000
    else:
        frame = 'fk5'
    
    model = Model.model(name='mymodel', eventsfile=args.ev, catalog=args.catalog, out=args.outputxml, roi=args.Region, frame=frame, unit=args.unit, allsky=args.allsky, model_type=args.model_type)
    model.loadCatalog(GDfile=args.galfile, GDname=args.galname, ISOfile=args.isofile, ISOname=args.isoname, normsOnly=args.normsonly, extDir=args.extDir, radLim=args.radLim, maxRad=args.maxRad, ExtraRad=args.ExtraRad, sigFree=args.sigFree, varFree=args.varFree, psForce=args.psForce, E2CAT=args.E2CAT, makeRegion=args.makeRegion, GIndexFree=args.GIndexFree, wd=args.writeDir, oldNames=args.oldNames, emin=args.emin, emax=args.emax, frame=frame, extSrcRes=args.extSrcRes, apply_mask=args.apply_mask, GDflux=args.galflux, ISOflux=args.isoflux)

    xmlpath = model.writeXML()
    print "Model XML has been written to {}".format(xmlpath)

    if args.model_type == 'simulation':
        srcpath = model.writeSrcList()
        "Source list has been written to {0}".format(srcpath)
    
if __name__=='__main__': cli()