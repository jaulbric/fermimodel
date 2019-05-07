import astropy.io.fits as pyfits
from astropy.coordinates import SkyCoord
from xml.dom import minidom
# from xml.dom.minidom import parseString as pS
import numpy as np
import os

import Tools
import BuildRegion
import LikelihoodSpectra
import SimulationSources
import maskFits

class likelihood:
    def __init__(self, **params):
        for key, val in params.items():
            setattr(self, key, val)

        self.varValue = 72.44

        extension = os.path.splitext(self.srcs)[-1]
        if extension == '.xml':
            self.model, self.Sources = self.xml()
        elif extension in ['.fits', '.fit']:
            self.model, self.Sources = self.fits()
        else:
            raise IOError("{0} is not a compatible catalog.".format(self.srcs))

    def xml(self):
        inputXml = minidom.parse(self.srcs)
        outputXml = minidom.getDOMImplementation().createDocument(None, 'source_library', None)
        outputXml.documentElement.setAttribute('title', 'source_library')
        catalog = inputXml.getElementsByTagName('source')
        Sources = {}
        prSrcNum = 0
        extSrcNum = 0

        for src in catalog:
            if src.getAttribute('type') == 'PointSource':
                for p in src.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
                    if p.getAttribute('name') == 'RA':
                        srcRA = float(p.getAttribute('value'))
                    elif p.getAttribute('name') == 'DEC':
                        srcDEC = float(p.getAttribute('value'))
            else:
                try:
                    srcDEC = float(src.getAttribute('DEC'))
                    srcRA = float(src.getAttribute('RA'))
                except:
                    for p in src.getElementsByTagName('spatialModel')[0].getElementsByTagName('parameter'):
                        if p.getAttribute('name') == 'RA':
                            srcRA=float(p.getAttribute('value'))
                        if p.getAttribute('name') == 'DEC':
                            srcDEC=float(p.getAttribute('value'))

            dist = Tools.angsep(self.roi[0], self.roi[1], srcRA, srcDEC) #check that source is within ROI radius + 10 degress of ROI center

            c = SkyCoord(ra=srcRA, dec=srcDEC, frame='fk5', unit='degree')
            srcGL = c.galactic.l.degree
            srcGB = c.galactic.b.degree
            
            if (srcRA == self.roi[0]) and (srcDEC == self.roi[1]):
                dist=0.0
            
            if dist <= self.roi[2] + self.ER:
                spec = src.getElementsByTagName('spectrum')
                specType = spec[0].getAttribute('type')
                specPars = spec[0].getElementsByTagName('parameter')
                Ext = (True if (src.getAttribute('type') == 'DiffuseSource' and not self.psF) else False)
                sname = src.getAttribute('name')
                #comment out the stuff below for now...likely add something back for official 4FGL release
                #fixAll=(True if str(sname) in ['3FGL J0534.5+2201i','3FGL J0833.1-4511e','3FGL J1514.0-5915e','3FGL J2021.0+4031e','3FGL J2028.6+4110e'] else False)#account for sources held fixed in 3FGL analysis
                
                if self.oldNames:#if you want the same naming convention as in make1FGLxml.py and make2FGLxml.py, e.g., preceeded by an underscore and no spaces
                    sn = '_' + sname.replace(' ', '')
                varIdx = float(src.getAttribute('Variability_Index'))
                Sources[sname] = {'ra':srcRA, 'dec':srcDEC, 'glon':srcGL, 'glat':srcGB, 'E':Ext, 'stype':str(specType)}
                specOut = outputXml.createElement('spectrum')
                
                if str(specType) == 'PLSuperExpCutoff2':
                    specOut.setAttribute('type', 'PLSuperExpCutoff')
                else:
                    specOut.setAttribute('type', specType)
                
                spatialOut = outputXml.createElement('spatialModel')
                srcOut = outputXml.createElement('source')
                srcOut.setAttribute('name', sname)
                srcOut.setAttribute('ROI_Center_Distance',"%.2f"%dist)
                
                if dist >= self.roi[2] or dist >= self.maxRad:
                    Sources[sname]['free'] = False
                    #specOut.setAttribute('apply_edisp','false')#source is fixed, so never apply edisp
                    for p in specPars:
                        specOut.appendChild(Tools.parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                elif dist > self.radLim:
                    if self.var and varIdx >= self.varValue:
                        Sources[sname]['free'] = True
                        #specOut.setAttribute('apply_edisp',ed)
                        for p in specPars:
                            freeFlag=("1" if p.getAttribute('name') == spec[0].getAttribute('normPar') else "0")
                            specOut.appendChild(Tools.parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                    else:
                        Sources[sname]['free'] = False
                        #specOut.setAttribute('apply_edisp','false')
                        for p in specPars:
                            specOut.appendChild(Tools.parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                elif float(src.getAttribute('TS_value')) >= self.sig:
                    Sources[sname]['free'] = True
                    #specOut.setAttribute('apply_edisp',ed)
                    for p in specPars:
                        freeFlag=("1" if p.getAttribute('name') == spec[0].getAttribute('normPar') or (not self.nO and p.getAttribute('free') == "1") else "0")
                        specOut.appendChild(Tools.parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                else:
                    if self.var and varIdx >= self.varValue:
                        Sources[sname]['free'] = True
                        #specOut.setAttribute('apply_edisp',ed)
                        for p in specPars:
                            freeFlag=("1" if p.getAttribute('name') == spec[0].getAttribute('normPar') else "0")
                            specOut.appendChild(Tools.parameter_element("%s"%freeFlag,"%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                    else:
                        Sources[sname]['free'] = False
                        #specOut.setAttribute('apply_edisp','false')
                        for p in specPars:
                            specOut.appendChild(Tools.parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                if Ext:
                    spatial = src.getElementsByTagName('spatialModel')
                    spatType = spatial[0].getAttribute('type')
                    spatPars = spatial[0].getElementsByTagName('parameter')
                    if str(spatType) == 'SpatialMap':
                        spatialOut.setAttribute('type','SpatialMap')
                        spatialOut.setAttribute('map_based_integral','true')
                        efile = os.path.join(self.extD, spatial[0].getAttribute('file'))
                        spatialOut.setAttribute('file',efile)
                        print 'Extended source {0} in ROI, make sure {1} is the correct path to the extended template.'.format(sname, efile)
                    else:#have to do above to get correct extended source template file localtion
                        spatialOut.setAttribute('type', str(spatType))
                        for p in spatPars:#for radial disks and gaussians, can just do the following
                            spatialOut.appendChild(Tools.parameter_element("0","%s"%str(p.getAttribute('name')),"%s"%str(p.getAttribute('max')),"%s"%str(p.getAttribute('min')),"%s"%str(p.getAttribute('scale')),"%s"%str(p.getAttribute('value'))))
                            print 'Extended source {0} in ROI, with {1} spatial model.'.format(sname, str(spatType))
                    
                    srcOut.setAttribute('type', 'DiffuseSource')
                    extSrcNum += 1
                    #print 'Extended source %s in ROI, make sure %s is the correct path to the extended template.'%(sname,efile)
                else:
                    spatialOut.setAttribute('type', 'SkyDirFunction')
                    spatialOut.appendChild(Tools.parameter_element("0", "RA", "360.0", "-360.0", "1.0", "%.4f"%srcRA))
                    spatialOut.appendChild(Tools.parameter_element("0", "DEC", "360.0", "-360.0", "1.0", "%.4f"%srcDEC))
                    srcOut.setAttribute('type','PointSource')
                    ptSrcNum += 1
                
                srcOut.appendChild(specOut)
                srcOut.appendChild(spatialOut)
                outputXml.documentElement.appendChild(srcOut)
        
        if self.GD is not None:
            gal = outputXml.createElement('source')
            gal.setAttribute('name', self.GDn)
            gal.setAttribute('type', 'DiffuseSource')
            galspec = outputXml.createElement('spectrum')
            galspec.setAttribute('type', 'PowerLaw')
            #galspec.setAttribute('apply_edisp','false')
            galspec.appendChild(Tools.parameter_element("1", "Prefactor", "10", "0", "1", "1"))
            if self.GIF:
                galspec.appendChild(Tools.parameter_element("1", "Index", "1", "-1", "1", "0"))
            else:
                galspec.appendChild(Tools.parameter_element("0", "Index", "1", "-1", "1", "0"))
            galspec.appendChild(Tools.parameter_element("0", "Scale", "1e6", "2e1", "1", "100"))
            galspatial = outputXml.createElement('spatialModel')
            galspatial.setAttribute('type', 'MapCubeFunction')
            galspatial.setAttribute('file', self.GD)
            galspatial.appendChild(Tools.parameter_element("0", "Normalization", "1e3", "1e-3", "1", "1"))
            gal.appendChild(galspec)
            gal.appendChild(galspatial)
            outputXml.documentElement.appendChild(gal)
    
        if self.ISO is not None:
            iso = outputXml.createElement('source')
            iso.setAttribute('name', self.ISOn)
            iso.setAttribute('type', 'DiffuseSource')
            isospec = outputXml.createElement('spectrum')
            isospec.setAttribute('type', 'FileFunction')
            isospec.setAttribute('file', self.ISO)
            isospec.setAttribute('apply_edisp','false')
            isospec.appendChild(Tools.parameter_element("1", "Normalization", "10", "0.01", "1", "1"))
            isospatial = outputXml.createElement('spatialModel')
            isospatial.setAttribute('type', 'ConstantValue')
            isospatial.appendChild(Tools.parameter_element("0", "Value", "10", "0", "1", "1"))
            iso.appendChild(isospec)
            iso.appendChild(isospatial)
            outputXml.documentElement.appendChild(iso)
        
        xmlStr = outputXml.toprettyxml(' ').splitlines(True)
        outStr = filter(lambda xmlStr: len(xmlStr) and not xmlStr.isspace(), xmlStr)
        outfile=open(self.out,'w')
        outfile.write(''.join(outStr))
        outfile.close()
        if not self.psF:
            print 'Added {0} point sources and {1} extended sources'.format(ptSrcNum, extSrcNum)
            if extSrcNum > 0:
                print 'If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument psForce=True'
        else:
            print 'Added {0} point sources, note that any extended sources in ROI were modeled as point sources becaue psForce option was set to True'.format(ptSrcNum)
        
        if self.reg:
            print "Building DS9 region file..."
            try:
                regFile = BuildRegion.BuildRegion(self.regFile, Sources, model_type='likelihood', frame=self.frame)
                print "Region built. File is located at {0}".format(regFile)
            except BuildRegionError as e:
                print e
        
        return outputXml, Sources


    def fits(self):
        model = minidom.getDOMImplementation().createDocument(None, 'source_library', None)
        model.documentElement.setAttribute('title', 'source_library')

        fgl = pyfits.open(self.srcs) #open source list file and access necessary fields, requires LAT source catalog definitions and names
        data = fgl['LAT_Point_Source_Catalog'].data
        extendedinfo = fgl['ExtendedSources'].data
        fgl.close()

        extName = extendedinfo['Source_Name']   
        extFile = extendedinfo['Spatial_Filename']
        extFunc = extendedinfo['Spatial_Function']
        extSize = extendedinfo['Model_SemiMajor']
        extRa = extendedinfo['RAJ2000']
        extDec = extendedinfo['DEJ2000']
        name = data['Source_Name']
        Sigvals = data['Signif_Avg']
        #VarIdx = data['Variability_Index']
        EName = data['Extended_Source_Name']
        ra = data['RAJ2000']
        dec = data['DEJ2000']
        glon = data['GLON']
        glat = data['GLAT']
        plflux = data['PL_Flux_Density']
        lpflux = data['LP_Flux_Density']
        coflux = data['PLEC_Flux_Density']
        pivot = data['Pivot_Energy']
        plIndex = data['PL_Index']
        lpIndex = data['LP_Index']
        lpbeta = data['LP_beta']
        plecIndex = data['PLEC_Index']
        plecexpFact = data['PLEC_Expfactor']
        plecexpIndex = data['PLEC_Exp_Index']
        spectype = data['SpectrumType']

        comment = model.createComment("Catalog Sources")
        model.documentElement.appendChild(comment)
        
        rad_bins = 5.
        step = (self.roi[2] + self.ER)/rad_bins
        radii = np.linspace(step, self.roi[2] + self.ER, num=rad_bins)
        ptSrcNum = 0
        extSrcNum = 0
        Sources={}#dictionary for sources, useful for creating region file later.

        # Calculate the distances first so we only loop over sources in each radii bin
        distances = Tools.angsep(self.roi[0], self.roi[1], ra, dec)

        for x in radii:
            if x == self.roi[2] + self.ER:
                comment = model.createComment('Sources between [{0},{1}] degrees of ROI center'.format(x - step, x))
                idx = np.logical_and(distances <= x, distances >= (x - step))
            else:
                comment = model.createComment('Sources between [{0},{1}) degrees of ROI center'.format(x-step , x))
                idx = np.logical_and(distances < x, distances >= (x - step))

            model.documentElement.appendChild(comment)

            for n, plf, lpf, cof, r, d, gl, gb, p, pli, lpi, lpb, pleci, plecef, plecei, t, TS, En, dist in zip(name[idx], plflux[idx], lpflux[idx], coflux[idx], ra[idx], dec[idx], glon[idx], glat[idx], pivot[idx], plIndex[idx], lpIndex[idx], lpbeta[idx], plecIndex[idx], plecexpFact[idx], plecexpIndex[idx], spectype[idx], Sigvals[idx], EName[idx], distances[idx]):
                source = model.createElement('source')
                vi = 0#remove later if we ever get a similar variability index thing going on
                E = (True if n[-1] == 'e' else False)
                if E and not self.psF:
                    Sources[En] = {'ra':r, 'dec':d, 'glon':gl, 'glat':gb, 'stype':t, 'E':E}
                    extSrcNum += 1
                    source.setAttribute('ROI_Center_Distance', "{0:.3f}".format(dist))
                    source.setAttribute('name', En)
                    source.setAttribute('type', "DiffuseSource")
                else:
                    if E and not self.E2C:#even if forcing all to point sources, use extended name except if E2CAT flag is set
                        Sources[En] = {'ra':r, 'dec':d, 'glon':gl, 'glat':gb, 'stype':t, 'E':E}
                        source.setAttribute('ROI_Center_Distance', "{0:.3f}".format(dist))
                        source.setAttribute('name', En)
                        source.setAttribute('type', "PointSource")
                    else:
                        Sources[n] = {'ra':r, 'dec':d, 'glon':gl, 'glat':gb, 'stype':t, 'E':E}
                        if oldNames:
                            srcname = '_' + srcname.replace(' ', '')
                            source.setAttribute('ROI_Center_Distance', "{0:.3f}".format(dist))
                            source.setAttribute('name', srcname)
                            source.setAttribute('type', "PointSource")
                        else:
                            source.setAttribute('ROI_Center_Distance', "{0:.3f}".format(dist))
                            source.setAttribute('name', n)
                            source.setAttribute('type', "PointSource")
                    ptSrcNum += 1
                
                if t == 'PowerLaw':
                    spec, free, comments = LikelihoodSpectra.PLspec(self.roi, self.radLim, self.maxRad, self.varValue, self.var, self.sig, self.nO, plf, pli, p, dist, TS, vi, False)
                elif t == 'PowerLaw2':#no value for flux from 100 MeV to 100 GeV in fits file
                    if pli != 1.:#so calculate it by integrating PowerLaw spectral model
                        F = plf*p**pli/(-pli+1.)*(1.e5**(-pli+1.)-1.e2**(-pli+1.))
                    else:
                        F=plf*p*np.log(1.e3)
                    spec, free, comments = LikelihoodSpectra.PL2spec(self.roi, self.radLim, self.maxRad, self.varValue, self.var, self.sig, self.nO, F, pli, dist, TS, vi)
                elif t == 'LogParabola':
                    spec, free, comments = LikelihoodSpectra.LPspec(self.roi, self.radLim, self.maxRad, self.varValue, self.var, self.sig, self.nO, lpf, lpi, p, lpb, dist, TS, vi)
                elif (t == 'PLSuperExpCutoff') or (t == 'PLSuperExpCutoff2'):
                    spec, free, comments = LikelihoodSpectra.CO2spec(self.roi, self.radLim, self.maxRad, self.varValue, self.var, self.sig, self.nO, cof, pleci, p, plecef, plecei, dist, TS, vi)
                else:
                    print "{0} has spectrum {1} which currently can't be modeled.".format(srcname, t)
                    spec = None
                    free = 0
                    comments = []
                    continue

                if E and not self.E2C:
                    Sources[En]['free'] = free
                else:
                    Sources[n]['free'] = free
                
                spatialModel = model.createElement('spatialModel')
                if E and not self.psF:
                    efile = None
                    efunc = None
                    eSize = None
                    eR = None
                    eD = None
                    for EXTNAME,EXTFILE,EXTFUNC,EXTSIZE,EXTRA,EXTDEC in zip(extName,extFile,extFunc,extSize,extRa,extDec):
                        if En == EXTNAME:
                            efunc = EXTFUNC
                            efunc = ('RadialGaussian' if efunc == 'RadialGauss' else efunc)
                            if efunc == 'SpatialMap':
                              efile = os.path.join(self.extD, EXTFILE)
                            else:
                              eSize = EXTSIZE
                              eR = EXTRA
                              eD = EXTDEC
                    if efunc == 'SpatialMap':
                        if efile == None:
                            print 'could not find a match for {0} in the list:'.format(En)
                            print extName
                            efile = ''
                            spatialModel.setAttribute('file', efile)
                            spatialModel.setAttribute('map_based_integral', "true")
                            spatialModel.setAttribute('type', "SpatialMap")
                            print 'Extended source {0} in ROI, make sure {1} is the correct path to the extended template.'.format(En, efile)
                        spatialModel.appendChild(Tools.parameter_element("0", "Prefactor", "1000.0", "0.001", "1.0", "1.0"))
                    else:
                        spatialModel.setAttribute('type', efunc)
                        spatialModel.appendChild(Tools.parameter_element("0", "RA", "360.0", "-360.0", "1.0", "{0}".format(eR)))
                        spatialMap.appendChild(Tools.parameter_element("0", "DEC", "90.0", "-90.0", "1.0", "{0}".format(eD)))
                        if efunc == 'RadialDisk':
                            spatialMap.appendChild(Tools.parameter_element("0", "Radius", "10.0", "0.0", "1.0", "{0}".format(eSize)))
                        else:
                            spatialMap.appendChild(Tools.parameter_element("0", "Sigma", "10.0", "0.0", "1.0", "{0}".format(eSize)))
                        
                        print 'Extended source {0} in ROI with {1} spatial model.'.format(En,efunc)
                else:
                    spatialModel.setAttribute('type', "SkyDirFunction")
                    spatialModel.appendChild(Tools.parameter_element("0", "RA", "360.0", "-360.0", "1.0", "{0}".format(r)))
                    spatialModel.appendChild(Tools.parameter_element("0", "DEC", "90.0", "-90.0", "1.0", "{0}".format(d)))

                for comment in comments:
                    source.appendChild(comment)
                
                source.appendChild(spec)
                source.appendChild(spatialModel)
                model.documentElement.appendChild(source)
        
        if not self.psF:
            print 'Added {0} point sources and {1} extended sources.'.format(ptSrcNum, extSrcNum)
            if extSrcNum > 0:
                print 'If using unbinned likelihood you will need to rerun gtdiffrsp for the extended sources or rerun the makeModel function with optional argument psForce=True.'
        else:
            print 'Added {0} point sources, note that any extended sources in ROI were modeled as point sources becaue psForce option was set to True.'.format(ptSrcNum)
        
        #add galactic diffuse with PL spectrum, fix index to zero for general use, those who want it to be free can unfreeze parameter manually
        comment = model.createComment("Diffuse Sources")
        model.documentElement.appendChild(comment)
        if self.GD is not None:
            source = model.createElement('source')
            source.setAttribute('name', self.GDn)
            source.setAttribute('type', 'DiffuseSource')
            spec = model.createElement('spectrum')
            spec.setAttribute('type', 'PowerLaw')
            spec.appendChild(Tools.parameter_element("1", "Prefactor", "10.0", "0.0", "1.0", "1.0"))
            if self.GIF:
                spec.appendChild(Tools.parameter_element("1", "Index", "1.0", "-1.0", "1.0", "0.0"))
            else:
                spec.appendChild(Tools.parameter_element("0", "Index", "1.0", "-1.0", "1.0", "0.0"))
        
            spec.appendChild(Tools.parameter_element("0", "Scale", "2e2", "5e1", "1.0", "1e2"))
            spatialModel = model.createElement('spatialModel')
            spatialModel.setAttribute('file', self.GD)
            spatialModel.setAttribute('type', "MapCubeFunction")
            spatialModel.appendChild(Tools.parameter_element("0", "Normalization", "1e3", "1e-3", "1.0", "1.0"))
            (src,) = (Name + spec + skydir,)
            galdiff = pS(src).getElementsByTagName('source')[0]
            galdiff.writexml(model)
            model.write('\n')
            source.appendChild(spec)
            source.appendChild(spatialModel)
            model.documentElement.appendChild(source)

        # Add isotropic diffuse
        if self.ISO is not None:
            source = model.createElement('source')
            source.setAttribute('name', self.ISOn)
            source.setAttribute('type', 'DiffuseSource')
            spec = model.createElement('spectrum')
            spec.setAttribute('type', 'FileFunction')
            spec.setAttribute('file', self.ISO)
            spec.setAttribute('apply_edisp', 'false')
            spec.appendChild(Tools.parameter_element("1", "Normalization", "10.0", "0.01", "1", "1"))
            spatialModel = model.createElement('spatialModel')
            spatialModel.setAttribute('type', 'ConstantValue')
            spatialModel.appendChild(Tools.parameter_element("0", "Value", "10.0", "0.0", "1.0", "1.0"))
            source.appendChild(spec)
            source.appendChild(spatialModel)
            model.documentElement.appendChild(source)
        
        if self.reg:
            print "Building DS9 region file..."
            try:
                regFile = BuildRegion.BuildRegion(self.regFile, Sources, model_type='likelihood', frame=self.frame)
                print "Region built. File is located at {0}".format(regFile)
            except BuildRegionError as e:
                print e
        
        return model, Sources

class simulation:
    def __init__(self, **params):
        for key, val in params.items():
            setattr(self, key, val)

        self.varValue = 72.44

        extension = os.path.splitext(self.srcs)[-1]
        if extension == '.xml':
            self.model, self.Sources = self.xml()
        elif extension in ['.fits', '.fit']:
            self.model, self.Sources = self.fits()
        else:
            raise IOError("{0} is not a compatible catalog.".format(self.srcs))

    def fits(self):
        """Generate the model from .fits catalog for use with gtobssim.
    
        Parameters
        ----------
        GD : str
            Path to galactic diffuse file
        GDn : str
            Galactic Diffuse name in model
        ISO : str
            Path to isotropic diffuse file
        ISOn : str
            Isotropic diffuse name in model
        oldNames : bool
            Name sources using old naming convention
        emin : float
            Minimum energy for integrated flux
        emax : float
            Maximum energy for integrated flux
        frame : str
            Coordinate frame for source directions
        apply_mask : bool
            Apply region of interest mask to diffuse emission models

        Returns
        -------
        model : file
            this model
        Sources : dict
            List of sources included in model

        """
        model = minidom.getDOMImplementation().createDocument(None,'source_library',None)
        model.documentElement.setAttribute('title', 'source_library')

        fgl = pyfits.open(self.srcs)
        data = fgl['LAT_Point_Source_Catalog'].data
        extendedinfo = fgl['ExtendedSources'].data
        fgl.close()

        extName = extendedinfo['Source_Name']
        extFile = extendedinfo['Spatial_Filename']
        extFunc = extendedinfo['Spatial_Function']
        extSemiMajor = extendedinfo['Model_SemiMajor']
        extSemiMinor = extendedinfo['Model_SemiMinor']
        extPosAng = extendedinfo['Model_PosAng']
        extRa = extendedinfo['RAJ2000']
        extDec = extendedinfo['DEJ2000']
        name = data['Source_Name']
        Sigvals = data['Signif_Avg']
        EName = data['Extended_Source_Name']
        ra = data['RAJ2000']
        dec = data['DEJ2000']
        glat = data['GLAT']
        glon = data['GLON']
        plflux = data['PL_Flux_Density']
        lpflux = data['LP_Flux_Density']
        coflux = data['PLEC_Flux_Density']
        pivot = data['Pivot_Energy']
        plIndex = data['PL_Index']
        lpIndex = data['LP_Index']
        lpbeta = data['LP_beta']
        plecIndex = data['PLEC_Index']
        plecexpFact = data['PLEC_Expfactor']
        plecexpIndex = data['PLEC_Exp_Index']
        spectype = data['SpectrumType']

        Sources = {} # Dictionary for sources
        ptSrcNum = 0 # Counter for point sources
        extSrcNum = 0 # Counter for extended sources

        rad_bins = 5    
        step = (self.roi[2] + self.ER)/float(rad_bins) # Divide ROI radius plus ExtraRadius degrees into 5 steps for ordering of sources
        radii = np.linspace(step, self.roi[2] + self.ER, num=rad_bins)

        # Calculate the distances first so we only loop over sources in each radii bin
        distances = Tools.angsep(self.roi[0], self.roi[1], ra, dec)
        for x in radii:
            if x == self.roi[2] + self.ER:
                comment = model.createComment("Sources between [{0},{1}] degrees of ROI center".format(x - step, x))
                idx = np.logical_and(distances <= x, distances >= x - step)
            else:
                comment = model.createComment("Sources between [{0},{1})".format(x - step, x))
                idx = np.logical_and(distances < x, distances >= x - step)
            model.documentElement.appendChild(comment)
            for n, plf, lpf, cof, r, d, gb, gl, p, pli, lpi, lpb, pleci, plecef, plecei, t, En, dist in zip(name[idx], plflux[idx], lpflux[idx], coflux[idx], ra[idx], dec[idx], glat[idx], glon[idx], pivot[idx], plIndex[idx], lpIndex[idx], lpbeta[idx], plecIndex[idx], plecexpFact[idx], plecexpIndex[idx], spectype[idx], EName[idx], distances[idx]):
                vi = 0 #remove later if we ever get a similar variability index thing going on
                E = (True if n[-1] == 'e' else False)

                # Determine the name of the source given the input flags and get the efile if it exists
                srcname = n
                efile = ''
                semimajor = None
                semiminor = None
                posang = None
                spatialfunc = ''
                if E and not sL.psF:
                    name_idx = np.where(extName == En)[0]
                    if len(name_idx) == 0:
                        print 'coult not find a match for {0} in the list:'.format(En)
                        print extName
                        print 'Skipping source...'
                        continue
                    elif len(name_idx) > 1:
                        print 'Found multiple extended sources with the name {0} in list:'.format(En)
                        print extName
                        print 'Skipping source...'
                        continue
                    else:
                    # efile = sL.extD + extFile[idx[0]]
                        efile = extFile[name_idx[0]]
                        semimajor = extSemiMajor[name_idx[0]]
                        semiminor = extSemiMinor[name_idx[0]]
                        posang = extPosAng[name_idx[0]]
                        spatialfunc = extFunc[name_idx[0]]

                    if not self.E2C: #even if forcing all to point sources, use extended name except if E2CAT flag is set
                        srcname = En.replace(" ", "_").replace(".", "d").replace("+", "p").replace("-", "m")
                    else:
                        if self.oldNames:
                            srcname = '_' + n.replace(" ","_").replace(".", "d").replace("+", "p").replace("-", "m")
                        else:
                            srcname = En.replace(" ", "_").replace(".", "d").replace("+", "p").replace("-", "m")
                else:
                    if self.oldNames:
                        srcname = '_' + srcname.replace(" ","_").replace(".","d").replace("+","p").replace("-","m")

                if E and not self.psF:
                    try:
                        source, modeled_extended = SimulationSources.AddExtendedSource(srcname, t, spatialfunc, directory=self.wd, extDir=self.extD, ra=r, dec=d, glon=gl, glat=gb, major_axis=semimajor, minor_axis=semiminor, position_angle=posang, efile=efile, emin=self.emin, emax=self.emax, frame=self.frame, resolution=self.extSrcRes, pivot_energy=p, pl_flux_density=plf, lp_flux_density=lpf, plec_flux_density=cof, pl_index=pli, lp_index=lpi, lp_beta=lpb, plec_index=pleci, plec_expfactor=plecef, plec_exp_index=plecei)
                        model.documentElement.appendChild(source)
                        if modeled_extended:
                            Sources[srcname] = {'ra': r, 'dec':d, 'glon':gl, 'glat':gb, 'stype':t, 'Spatial_Function':spatialfunc, 'extFile':efile, 'E':True}
                            extSrcNum += 1
                        else:
                            Sources[srcname] = {'ra': r, 'dec':d, 'glon':gl, 'glat':gb, 'stype':t, 'Spatial_Function':spatialfunc, 'extFile':efile, 'E':False}
                            ptSrcNum += 1
                    except AddSourceError as e:
                        print "Error encountered when adding extended source {0}.".format(srcname)
                        print e
                        print "Skipping source..."
                        continue
                else:
                    try:
                        source = SimulationSources.AddPointSource(srcname, t, self.emin, self.emax, self.wd, ra=r, dec=d, glon=gl, glat=gb, frame=self.frame, pivot_energy=p, pl_flux_density=plf, lp_flux_density=lpf, plec_flux_density=cof, pl_index=pli, lp_index=lpi, lp_beta=lpb, plec_index=pleci, plec_expfactor=plecef, plec_exp_index=plecei)
                        model.documentElement.appendChild(source)
                        Sources[srcname] = {'ra':r, 'dec':d, 'glon':gl, 'glat':gb, 'stype':t, 'E':E}
                        ptSrcNum += 1
                    except AddSourceError as e:
                        print "Error encountered when adding point source {0}.".format(srcname)
                        print e
                        print "Skipping source..."
                        continue

        # Add galactic diffuse emission
        if self.GDn is not None:
            headers_present = False
            try:
                headers_present = Tools.header_check(self.GD)
            except DiffuseHeaderError as e:
                print e

            if headers_present:
                if self.GDflux is None: 
                    if os.path.basename(GD) == 'gll_iem_v06.fits':
                        self.GDflux = 0.0006728539887251224
                    elif (os.path.basename(GD) == 'gll_iem_v07.fits') or (os.path.basename(GD) == 'gll_iem_v07_revised.fits'):
                        self.GDflux = 0.0008463167432920544
                    else:
                        self.GDflux = 0.

                if apply_mask:
                    gd, ext = os.path.splitext(os.path.basename(self.GD))
                    self.GD, self.GDflux = maskFits.MaskFits(self.GD, out=os.path.join(self.wd, gd + "_masked" + ext), mask_type='radial', radius=self.roi[2] + self.ER, radius2=None, angle=0., center=(self.roi[0], self.roi[1]), clobber=True)
                    print "Masking Galactic Diffuse model..."
                    print "Applying mask of {0} degrees around ({1}, {2})".format(self.roi[2], self.roi[0], self.roi[1])
                
                source = model.createElement('source')
                source.setAttribute('name', self.GDn)
                spec = model.createElement('spectrum')
                spec.setAttribute('escale', "MeV")

                spectrumClass = model.createElement('SpectrumClass')
                spectrumClass.setAttribute('name', "MapCube")
                spectrumClass.setAttribute('params', "flux={0},fitsFile={1}".format(self.GDflux*1.e4, self.GD))

                use_spectrum = model.createElement('use_spectrum')
                use_spectrum.setAttribute('frame', "galaxy")

                spec.appendChild(spectrumClass)
                spec.appendChild(use_spectrum)

                comment = model.createComment("This is v07 of the diffuse emission model. Integrated flux from mapcube is {0} (#/m^2/s)".format(GDflux*1.e4))
                source.appendChild(comment)

                source.appendChild(spec)
                model.documentElement.appendChild(source)
                Sources[self.GDn] = {'flux':self.GDflux, 'SpatialFunction':'MapCube'}
                extSrcNum += 1

        # Add isotropic diffuse model
        if self.ISOn is not None:
            source = model.createElement('source')
            source.setAttribute('name', self.ISOn)
            spec = model.createElement('spectrum')
            spec.setAttribute('escale', "MeV")
            spectrumClass = model.createElement('SpectrumClass')
            spectrumClass.setAttribute('name', "FileSpectrumMap")

            ISOpath = "$(FERMI_DIR)/refdata/fermi/galdiffuse/isotropic_allsky.fits"
            if self.ISOflux is None:
                self.ISOflux = 0.0

            spectrumClass.setAttribute('params', "flux={0},fitsFile={1},specFile={2}".format(self.ISOflux*1.e4, ISOpath, self.ISO))

            use_spectrum = model.createElement('use_spectrum')
            use_spectrum.setAttribute('frame', "galaxy")

            comment = model.createComment("This is the isotropic diffuse spectrum. Integrated flux is {0} (#/m^2/s) [Note units]. When a flux of 0 is given the integral is calculated automatically and used.".format(self.ISOflux*1.e4))

            spec.appendChild(spectrumClass)
            spec.appendChild(use_spectrum)

            source.appendChild(comment)
            source.appendChild(spec)
            model.documentElement.appendChild(source)
            Sources[self.ISOn] = {'flux':self.ISOflux, 'SpatialFunction':'FileSpectrumMap'}
            extSrcNum += 1

        if not sL.psF:
            print 'Added {0} point sources and {1} extended sources.'.format(ptSrcNum, extSrcNum)
        else:
            print 'Added {0} point sources, note that any extended sources in ROI were modeled as point sources because psForce option was set to True.'.format(ptSrcNum)

        if self.reg:
            print "Building DS9 region file..."
            try:
                regFile = BuildRegion.BuildRegion(self.regFile, Sources, model_type='simulation', frame=self.frame)
                print "Region built. File is located at {0}".format(regFile)
            except BuildRegionError as e:
                print e

        return model, Sources

    def xml(self):
        print "I can't do this yet"