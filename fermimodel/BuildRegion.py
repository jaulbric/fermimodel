from astropy.coordinates import SkyCoord

from Exceptions import BuildRegionError
# class BuildRegionError(Exception):
#     """Raise when the region cannot be built"""
#     pass

class Region:
    def __init__(self, regFile, Sources, frame='fk5', model_type='likelihood'):
        self.regFile = regFile
        self.Sources = Sources

    def build(self, frame, model_type):
        """Build a DS9 region file

        Parameters
        ----------
        model_type : str
            Model type that was built. Must be 'likelihood' or 'simulation'

        Returns
        -------
        regFile : str
            Path to region file

        """ 
        myreg = open(self.regFile,'w')#note that this will overwrite previous region files of the same name
        myreg.write('# Region File format: DS9 version 8.0.1')#I don't actually know, but I think it's one of the later ones, need to verify
        myreg.write('\n# Created by fermimodel')
        myreg.write('\nglobal font="roman 10 normal" move =0')
        for k, src in self.Sources.items():
            if src['diffuse']:
                continue

            if frame in ['galactic', 'GALACTIC']:
                xcoord = src['glon']
                ycoord = src['glat']
                use_frame = 'GALACTIC'
            elif frame in ['fk5', 'FK5', 'J2000']:
                xcoord = src['ra']
                ycoord = src['dec']
                use_frame = 'FK5'
            elif frame == 'icrs':
                c = SkyCoord(ra=src['ra'], dec=src['dec'], frame='fk5', unit='degree')
                xcoord = c.icrs.ra.degree
                ycoord = c.icrs.dec.degree
                use_frame = 'ICRS'
            else:
                c = SkyCoord(ra=src['ra'], dec=src['dec'], frame='fk5', unit='degree')
                xcoord = float(c.transform_to(frame).to_string().split(' ')[0])
                ycoord = float(c.transfrom_to(frame).to_string().split(' ')[1])
                use_frame = frame

            if model_type == 'likelihood':
                #get color based on if the source is free or not
                color = ('green' if src['free'] else 'magenta')
                if src['E']:#if the source is extended, always have the point be a "big" box
                    myreg.write('\n{0};point({1:.3f},{2:.3f}) # point = box 18 color = {3} text={{4}}'.format(use_frame, xcoord, ycoord, color, k))
                else:#if the source is a point source, choose the point type based on spectral model
                    ptype = ('cross' if (src['stype'] == 'PLSuperExpCutoff' or src['stype'] == 'PLSuperExpCutoff2') else 'diamond' if src['stype'] == 'LogParabola' else 'circle')
                    myreg.write('\n{0};point({1:.3f},{2:.3f}) # point = {3} 15 color = {4} text={{5}}'.format(use_frame, xcoord, ycoord, ptype, color, k))
            elif model_type == 'simulation':
                color = 'green'
                if src['E']:
                    myreg.write('\n{0};point({1:.3f},{2:.3f}) # point = box 18 color = {3} text={{4}}'.format(use_frame, xcoord, ycoord, color, k))
                else:
                    ptype = ('cross' if (src['stype'] == 'PLSuperExpCutoff' or src['stype'] == 'PLSuperExpCutoff2') else 'diamond' if src['stype'] == 'LogParabola' else 'circle')
                    myreg.write('\n{0};point({1:.3f},{2:.3f}) # point = {3} 15 color = {4} text={{5}}'.format(use_frame, xcoord, ycoord, ptype, color, k))
            else:
                raise BuildRegionError("Cannot build region for model type {0}".format(model_type))
        
        myreg.close()
        return self.regFile

def buildRegion(regFile, Sources, frame='fk5', model_type='likelihood'):
    region = Region(regFile, Sources, frame='fk5', model_type='likelihood')
    return region.build(frame, model_type)