# import Model
# import AddCatalogSources
# import BuildRegion
# import LikelihoodSpectra
# import SimulationSources
# import Tools

class GetFluxError(Exception):
    """Raise when we cannot calculate the flux of the source from parameters"""
    pass

class WriteSpectrumError(Exception):
    """Raised when spectrum type cannot be written to file."""
    pass

class HeaderCheckError(Exception):
    """Raised when diffuse emission is missing essential headers."""
    pass

class ExtendedTemplateError(Exception):
    """Raise when an extended template cannot be found."""
    pass

class BuildRegionError(Exception):
    """Raise when the region cannot be built"""
    pass

class AddSourceError(Exception):
    """Raised when model cannot add a source."""
    pass

class SpectrumError(Exception):
    """Raised when input spectrum type does not match allowed spectrum types"""
    pass