class FermiModelError(Exception):
	"""General exception for the package"""
	pass

class GetFluxError(FermiModelError):
    """Raise when we cannot calculate the flux of the source from parameters"""
    pass

class WriteSpectrumError(FermiModelError):
    """Raised when spectrum type cannot be written to file."""
    pass

class HeaderCheckError(FermiModelError):
    """Raised when diffuse emission is missing essential headers."""
    pass

class ExtendedTemplateError(FermiModelError):
    """Raise when an extended template cannot be found."""
    pass

class BuildRegionError(FermiModelError):
    """Raise when the region cannot be built"""
    pass

class AddSourceError(FermiModelError):
    """Raised when model cannot add a source."""
    pass

class SpectrumError(FermiModelError):
    """Raised when input spectrum type does not match allowed spectrum types"""
    pass