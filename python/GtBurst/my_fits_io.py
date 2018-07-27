# If pywcs is installed, use pyfits, otherwise use
# astropy.io.fits

try:

    import astropy.io.fits as pyfits

except:
    
    import pywcs
    import pyfits
