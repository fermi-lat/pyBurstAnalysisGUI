try:
    # Try and load the aplpy from the environment
    from aplpy import *

except ImportError:
    
    # Use internal aplpy
    
    from my_aplpy import *
