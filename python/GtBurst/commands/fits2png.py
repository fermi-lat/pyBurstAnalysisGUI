#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
from GtBurst import aplpy
import matplotlib.pyplot as plt
import sys
    
def fitsToPNG(fitsfile,pngfile,vmin=None,vmax=None):
    #Display the TS map    
    #figure                      = plt.figure()
    tsfig                       = aplpy.FITSFigure(fitsfile,convention='calabretta')
    tsfig.set_tick_labels_font(size='small')
    tsfig.set_axis_labels_font(size='small')
    if(vmin!=None and vmax!=None):
      tsfig.show_colorscale(cmap='gist_heat',aspect='auto',vmin=float(vmin),vmax=float(vmax))
    else:
      tsfig.show_colorscale(cmap='gist_heat',aspect='auto')
    # Modify the tick labels for precision and format
    tsfig.tick_labels.set_xformat('ddd.dd')
    tsfig.tick_labels.set_yformat('ddd.dd')
    
    # Display a grid and tweak the properties
    tsfig.show_grid()
    tsfig.add_colorbar()
    
    #figure.canvas.draw()
    tsfig.save(pngfile)
pass

def fitsToPNGembedded(fitsfile,pngfile='__png',vmin=None,vmax=None):
    fitsToPNG(fitsfile,pngfile,vmin,vmax)
    data_uri                  = open(pngfile, 'rb').read().encode('base64').replace('\n', '')
    #img_tag                   = '<img src="data:image/png;base64,{0}">'.format(data_uri)
    return data_uri
    
pass

if __name__=='__main__':
  #Get all key=value pairs as a dictionary
  fitsfile                    = sys.argv[1]
  pngfile                     = sys.argv[2]
  fitsToPNG(fitsfile,pngfile)
