#! /usr/bin/env python

import sys
import os
from GtBurst import commandDefiner
import pyfits, numpy,math
import re

################ Command definition #############################
executableName                = "gtdotsmap"
version                       = "1.0.0"
shortDescription              = "Produce a TS map"
author                        = "G.Vianello, giacomov@slac.stanford.edu"
thisCommand                   = commandDefiner.Command(executableName,shortDescription,version,author)

#Define the command parameters
thisCommand.addParameter("filteredeventfile","Input event list (FT1 file)",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="fit")
thisCommand.addParameter("rspfile","LAT response (RSP file)",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="rsp")
thisCommand.addParameter("ft2file","Spacecraft file (FT2)",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("tsltcube","Pre-computed livetime cube",commandDefiner.OPTIONAL,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("tsexpomap","Pre-computed exposure map",commandDefiner.OPTIONAL,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("xmlmodel","XML model",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("tsmap","Name for the output file for the TS map",commandDefiner.MANDATORY,partype=commandDefiner.OUTPUTFILE,extension="fits")
#thisCommand.addParameter("skymap","Name for the for the sky map (needed only if you want to plot your results)",commandDefiner.OPTIONAL,partype=commandDefiner.INPUTFILE,extension="fit")
#thisCommand.addParameter("tsmin","Minimum TS to consider a source detected",commandDefiner.OPTIONAL,20)
#thisCommand.addParameter("optimizeposition","Optimize position?",commandDefiner.OPTIONAL,"yes",possiblevalues=['yes','no'])
#thisCommand.addParameter("showmodelimage","Show an image representing the best fit likelihood model?",commandDefiner.OPTIONAL,"yes",possiblevalues=['yes','no'])
thisCommand.addParameter("step","Size of the grid step (deg)",commandDefiner.OPTIONAL,0.8)
thisCommand.addParameter("side","Number of steps per side of the TS map",commandDefiner.OPTIONAL,'auto',partype=commandDefiner.HIDDEN)
thisCommand.addParameter("clobber","Overwrite output file? (possible values: 'yes' or 'no')",commandDefiner.OPTIONAL,"yes")
thisCommand.addParameter("verbose","Verbose output (possible values: 'yes' or 'no')",commandDefiner.OPTIONAL,"yes")
thisCommand.addParameter("figure","Matplotlib figure for the interactive mode",commandDefiner.OPTIONAL,None,partype=commandDefiner.PYTHONONLY)

GUIdescription                = "Here you will build a TS map of the ROI you defined in the first step,"
GUIdescription               += " using the model you selected in the 2nd step. A likelihood is performed in each"
GUIdescription               += " point of a grid of coordinates, then the data are tested for a source at that"
GUIdescription               += " coordinates. The TS is the difference in log Likelihood between the model without"
GUIdescription               += " the source and the model with the source."
GUIdescription               += "TIP The TS map should take between 5 and 10 minutes to complete."
thisCommand.setGUIdescription(GUIdescription)

##################################################################

def _yesOrNoToBool(value):      
  if(value.lower()=="yes"):
    return True
  elif(value.lower()=="no"):
    return False
  else:
    raise ValueError("Unrecognized clobber option. You can use 'yes' or 'no'")    
  pass
pass

class Message(object):
  def __init__(self,verbose):
    self.verbose              = bool(verbose)
  pass
  
  def __call__(self,string):
    if(self.verbose):
      print(string)
pass   

def gtdotsmap(**kwargs):
  run(**kwargs)
pass

def run(**kwargs):
  if(len(kwargs.keys())==0):
    #Nothing specified, the user needs just help!
    thisCommand.getHelp()
    return
  pass
  
  #Get parameters values
  thisCommand.setParValuesFromDictionary(kwargs)
  try:
    eventfile                   = thisCommand.getParValue('filteredeventfile')
    rspfile                     = thisCommand.getParValue('rspfile')
    ft2file                     = thisCommand.getParValue('ft2file')
    xmlmodel                    = thisCommand.getParValue('xmlmodel')
    tsltcube                    = thisCommand.getParValue('tsltcube')
    tsexpomap                   = thisCommand.getParValue('tsexpomap')
    tsmap                       = thisCommand.getParValue('tsmap')
    step                        = float(thisCommand.getParValue('step'))
    side                        = thisCommand.getParValue('side')
    if(side=='auto'):
      side                      = None
#    showmodelimage              = thisCommand.getParValue('showmodelimage')
#    optimize                    = thisCommand.getParValue('optimizeposition')
#    tsmin                       = float(thisCommand.getParValue('tsmin'))
#    skymap                      = thisCommand.getParValue('skymap')
    clobber                     = _yesOrNoToBool(thisCommand.getParValue('clobber'))
    verbose                     = _yesOrNoToBool(thisCommand.getParValue('verbose'))
    figure                      = thisCommand.getParValue('figure')
  except KeyError as err:
    print("\n\nERROR: Parameter %s not found or incorrect! \n\n" %(err.args[0]))
    
    #Print help
    print thisCommand.getHelp()
    return
  pass
  
  from GtBurst import dataHandling
  from GtBurst.angularDistance import getAngularDistance
  
  origra                      = float(dataHandling._getParamFromXML(xmlmodel,'RA'))
  origdec                     = float(dataHandling._getParamFromXML(xmlmodel,'DEC'))
  sourceName                  = dataHandling._getParamFromXML(xmlmodel,'OBJECT')
  
  LATdata                     = dataHandling.LATData(eventfile,rspfile,ft2file)
  tsmap                       = LATdata.makeTSmap(xmlmodel,sourceName,step,side,tsmap,tsltcube,tsexpomap)
  tsltcube                    = LATdata.livetimeCube
  tsexpomap                   = LATdata.exposureMap
  
  ra,dec,tsmax                = findMaximumTSmap(tsmap,tsexpomap)
  
  print("\nCoordinates of the maximum of the TS map in the allowed region (TS = %.1f):" %(tsmax))
  print("(R.A., Dec.)              = (%6.3f, %6.3f)\n" %(ra,dec))
  print("Distance from ROI center  = %6.3f\n\n" %(getAngularDistance(origra,origdec,ra,dec)))

  if(figure!=None):
    from GtBurst import aplpy   
    #Display the TS map    
    figure.clear()
    tsfig                       = aplpy.FITSFigure(tsmap,convention='calabretta',figure=figure)
    tsfig.set_tick_labels_font(size='small')
    tsfig.set_axis_labels_font(size='small')
    tsfig.show_colorscale(cmap='gist_heat',aspect='auto')
    tsfig.show_markers([ra], [dec], edgecolor='green', facecolor='none', marker='o', s=10, alpha=0.5)
    # Modify the tick labels for precision and format
    tsfig.tick_labels.set_xformat('ddd.dd')
    tsfig.tick_labels.set_yformat('ddd.dd')
    
    # Display a grid and tweak the properties
    tsfig.show_grid()
    
    figure.canvas.draw()
  pass
  
  return 'tsmap', tsmap, 'tsmap_ra', ra, 'tsmap_dec', dec, 'tsmap_maxTS', tsmax, 'tsltcube', tsltcube, 'tsexpomap', tsexpomap
pass

def findMaximumTSmap(tsmap,tsexpomap):
  #Find the maximum of the TS map
  import pywcs
  f                           = pyfits.open(tsmap)
  image                       = f[0].data
  wcs                         = pywcs.WCS(f[0].header)
  f.close()
    
  #Position of the maximum
  idxs                        = numpy.unravel_index(image.argmax(), image.shape)
  #R.A., Dec of the maximum (the +1 is due to the FORTRAN Vs C convention
  ra,dec                      = wcs.wcs_pix2sky(idxs[1]+1,idxs[0]+1,1)
  ra,dec                      = ra[0],dec[0]
  
  #Now check that the value in the exposure map for this ra,dec is not too small,
  #nor that this Ra,Dec is at the margin of an excluded zones
  #(this avoid triggering on the Earth limb when strategy=events)
  fexp                        = pyfits.open(tsexpomap)
  expmap                      = fexp[0].data[0]
  wcsexp                      = pywcs.WCS(fexp[0].header)
  fexp.close()
  
  while(1==1):
    #Note that the value of one pixel is valid from .5 to 1.5
    #Note also that the exposure map is larger than the TS map
    pixels                    = wcsexp.wcs_sky2pix([[ra,dec,1]],1)[0]
    exposureHere              = expmap[pixels[1]-0.5,pixels[0]-0.5]
    exposureUp                = expmap[pixels[1]-0.5-1,pixels[0]-0.5]
    exposureDown              = expmap[pixels[1]-0.5+1,pixels[0]-0.5]
    exposureRight             = expmap[pixels[1]-0.5,pixels[0]-0.5-1]
    exposureLeft              = expmap[pixels[1]-0.5,pixels[0]-0.5+1]
    #print("Exposure: %s" %(exposureHere))
    if(exposureHere > 0 and 
       exposureUp > 0 and 
       exposureDown > 0 and 
       exposureLeft > 0 and 
       exposureRight > 0):
      break
    else:
      #Mask out this value
      print("Neglecting maximum at %s,%s because of low exposure there..." %(ra,dec))
      image[idxs[0],idxs[1]]  = 0.0
      idxs                        = numpy.unravel_index(image.argmax(), image.shape)
      #R.A., Dec of the maximum (the +1 is due to the FORTRAN Vs C convention
      ra,dec                      = wcs.wcs_pix2sky(idxs[1]+1,idxs[0]+1,1)
      ra,dec                      = ra[0],dec[0]
      continue
    pass
  pass
  
  tsmax                       = image.max()
  return ra,dec, tsmax
pass

thisCommand.run = run

if __name__=='__main__':
  #Get all key=value pairs as a dictionary
  args                           = dict(arg.split('=') for arg in sys.argv[1:])
  gtdotsmap(**args)
pass
