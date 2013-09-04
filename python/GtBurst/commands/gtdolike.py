#! /usr/bin/env python

import sys
import os
from GtBurst import commandDefiner
from GtBurst.GtBurstException import GtBurstException

import pyfits, numpy,math
import re

################ Command definition #############################
executableName                = "gtdolike"
version                       = "1.0.0"
shortDescription              = "Perform an unbinned likelihood analysis"
author                        = "G.Vianello, giacomov@slac.stanford.edu"
thisCommand                   = commandDefiner.Command(executableName,shortDescription,version,author)

#Define the command parameters
thisCommand.addParameter("filteredeventfile","Input event list (FT1 file)",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="fit")
thisCommand.addParameter("rspfile","LAT response (RSP file)",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="rsp")
thisCommand.addParameter("ft2file","Spacecraft file (FT2)",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("expomap","pre-computed exposure map",commandDefiner.OPTIONAL,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("ltcube","pre-computed livetime cube",commandDefiner.OPTIONAL,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("xmlmodel","XML model",commandDefiner.MANDATORY,partype=commandDefiner.DATASETFILE,extension="fits")
thisCommand.addParameter("skymap","Name for the sky map (needed only if you want to plot your results)",commandDefiner.OPTIONAL,partype=commandDefiner.INPUTFILE,extension="fit")
thisCommand.addParameter("tsmin","Minimum TS to consider a source detected",commandDefiner.OPTIONAL,20)
thisCommand.addParameter("optimizeposition","Optimize position?",commandDefiner.OPTIONAL,"yes",possiblevalues=['yes','no'])
thisCommand.addParameter("showmodelimage","Show an image representing the best fit likelihood model?",commandDefiner.OPTIONAL,"yes",possiblevalues=['yes','no'])
thisCommand.addParameter("spectralfiles","Produce spectral files for XSPEC?",commandDefiner.OPTIONAL,"yes",possiblevalues=['yes','no'])

#thisCommand.addParameter("irf","Data class (TRANSIENT or SOURCE)",commandDefiner.MANDATORY,'TRANSIENT',possiblevalues=['TRANSIENT','SOURCE'])
thisCommand.addParameter("clobber","Overwrite output file? (possible values: 'yes' or 'no')",commandDefiner.OPTIONAL,"yes")
thisCommand.addParameter("verbose","Verbose output (possible values: 'yes' or 'no')",commandDefiner.OPTIONAL,"yes")
thisCommand.addParameter("figure","Matplotlib figure for the interactive mode",commandDefiner.OPTIONAL,None,partype=commandDefiner.PYTHONONLY)

GUIdescription                = "Here you will perform an unbinned likelihood analysis on the data you selected in the first step,"
GUIdescription               += " using the model you selected in the 2nd step."
GUIdescription               += "TIP The likelihood analysis should take between 5 and 10 minutes to complete."
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

def gtdolike(**kwargs):
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
    expomap                     = thisCommand.getParValue('expomap')
    ltcube                      = thisCommand.getParValue('ltcube')
    xmlmodel                    = thisCommand.getParValue('xmlmodel')
    showmodelimage              = thisCommand.getParValue('showmodelimage')
    optimize                    = thisCommand.getParValue('optimizeposition')
    spectralfiles               = thisCommand.getParValue('spectralfiles')
    tsmin                       = float(thisCommand.getParValue('tsmin'))
    skymap                      = thisCommand.getParValue('skymap')
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
  
  LATdata                     = dataHandling.LATData(eventfile,rspfile,ft2file)
  try:
    outfilelike, sources        = LATdata.doUnbinnedLikelihoodAnalysis(xmlmodel,tsmin,expomap=expomap,ltcube=ltcube)
  except GtBurstException as gt:
    raise gt
  except:
    raise
  
  #Transfer information on the source from the input to the output XML
  irf                         = dataHandling._getParamFromXML(xmlmodel,'IRF')
  ra                          = dataHandling._getParamFromXML(xmlmodel,'RA')
  dec                         = dataHandling._getParamFromXML(xmlmodel,'DEC')
  name                        = dataHandling._getParamFromXML(xmlmodel,'OBJECT')
  
  grb                         = filter(lambda x:x.name.lower().find(name.lower())>=0,sources)[0]
  grb_TS                      = grb.TS
  
  if(irf==None):
    print("\n\nWARNING: could not read IRF from XML file. Be sure you know what you are doing...")
  else:
    dataHandling._writeParamIntoXML(outfilelike,IRF=irf,OBJECT=name,RA=ra,DEC=dec)
  pass

  
  if(spectralfiles=='yes'):
    phafile,rspfile,bakfile   = LATdata.doSpectralFiles(outfilelike)
  pass
  
  localizationMessage         = ''
  bestra                      = ''
  bestdec                     = ''
  poserr                      = ''
  distance                    = ''
  if(optimize=='yes'):
    sourceName                = name
    
    #If the TS for the source is above tsmin, optimize its position
    grb                       = filter(lambda x:x.name.lower().find(sourceName.lower())>=0,sources)[0]
    if(math.ceil(grb.TS) >= tsmin):
      try:
        bestra,bestdec,poserr = LATdata.optimizeSourcePosition(outfilelike,sourceName)
      except:
        raise GtBurstException(207,"gtfindsrc execution failed. Were the source detected in the likelihood step?")
      else:
        localizationMessage += "\nNew localization from gtfindsrc:\n\n"
        localizationMessage += "(R.A., Dec)                     = (%6.3f, %6.3f)\n" %(bestra,bestdec)
        localizationMessage += "68 %s containment radius        = %6.3f\n" %('%',poserr)
        distance             = getAngularDistance(float(ra),float(dec),float(bestra),float(bestdec))
        localizationMessage += "Distance from initial position  = %6.3f\n\n" %(distance)
        localizationMessage += "NOTE: this new localization WILL NOT be used by default. If you judge"
        localizationMessage += " it is a better localization than the one you started with, update the"
        localizationMessage += " coordinates yourself and re-run the likelihood\n"
    pass
  pass
  
  if(figure!=None and skymap!=None and showmodelimage=='yes'):
    
    #Now produce the binned exposure map (needed in order to display the fitted model as an image)
    modelmapfile              = LATdata.makeModelSkyMap(outfilelike) 
    
    #Display point sources in the image, and report in the table
    #all 2FGL sources with TS > 9 + always the GRB, independently of its TS
    detectedSources           = []
    grbFlux                   = 1e-13
    for src in sources:
      if(src.TS > 20):
        weight                = 'bold'
      else:
        weight                = 'regular'
      pass
      
      if(src.type=='PointSource'):
        if(src.TS > 4):
          detectedSources.append(src)
          if(src.name.find('2FGL')<0):
            #GRB
            grbFlux           = src.flux
        pass
      pass
    pass
    
    #Display the counts map
    from GtBurst import aplpy
    import matplotlib.pyplot as plt
    
    figure.clear()
    orig                     = aplpy.FITSFigure(skymap,convention='calabretta',
                                                figure=figure,subplot=[0.1,0.1,0.45,0.7])
    vmax                     = max(pyfits.open(skymap)[0].data.flatten())
    nEvents                  = numpy.sum(pyfits.open(skymap)[0].data)
    telapsed                 = pyfits.open(eventfile)[0].header['TSTOP']-pyfits.open(eventfile)[0].header['TSTART']
    orig.set_tick_labels_font(size='small')
    orig.set_axis_labels_font(size='small')
    orig.show_colorscale(cmap='gist_heat',vmin=0.1,vmax=vmax,stretch='log',aspect='auto')
    # Modify the tick labels for precision and format
    orig.tick_labels.set_xformat('ddd.dd')
    orig.tick_labels.set_yformat('ddd.dd')
    
    # Display a grid and tweak the properties
    orig.show_grid()
    
    #Renormalize the modelmapfile to the flux of the grb
    f                              = pyfits.open(modelmapfile,'update')
    f[0].data                      = f[0].data/numpy.max(f[0].data)*nEvents/telapsed
    print(numpy.max(f[0].data))
    f.close()
    
    img                      = aplpy.FITSFigure(modelmapfile,convention='calabretta',
                                                figure=figure,subplot=[0.55,0.1,0.4,0.7])
    img.set_tick_labels_font(size='small')
    img.set_axis_labels_font(size='small')
    vmax                     = max(pyfits.open(modelmapfile)[0].data.flatten())
    img.show_colorscale(cmap='gist_heat',aspect='auto',stretch='log',vmax=0.5) 
    
    for src in detectedSources:
      img.add_label(float(src.ra),float(src.dec),
                    "%s\n(ts = %i)" %(src.name,int(math.ceil(src.TS))),
                    relative=False,weight=weight,
                    color='green', size='small')
    pass
    
    # Modify the tick labels for precision and format
    img.tick_labels.set_xformat('ddd.dd')
    img.tick_labels.set_yformat('ddd.dd')
    
    # Display a grid and tweak the properties
    img.show_grid()
    
    ax                        = figure.add_axes([0.1,0.72,0.85,0.25],frame_on=False)
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False)
    col_labels                =['Source Name','TS','Energy flux','Photon index']
    table_vals                = map(lambda x:[x.name,"%i" %(int(math.ceil(x.TS))),
                                              "%s +/- %s" % (x.flux,x.fluxError),
                                              "%s +/- %s" % (x.photonIndex,x.photonIndexError)],detectedSources)
    
    if(len(table_vals)>0):
      the_table                 = ax.table(cellText=table_vals,
                                            colLabels=col_labels,
                                            loc='upper center')
    figure.canvas.draw()
    figure.savefig("likelihood_results.png")
  pass
  
  print(localizationMessage)
    
  return 'likexmlresults', outfilelike, 'TS', grb_TS, 'bestra', bestra, 'bestdec', bestdec, 'poserr', poserr, 'distance', distance
pass

thisCommand.run = run

if __name__=='__main__':
  #Get all key=value pairs as a dictionary
  args                           = dict(arg.split('=') for arg in sys.argv[1:])
  gtdolike(**args)
pass
