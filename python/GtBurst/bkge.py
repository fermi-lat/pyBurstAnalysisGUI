#author Vlasios Vasileiou 2009--

import sys,os,glob,math
import pyfits


import contextlib
import sys

systematicError = 0.15 # 15%

class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in 
    Python, i.e. will suppress all print, even if the print originates in a 
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).      

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])
#Verify if there exist the BKGE library
with suppress_stdout_stderr():
  try:
    import pyIrfLoader
    import ROOT
    ROOT.gStyle.SetOptStat(0)
    ROOT.gSystem.Load('librootIrfLoader')
    ROOT.gSystem.Load('libBackgroundEstimator')
    ROOT.gSystem.Load('libDurationEstimator')
    from ROOT import TOOLS, BackgroundEstimator
    from ROOT import BKGE_NS
    TOOLS.Set("GRB_NAME","") #this is to create the datadir directory of the bkg estimator (kluge)
    TOOLS.Set("BASEDIR",os.environ['BASEDIR'])
    active = True
  except:
    active = False
  pass
  sys.stderr.flush()
  sys.stdout.flush()
pass

class BKGE(object):
  def __init__(self,filteredFt1,ft2,triggerName,triggertime,outdir,inputDirectory,locError=0.0):
    global active
    if(active):
      self.ft1                = filteredFt1
      self.ft2                = ft2
      
      #Read the cuts from the filtered ft1 file
      f                       = pyfits.open(self.ft1)
      h                       = f[0].header
      self.ra                 = float(h['_ROI_RA'])
      self.dec                = float(h['_ROI_DEC'])
      self.roi                = float(h['_ROI_RAD'])
      self.emin               = float(h['_EMIN'])
      self.emax               = float(h['_EMAX'])
      self.zmax               = float(h['_ZMAX'])
      self.irf                = str(h['_IRF'])
      f.close()
      
      self.triggertime        = triggertime
      
      self.grbname            = triggerName
      self.outdir             = outdir
      self.locError           = locError
      
      TOOLS.Set("GRB_NAME",self.grbname)
      TOOLS.Set("GRB_RA",float(self.ra))
      TOOLS.Set("GRB_DEC",float(self.dec))
      TOOLS.Set("DATA_CLASS",self.irf) #do not move this down in initialize part
      TOOLS.Set("OUTPUT_DIR",self.outdir)
      TOOLS.Set("ROI_LOCALIZATION_ERROR",self.locError)
      
      TOOLS.Set("CALCULATE_ROI",1)
      TOOLS.Set("SEPARATE_FB",0)
      TOOLS.Set("BKG_ESTIMATE_ERROR",systematicError)
      TOOLS.Set("ROI_CONTAINMENT",0.95)
      TOOLS.Set("INPUT_DIR",inputDirectory)
      TOOLS.Set("ROI_MAX_RADIUS",12)
    else:
      raise RuntimeError("BKGE is not available in your system, you cannot use it!")
  pass
  
  def makeLikelihoodTemplate(self,start,stop,chatter=1):
    #Transform in MET
    start                     = float(start)+int(float(start) < 231292801.000)*self.triggertime
    stop                      = float(stop)+int(float(stop) < 231292801.000)*self.triggertime
    
    duration                  = stop-start
    self.CalculateBackground(start,stop,chatter)
    CR_EGAL_BKG               = ROOT.Double(0)
    GALGAMMAS_BKG             = ROOT.Double(0)
    
    outdir                    = self.outdir+'/Bkg_Estimates/%.2f_%.2f/' %(start-self.triggertime,stop-self.triggertime)
    
    BKGE_NS.MakeGtLikeTemplate(self.roi, outdir, self.irf, 
                               self.zmax, GALGAMMAS_BKG,  CR_EGAL_BKG,chatter)
    bkgetemplate              = os.path.abspath(os.path.join(outdir,'%s_gtlike_%s_CR_EGAL.txt' %(self.irf,self.grbname)))
    statisticalError          = 1.0/math.sqrt(CR_EGAL_BKG)
    return bkgetemplate, statisticalError, systematicError
  pass
  
  def CalculateBackground(self,start,stop,emin=50,emax=300000,ebins=20,chatter=1):
    #Transform in MET
    start                     = float(start)+int(float(start) < 231292801.000)*self.triggertime
    stop                      = float(stop)+int(float(stop) < 231292801.000)*self.triggertime

    duration                  = stop - start
    myname                    = sys._getframe().f_code.co_name
        
    if chatter>2: 
      TOOLS.PrintConfig()
    if chatter>1: 
      print "%s: Calculating Background..." %(myname)
    
    try:
      print("\nCalculating background between %.2f and %.2f from the trigger time %s...\n" %(start-self.triggertime,stop-self.triggertime,self.triggertime))
      BKGE_NS.CalculateBackground("%.2f_%.2f" %(start-self.triggertime,stop-self.triggertime), 
                                start, 
                                duration-1e-3,
                                str(self.ft1),
                                str(self.ft2),
                                str(self.irf),
                                emin,
                                emax,
                                ebins,
                                self.zmax,
                                chatter,
                                True)    
    except:
      print("BKGE failed!")
      raise
    pass
  pass
pass
