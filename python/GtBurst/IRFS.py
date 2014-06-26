import collections

class IRF(object):
  def __init__(self,shortname,name,reprocVer,evclass,galacticTemplate='',isotropicTemplate=''):
    self.shortname            = shortname
    self.name                 = name
    self.reprocessingVersion  = reprocVer
    self.evclass              = evclass
    self.galacticTemplate     = galacticTemplate
    self.isotropicTemplate    = isotropicTemplate
  pass
  
  def validateReprocessing(self,reproc):
    if(str(reproc) in self.reprocessingVersion.split(",")):
      return True
    else:
      return False
    pass
  pass
  
pass

class CaseInsensitiveDict(collections.OrderedDict):
  def __setitem__(self, key, value):
      super(CaseInsensitiveDict, self).__setitem__(key.lower(), value)
  pass
  
  def __getitem__(self, key):
      try:
        item                    = super(CaseInsensitiveDict, self).__getitem__(key.lower())
      except:
        #Try with the short name
        item                    = filter(lambda x:x.name==key,self.values())[0]
      pass
      return item
  pass
pass

def isOfClass(cl,evclass):
  bitmask                     = 1 << evclass
  if(cl & bitmask != 0):
    return True
  else:
    return False
  pass
pass


def fromEvclassToIRF(rep,event_class):
  #This assumes that classes are one into the other
  for i in range(1,6):
    
    if(i==1 and str(rep)=='202'):
      #In P7, evclass=1 is random (and must not be used, nor tested)
      continue
    pass
    
    if(isOfClass(event_class,i)==False):
      break
    pass
  pass
  
  evclass                     = i-1
  if(evclass==1 and str(rep)=='202'):
    evclass                   = 0
  pass
  
  #print("Evclass is %s for event_class %s" %(evclass,event_class))
  irfsForThisRepr             = PROCS[rep]
  for irf in irfsForThisRepr:
    if(IRFS[irf].evclass==evclass):
      return irf
    pass
  pass
pass


IRFS                          = CaseInsensitiveDict()

#P7REP
IRFS['P7REP_TRANSIENT']   = IRF('P7REP_TRANSIENT','P7REP_TRANSIENT_V15','202,203',0,'gll_iem_v05_rev1.fit,gll_iem_v05_rev1.fits,gll_iem_v05.fits,gll_iem_v05.fit,template_4years_P7_v15_repro_v2_trim.fits','iso_transient_v05.txt')
IRFS['P7REP_SOURCE']      = IRF('P7REP_SOURCE','P7REP_SOURCE_V15','202,203',2,'gll_iem_v05_rev1.fit,gll_iem_v05_rev1.fits,gll_iem_v05.fits,gll_iem_v05.fit,template_4years_P7_v15_repro_v2_trim.fits','iso_source_v05.txt,iso_source_v05_rev1.txt')
IRFS['P7REP_CLEAN']       = IRF('P7REP_CLEAN','P7REP_CLEAN_V15','202,203',3,'gll_iem_v05_rev1.fit,gll_iem_v05_rev1.fits,gll_iem_v05.fits,gll_iem_v05.fit,template_4years_P7_v15_repro_v2_trim.fits','iso_clean_v05.txt')
IRFS['P7REP_ULTRACLEAN']  = IRF('P7REP_ULTRACLEAN','P7REP_ULTRACLEAN_V15','202,203',4,'gll_iem_v05_rev1.fit,gll_iem_v05_rev1.fits,gll_iem_v05.fits,gll_iem_v05.fit,template_4years_P7_v15_repro_v2_trim.fits','iso_clean_v05.txt')


#P8        
IRFS['P8_TRANSIENT_R100'] = IRF('P8_TRANSIENT_R100','P8_TRANSIENT_R100_V1','300',1,'gll_iem_v05.fits','p8_p300x_isotropic_source_v1.txt')
IRFS['P8_TRANSIENT_R050'] = IRF('P8_TRANSIENT_R050','P8_TRANSIENT_R050_V1','300',2,'gll_iem_v05.fits','p8_p300x_isotropic_source_v1.txt')
IRFS['P8_TRANSIENT_R020'] = IRF('P8_TRANSIENT_R020','P8_TRANSIENT_R020_V1','300',3,'gll_iem_v05.fits','p8_p300x_isotropic_source_v1.txt')
IRFS['P8_SOURCE']         = IRF('P8_SOURCE','P8_SOURCE_V1','300',4,'gll_iem_v05.fits','p8_p300x_isotropic_source_v1.txt')                

PROCS                         = collections.OrderedDict()
for k,v in IRFS.iteritems():
  thisReprocessings       = v.reprocessingVersion.split(",")
  for repro in thisReprocessings:
    if(repro in PROCS.keys()):
      PROCS[repro].append(k)
    else:
      PROCS[repro] = [k]
    pass
  pass
pass
  
