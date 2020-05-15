import os, errno
from GtBurst import commandDefiner
import multiprocessing

packageName                   = 'pyBurstAnalysisGUI'


class Configuration(object):

    def __init__(self):
        
      self.configuration              = {}
      
      # Data directory
      if os.environ.get("GTBURST_DATA") is not None:        
          
          self.configuration['dataRepository'] = os.environ.get("GTBURST_DATA")
      
      else:
          
          self.configuration['dataRepository'] = os.path.join(os.path.expanduser('~'),'FermiData')
      
      # FTP website
      
      if os.environ.get("GTBURST_FTP") is not None:        
          
          self.configuration['ftpWebsite']    = os.environ.get("GTBURST_FTP")
      
      else:
          
          self.configuration['ftpWebsite']    = "ftps://legacy.gsfc.nasa.gov/fermi/data"
      
      # Number of CPUs to use
      
      if os.environ.get("GTBURST_NCPUS") is not None:        
          
          self.configuration['maxNumberOfCPUs'] = os.environ.get("GTBURST_NCPUS")
      
      else:
          
          self.configuration['maxNumberOfCPUs'] = multiprocessing.cpu_count()
      
      if not os.path.exists(self.configuration['dataRepository']):
        
        #Create it!
        os.makedirs(os.path.abspath(self.configuration['dataRepository']))
            
      self.description = {}
      self.description['dataRepository'] = 'Directory for the storage of data'
      self.description['ftpWebsite']     = 'FTP data repository'
      self.description['maxNumberOfCPUs']= 'Max. number of CPUs to use'
    
    def set(self,key,value):

      if(key in list(self.description.keys())):

        self.configuration[key] = value

      else:

        raise RuntimeError("Got a unknown key!")
    
    def get(self,key):
    
      if key in list(self.description.keys()):
      
        return self.configuration[key]
        
      else:
      
        raise RuntimeError("Got a unknown key!")
    
    def getDescription(self,key):
    
      if(key in list(self.description.keys())):
    
        return self.description[key]

      
    def keys(self):

      return list(self.configuration.keys())

