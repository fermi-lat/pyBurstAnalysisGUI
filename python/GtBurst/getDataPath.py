import GtBurst
import os

def getDataPath():
  path                 = GtBurst.__file__
  installationPath     = os.path.join(os.path.sep.join(path.split(os.path.sep)[0:-3]))
  dataPath             = os.path.join(installationPath,'data')
  if(not os.path.exists(os.path.join(dataPath,'glast_logo.png'))):
    #In the version within the Science Tools, data are saved in data/GtBurst/
    dataPath           = os.path.join(dataPath,'pyBurstAnalysisGUI')
  pass
  return dataPath
