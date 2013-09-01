
codes = {}
codes[1] = 'No network or generic network failure'
codes[11] = 'LAT server time out'
codes[12] = 'LAT data server is probably down'
codes[13] = 'Corrupted downloaded file'
codes[14] = 'No data coverage'
codes[2] = 'No events survived the data cuts'
codes[21] = 'Too few counts for the likelihood analysis'
codes[22] = "gtmktime failed"
codes[23] = "gtselect failed"
codes[24] = "gtbin failed during skymap production"
codes[25] = "gtbin failed during skycube production"
codes[26] = "gtltcube failed"
codes[27] = "gtexpmap failed"
codes[28] = "gtexpcube2 failed"
codes[29] = "gtsrcmaps failed"
codes[201] = "gtmodel failed"
codes[202] = "gtdiffrsp failed"
codes[203] = "gtrspgen failed"
codes[204] = "gtbin failed during PHA1 production"
codes[205] = "gtbkg failed"
codes[206] = "gttsmap failed"
codes[207] = "gtfindsrc failed"

codes[3] = "I/O error"
codes[31] = "Error opening filtered event file"
codes[4] = "Could not update"

class GtBurstException(RuntimeError):
  def __init__(self,code,message):
    RuntimeError.__init__(self,message)
    self.shortMessage         = codes[code]
    self.longMessage          = message
    self.code                 = code
  pass
###
