# -*- python -*-
#
# $Id: SConscript,v 1.0 2013/01/25 12:00:12 giacomov Exp $
# Authors: Giacomo Vianello <giacomov@slac.stanford.edu>
# Version: pyBurstAnalysisGUI-01-00-00

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()

gtburstBin = progEnv.Program('gtburst', 'src/gtburst.cxx')

progEnv.Tool('registerTargets', package = 'pyBurstAnalysisGUI', 
             data = listFiles(['data/*']),
             python = listFiles(['python/*.py'])+listFiles(['python/GtBurst/*.py']),
             binaryCxts=[[gtburstBin, progEnv]])

