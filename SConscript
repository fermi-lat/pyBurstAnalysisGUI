# -*- python -*-
#
# $Id: SConscript,v 1.18 2014/05/23 01:08:16 giacomov Exp $
# Authors: Giacomo Vianello <giacomov@slac.stanford.edu>
# Version: pyBurstAnalysisGUI-02-00-00

Import('baseEnv')
Import('listFiles')
progEnv = baseEnv.Clone()

gtburstBin = progEnv.Program('gtburst', 'src/gtburst.cxx')

progEnv.Tool('registerTargets', package = 'pyBurstAnalysisGUI', 
             data = (listFiles(['data/*']) +
                     listFiles(['data/templates/*']) +
                     listFiles(['data/tcl_extensions/*']) +
                     listFiles(['data/tcl_extensions/fsdialog/*']) +
                     listFiles(['data/tcl_extensions/msgcat/*'])),
             python = (listFiles(['python/*.py']) +
                       listFiles(['python/GtBurst/*.py']) +
                       listFiles(['python/GtBurst/aplpy/*.py']) +
                       listFiles(['python/GtBurst/commands/*.py']) +
                       listFiles(['python/GtBurst/gtapps_mp/*.py']) +
                       listFiles(['python/GtBurst/scripts/*.py']) +
                       listFiles(['python/GtBurst/tasks/*.py'])),
             binaryCxts=[[gtburstBin, progEnv]])
