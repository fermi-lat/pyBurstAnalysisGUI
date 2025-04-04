#!/bin/bash
PYTHON=3.11
diff -t ${CONDA_PREFIX}/lib/python${PYTHON}/site-packages/fermitools/GtBurst python/GtBurst
diff -t ${CONDA_PREFIX}/lib/python${PYTHON}/site-packages/fermitools/gtburst.py python/gtburst.py
rm -rf ${CONDA_PREFIX}/lib/python${PYTHON}/site-packages/fermitools/GtBurst
rm ${CONDA_PREFIX}/lib/python${PYTHON}/site-packages/fermitools/gtburst.py
cp -r python/* ${CONDA_PREFIX}/lib/python${PYTHON}/site-packages/fermitools/.
