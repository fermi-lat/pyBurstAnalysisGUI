#!/bin/bash
if [ -z "$CONDA_PREFIX" ]; then
    echo "Please run this script in a conda environment."
    exit 1
fi
PVER=$(python3 -c "import sys; print('{}.{}'.format(sys.version_info[0], sys.version_info[1]))")
echo "Python version: $PVER"
echo "Conda prefix: $CONDA_PREFIX"
echo "Installing GtBurst in conda environment..."

diff -t ${CONDA_PREFIX}/lib/python${PVER}/site-packages/fermitools/GtBurst python/GtBurst
diff -t ${CONDA_PREFIX}/lib/python${PVER}/site-packages/fermitools/gtburst.py python/gtburst.py
rm -rf ${CONDA_PREFIX}/lib/python${PVER}/site-packages/fermitools/GtBurst
rm ${CONDA_PREFIX}/lib/python${PVER}/site-packages/fermitools/gtburst.py
cp -r python/* ${CONDA_PREFIX}/lib/python${PVER}/site-packages/fermitools/.
