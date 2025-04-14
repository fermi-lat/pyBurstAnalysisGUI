#### gtburst ####
install(DIRECTORY data/ 
        DESTINATION ${FERMI_INSTALL_DATADIR}/pyBurstAnalysisGUI)
install(DIRECTORY python/ 
        DESTINATION ${FERMI_INSTALL_PYTHON}
        PATTERN "*.py"
        PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE 
        )
