#### gtburst ####
install(DIRECTORY data/ 
        DESTINATION ${FERMI_INSTALL_DATADIR}/pyBurstAnalysisGUI)
install(DIRECTORY python/ 
        DESTINATION ${FERMI_INSTALL_PYTHON}
        PATTERN "gtapps_mp/*" EXCLUDE
        PATTERN "README"
        PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE 
        )
