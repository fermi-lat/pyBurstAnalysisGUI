#### gtburst ####
install(DIRECTORY data/ 
        DESTINATION ${FERMI_INSTALL_DATADIR}/pyBurstAnalysisGUI)
install(DIRECTORY python/
        DESTINATION ${FERMI_INSTALL_PYTHON}
        PATTERN "*" EXCLUDE
        PATTERN "GtBurst/gtapps_mp/*"
        )
install(DIRECTORY python/ 
        DESTINATION ${FERMI_INSTALL_PYTHON}
        PATTERN "GtBurst/gtapps_mp/*" EXCLUDE
        PATTERN "README"
        PERMISSIONS OWNER_EXECUTE GROUP_EXECUTE 
        )
