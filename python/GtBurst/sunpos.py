import skymaps
import os

#This converts from MET to Julian Date
def _jd_from_MET( met ):
    jd=skymaps.JulianDate((skymaps.JulianDate_missionStart().seconds()
                           + met)/skymaps.JulianDate.secondsPerDay)
    return jd

def getSunPosition( met ):
    #environment variable (you need FTOOLS installed and configured)
    os.environ['TIMING_DIR']=os.path.join(os.environ['HEADAS'],"refdata")
    #Get the sun direction
    sun=skymaps.SolarSystem(skymaps.SolarSystem.SUN)
    SunSkyDir = sun.direction(_jd_from_MET(met))
    return SunSkyDir
