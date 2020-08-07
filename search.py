import os,sys,string,pandas
from numpy import *
import matplotlib.pyplot as plt
import astropy.time



database = "database.csv"
database = pandas.read_csv(database)


def reportstar(starname):
    #print "#######################################################"
    print "\n"
    print starname
    star = database[database["objectname"] == starname]

    #    print len(star)
    for i in range(len(star)):

        bjdtime = astropy.time.Time(star["bjd"].iloc[i], format='jd', scale='tdb')
        dateobs = bjdtime.fits        
        report = ""
        report += "BJD="+str(star["bjd"].iloc[i])+" , "
        report += "DATEOBS="+str(dateobs)+" , "
        report += "RV="+str(round(star["lsdRV"].iloc[i],3))+" , "
        report += "Teff="+str(star["teff"].iloc[i])+" , "
        report += "logg="+str(star["logg"].iloc[i])+" , "
        report += "vsini="+str(round(star["vsini"].iloc[i],1))+" , "
        report += "[Fe/H]="+str(star["feh"].iloc[i])+" , "
        #report += "SNRe="+str(round(star["snr"].iloc[i],1))+" , "
        report += "Decker="+star["decker"].iloc[i]
        print report



def runreport(utstart,utend):

    
    utstart = astropy.time.Time(utstart, format='fits', scale='utc')
    utend = astropy.time.Time(utend, format="fits", scale="utc")
    jdstart = float(utstart.jd)
    jdend = float(utend.jd)

    print jdstart,jdend
    mask = database["bjd"]>jdstart
    mask *= database["bjd"]<jdend
    run = database[mask]
    obslist = sort(unique(run["objectname"]))

    print obslist

    for objectname in obslist:
        reportstar(objectname)

if __name__ == "__main__":
    #reportstar(sys.argv[1])

    utstart = "2018-11-01"
    utend = "2019-09-01"
    runreport(utstart,utend)
