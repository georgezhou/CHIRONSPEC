import os,sys,string,pyfits,pickle,glob
from numpy import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.pyplot.switch_backend('agg')

import astropy.time
from PyAstronomy import pyasl


# Convert HH:MM:SS.SSS into Degrees :
def convHMS(ra):
   try :
      sep1 = ra.find(':')
      hh=int(ra[0:sep1])
      sep2 = ra[sep1+1:].find(':')
      mm=int(ra[sep1+1:sep1+sep2+1])
      ss=float(ra[sep1+sep2+2:])
   except:
      raise
   else:
      pass
   
   return(hh*15.+mm/4.+ss/240.)

# Convert +DD:MM:SS.SSS into Degrees :
def convDMS(dec):

   Csign=dec[0]
   if Csign=='-':
      sign=-1.
      off = 1
   elif Csign=='+':
      sign= 1.
      off = 1
   else:
      sign= 1.
      off = 0

   try :
      sep1 = dec.find(':')
      deg=int(dec[off:sep1])
      sep2 = dec[sep1+1:].find(':')
      arcmin=int(dec[sep1+1:sep1+sep2+1])
      arcsec=float(dec[sep1+sep2+2:])
   except:
      raise
   else:
      pass

   return(sign*(deg+(arcmin*5./3.+arcsec*5./180.)/100.))

def bjdcalc(hjd,ra,dec):

    inlc = array([[hjd,0,0]])
    fname = str(random.uniform(0,1.))[2:]
    savetxt(fname,inlc,fmt="%.10f")
    

        
    ### Use this for everything else
    command = "vartools -i "+fname+" -converttime input hjd inputsys-utc output bjd outputsys-tdb radec fix "+str(ra)+" "+str(dec)+" -o "+fname+".o"
    os.system(command)

    olc = loadtxt(fname+".o")
    bjd = olc[0]
    os.system("rm "+fname+"*")
    return bjd

def setjd(header,ra,dec):

    try:
       dateobs = header["EMMNWB"]
       t = astropy.time.Time(dateobs, format='fits', scale='utc')
    except ValueError:
       dateobs = header["DATE-OBS"]
       t = astropy.time.Time(dateobs, format='fits', scale='utc')
    

    try:
       jd = float(header["JD"])
    except KeyError:
       jd = t.jd

    # Coordinates of European Southern Observatory
    # (Coordinates of UT1)
    latitude = -30.169661
    longitude = -70.806525
    altitude = 2207.
    

    # Coordinates of HD 12345 (J2000)
    ra2000 = convHMS(ra)
    dec2000 = convDMS(dec)

    # Calculate barycentric correction (debug=True show
    # various intermediate results)
    corr, hjd = pyasl.helcorr(longitude, latitude, altitude, \
                              ra2000, dec2000, jd, debug=False)


    # #if abs(float(header["JD"])-jd) > 0.003:
    #    print "ERROR, check JD of exposures! mismatch between DATE-OBS or EMMNWB and JD header"
    #    #sys.exit(1)
    #    corr = header["BCV"]

    #    strout = header["OBJECT"] + " " + header["DATE-OBS"]
    #    strout = open("badfits_JD").read() + strout+"\n"
    #    o = open("badfits_JD","w")
    #    o.write(strout)
    #    o.close()
       

    bjd = bjdcalc(hjd,ra2000,dec2000)
    print corr,jd,hjd,bjd

    return corr,jd,hjd,bjd
    

def setjd_fits(fits):
    fitsdata = pyfits.open(fits)
    header = fitsdata[0].header
    ra = header["RA"]
    dec = header["DEC"]
    corr,jd,hjd,bjd = setjd(header,ra,dec)

    fitsdata[0].header["BCV"] = str(corr)
    fitsdata[0].header["JD"] = str(jd)
    fitsdata[0].header["HJD"] = str(hjd)
    fitsdata[0].header["BJD"] = str(bjd)
    del fitsdata[0].header["GAINS"]

    fitsdata.writeto(fits,clobber=True)

    return corr,jd,hjd,bjd


def setjd_folder(folder):

    fitslist = sort(glob.glob(folder+"achi*.fits"))
    for fits in fitslist:
       print fits
       header = pyfits.getheader(fits)
       try:
          bjd = header["BJD"]
       except KeyError:
          try:
             setjd_fits(fits)
          except:
             print "removing ",fits
             os.system("mv "+fits+" /data/tfop/chiron_data/rmed/")

   


if __name__ == "__main__":
    #fits = "/media/Onion/Data/Chiron/171019_planid_418/achi171019.1118.fits"
    #header = pyfits.getheader(fits)
    #setjd(header)

    folder = "/media/Onion/Data/Chiron/190222_planid_500/"
    setjd_folder(folder)
