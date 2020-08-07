import os,sys,glob,string,pyfits,pandas
from numpy import *
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
import setjd
from astroquery.mast import Catalogs
from astropy.time import Time

def is_number(s):
   try:
      float(s)
      return True
   except ValueError:
      return False



# Convert RA (deg) to H.M.S:
def deg2HMS( RAin ):

   if(RAin<0):
      sign = -1
      ra   = -RAin
   else:
      sign = 1
      ra   = RAin

   h = int( ra/15. )
   ra -= h*15.
   m = int( ra*4.)
   ra -= m/4.
   s = ra*240.

   if(sign == -1):
      out = '-%02d:%02d:%06.3f'%(h,m,s)
   else: out = '+%02d:%02d:%06.3f'%(h,m,s)
   
   return out
   
# Convert Decl. (deg) to D.M.S:
def deg2DMS( Decin ):

   if(Decin<0):
      sign = -1
      dec  = -Decin
   else:
      sign = 1
      dec  = Decin

   d = int( dec )
   dec -= d
   dec *= 100.
   m = int( dec*3./5. )
   dec -= m*5./3.
   s = dec*180./5.

   if(sign == -1):
      out = '-%02d:%02d:%06.3f'%(d,m,s)
   else: out = '+%02d:%02d:%06.3f'%(d,m,s)

   return out


     
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




def querytic(fitsheader):

    tic = fitsheader["OBJECT"]

    if tic[0] == "T" and tic[1] != "O" and tic[1] != "I":
       try:
          tic = int(tic[1:])
          tic = "TIC"+str(tic)

          print tic

          #catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="Tic",version=8,objType="STAR",Tmag=[-5,18])
          catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="TIC")

          indx = argsort(catalog_data["Tmag"])
          star = catalog_data[indx[0]]

          ra = star["ra"]
          dec = star["dec"]
          print deg2HMS(ra),deg2DMS(dec)
          pmra = star["pmRA"]
          pmdec = star["pmDEC"]
          nyear = Time(fitsheader["DATE-OBS"]).jyear-2000

          ra = ra+cos(dec*pi/180.)*(pmra*0.001/(60*60))*nyear
          dec = dec+(pmdec*0.001/(60*60))*nyear

          
          ra = deg2HMS(ra)
          dec = deg2DMS(dec) 
          print ra,dec
         
          #ra = deg2HMS(star["ra"])
          #dec = deg2DMS(star["dec"])

          print "ticquery NEW RA DEC",ra,dec,"old RA DEC",fitsheader["RA"],fitsheader["DEC"]
       except:

          ra = fitsheader["RA"]
          dec = fitsheader["DEC"]


    elif is_number(tic):
       try:
          tic = int(tic)
          tic = "TIC"+str(tic)

          print tic

          #catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="Tic",version=8,objType="STAR",Tmag=[-5,18])
          catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="TIC")

          indx = argsort(catalog_data["Tmag"])
          star = catalog_data[indx[0]]
         

          ra = star["ra"]
          dec = star["dec"]
          print deg2HMS(ra),deg2DMS(dec)
          pmra = star["pmRA"]
          pmdec = star["pmDEC"]
          nyear = Time(fitsheader["DATE-OBS"]).jyear-2000

          ra = ra+cos(dec*pi/180.)*(pmra*0.001/(60*60))*nyear
          dec = dec+(pmdec*0.001/(60*60))*nyear

          
          ra = deg2HMS(ra)
          dec = deg2DMS(dec) 
          print ra,dec

          print "ticquery NEW RA DEC",ra,dec,"old RA DEC",fitsheader["RA"],fitsheader["DEC"]

       except:

          ra = fitsheader["RA"]
          dec = fitsheader["DEC"]
  

          
    else:

       # try:
       #    catalog_data = Simbad.query_object(tic)
       #    star = catalog_data[0]
       #    ra = string.replace(str(star["RA"])," ",":")
       #    dec = string.replace(str(star["DEC"])," ",":")

       #    print "ticquery NEW RA DEC",ra,dec,"old RA DEC",fitsheader["RA"],fitsheader["DEC"]

       # except:

       #    ra = fitsheader["RA"]
       #    dec = fitsheader["DEC"]



        
       ra = fitsheader["RA"]
       dec = fitsheader["DEC"]

       
    return ra,dec

          
if __name__ == "__main__":

    database = pandas.read_csv("database.csv")

    for i in range(len(database)):

        if os.path.exists(database.loc[i,"filepath"]):
            fitsheader = pyfits.getheader(database.loc[i,"filepath"])
            print database.loc[i,"filepath"]
            
            ra,dec = querytic(fitsheader)
            bcorr,dummy,dummy,dummy = setjd.setjd(fitsheader,ra,dec)

            print database.loc[i,"bcorr"],bcorr

            database.loc[i,"lsdRV"] = database.loc[i,"lsdRV"]-database.loc[i,"bcorr"]+bcorr

            database.loc[i,"bcorr"] = bcorr


    database.to_csv("database.csv",index=False)
            

    # database = pandas.read_csv("database.csv")

    # for i in range(len(database)):

    #     if os.path.exists(database.loc[i,"filepath"]):

    #         if sys.argv[1] in database.loc[i,"filepath"]:
    #            fitsheader = pyfits.getheader(database.loc[i,"filepath"])

    #            print database.loc[i,"filepath"]

    #            ra,dec = querytic(fitsheader)
    #            bcorr,dummy,dummy,dummy = setjd.setjd(fitsheader,ra,dec)

    #            print database.loc[i,"bcorr"],bcorr

    #            database.loc[i,"lsdRV"] = database.loc[i,"lsdRV"]-database.loc[i,"bcorr"]+bcorr

    #            database.loc[i,"bcorr"] = bcorr
