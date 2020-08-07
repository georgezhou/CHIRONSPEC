import os,sys,glob,string,pyfits
from numpy import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.pyplot.switch_backend('agg')
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

import setjd
import cc_rv
import lsd
import plot_spectra_anl
import spectype

import run_specml_classifyfirst as run_specml
#from astroquery.vizier import Vizier
from astroquery.mast import Catalogs
from astroquery.simbad import Simbad

import warnings
warnings.filterwarnings("ignore")

speclib = "spectral_library/"


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





def roundnearest(val,r):
    val = val / r
    val = round(val)
    val = val*r
    return val



def querytic(fitsheader):

    tic = fitsheader["OBJECT"]

    ### TRY the T0XXX format
    if tic[0] == "T":
       try:
          tic = int(tic[1:])
          tic = "TIC"+str(tic)

          print tic

          #catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="Tic",version=8,objType="STAR",Tmag=[-5,18])
          catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="TIC")

          indx = argsort(catalog_data["Tmag"])
          star = catalog_data[indx[0]]
          ra = deg2HMS(star["ra"])
          dec = deg2DMS(star["dec"])

          print "ticquery NEW RA DEC",ra,dec,"old RA DEC",fitsheader["RA"],fitsheader["DEC"]

       except:

          ra = fitsheader["RA"]
          dec = fitsheader["DEC"]


    ### Try the no "TIC" string format that Avi uses
    elif is_number(tic):
       try:
          tic = int(tic)
          tic = "TIC"+str(tic)

          print tic

          #catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="Tic",version=8,objType="STAR",Tmag=[-5,18])
          catalog_data = Catalogs.query_object(tic, radius=0.0002,catalog="TIC")

          indx = argsort(catalog_data["Tmag"])
          star = catalog_data[indx[0]]
          ra = deg2HMS(star["ra"])
          dec = deg2DMS(star["dec"])

          print "ticquery NEW RA DEC",ra,dec,"old RA DEC",fitsheader["RA"],fitsheader["DEC"]

       except:

          ra = fitsheader["RA"]
          dec = fitsheader["DEC"]
  


    ### Try a simbad name query, works for the HD stuff
          
    else:

       try:
          catalog_data = Simbad.query_object(tic)
          star = catalog_data[0]
          ra = string.replace(str(star["RA"])," ",":")
          dec = string.replace(str(star["DEC"])," ",":")

          print "ticquery NEW RA DEC",ra,dec,"old RA DEC",fitsheader["RA"],fitsheader["DEC"]

       except:

          ra = fitsheader["RA"]
          dec = fitsheader["DEC"]
     
    return ra,dec



 

# def getteff(c):
#     gaia = Vizier(columns=["*"], catalog="I/345")
#     radius = 5
#     found = False
#     while not found:
#        results = gaia.query_region(c, radius=str(radius)+"s")
#        if len(results) > 0:
#           if len(results[0]) > 0:
#               print "Nstars found",len(results[0])
#               gmag = []
#               for i in range(len(results[0])):
#                   gmag.append(results[0][i]["Gmag"])

#               if min(gmag) < 15:
#                  print gmag

#                  gmag = array(gmag)
#                  indx = argmin(gmag)

#                  results = results[0][indx]

#                  teff = results["Teff"]
#                  rstar = results["Rad"]
#                  gmag = results["Gmag"]
#                  found = True
#               else:
#                  radius += 5
#           else:
#               radius += 5
                 

#        else:
#           radius += 5
              
#     print teff,rstar,gmag 
#     return teff,rstar,gmag






def getteff(c,cone=10):

    G = nan
    while G != G and cone < 180:
       print "cone search radius",cone
       #if True:
       try:
           width = u.Quantity(cone, u.arcsec)
           height = u.Quantity(cone, u.arcsec)
           result = Gaia.cone_search_async(c,width)
           result =  result.get_results()
           #print result.columns

           useindx = 0
           print "len result",len(result)
           if len(result) > 1:
               Gmag = []
               for i in range(len(result)):
                   Gmag.append(result[i]["phot_g_mean_mag"])
               Gmag = array(Gmag)
               print Gmag
               useindx = argmin(Gmag)

           G =  result[useindx]["phot_g_mean_mag"]
           teff =  result[useindx]["teff_val"]
           rstar = result[useindx]["radius_val"]

           if not is_number(G):
               G = nan
               teff = nan
               rstar = nan

           if G > 15:
               G = nan
               teff = nan
               rstar = nan
           
       except IndexError:
           G = nan
           teff = nan
           rstar = nan

       cone += 20
        
    return teff,rstar,G





 
def writerv(rvlist,delim=","):
    o = ""
    for i in rvlist:
        o += str(i)+delim

    return o

def measure_snr(fitsname):
   fits = pyfits.open(fitsname)
   data = fits[0].data
   snrlist = []
   for f in data:
      order = f[:,1]      
      snr_order = nanmedian(order[len(order)/2-100:len(order)/2+100])
      snr_order = snr_order/sqrt(snr_order)
      snr_order *= 2.5
      snrlist.append(snr_order)

   return nanmax(array(snrlist))

 
def main(folder):
    fitslist = sort(glob.glob(folder+"/*.fits"))
    for fits in fitslist:
        if not os.path.exists(fits+".rv"):

           fitsheader = pyfits.getheader(fits)

           if not fitsheader["OBJECT"] == "ThAr" and not os.path.exists(fits+".rv"):
              print fits
              objectname = fitsheader["OBJECT"]
              ra,dec = querytic(fitsheader) #fitsheader["RA"],fitsheader["DEC"]
              exptime = fitsheader["EXPTIME"]
              #bcorr,jd,hjd,bjd = setjd.setjd(fitsheader)
              #bcorr,jd,hjd,bjd = float(fitsheader["BCV"]),float(fitsheader["JD"]),float(fitsheader["HJD"]),float(fitsheader["BJD"])

              bcorr,dummy,dummy,dummy = setjd.setjd(fitsheader,ra,dec)
              dummy,jd,hjd,bjd = float(fitsheader["BCV"]),float(fitsheader["JD"]),float(fitsheader["HJD"]),float(fitsheader["BJD"])
              
              print ra,dec,bcorr,float(fitsheader["BCV"])

              c = SkyCoord(ra+" "+dec, frame='icrs', unit=(u.hourangle, u.deg))

              try:
              #if True:
                  teff,rstar,gmag = getteff(c)
                  print "Gaia teff,rstar",teff,rstar
                  teff = int(roundnearest(teff,250))
              except:
                  print "Error on teff"
                  teff = 6000
                  rstar = 1.
                  gmag = -99.


              if teff < 4000:
                  teff = 4000
              if teff > 10000:
                  teff = 10000

              if rstar < 3:
                 logg = 3.5
              else:
                 logg = 2.0 ### initial logg minimum. If rstar < 3, allow all logg, else only allow dwarfish logg


              template = speclib+"template_"+str(teff)+"_4.0_0.0.dat"

              snr = measure_snr(fits)
              print fits
              print template


              ccf,vsini,lsdshift = lsd.run_spectrum(fits,template)
              
              if teff < 10000:
                 teff,logg,feh,vsini,dummy = run_specml.runML(fits)
                 vsini,lsdshift,rvlist = spectype.run_rv(fits,folder+"/lsd_"+os.path.basename(fits)+".pkl")

              else:
                 teff,logg,feh,vsini,lsdshift,rvlist = spectype.get_best_template(fits,folder+"/lsd_"+os.path.basename(fits)+".pkl",teff,logg)

                 
              template = speclib+"template_"+str(int(roundnearest(teff,250)))+"_"+str(roundnearest(logg,0.5))+"_"+str(roundnearest(feh,0.5))+".dat"
              if not os.path.exists(template):
                 print "Template doesn't exist, trying a logg=4.5,feh=0.0 template"
                 #feh,logg = 0.0,4.5
                 #template = speclib+"template_"+str(int(teff))+"_"+str(logg)+"_"+str(feh)+".dat"
                 teffround = int(roundnearest(teff,250))
                 if teffround < 4000:
                    teffround = 4000
                 if teffround > 10000:
                    teffround = 10000
                 template = speclib+"template_"+str(teffround)+"_"+str(4.5)+"_"+str(0.0)+".dat"

                 
              print "using template",template


              #rvlist,ccfrv = cc_rv.main(fits,template,vsini=vsini)
              rvlist[:,1] += bcorr
              ccfrv = 0#+= bcorr
              
              telrv = cc_rv.measure_telluric_rv(fits,template,vsini,lsdshift)
              #telrv = 0

              print "CCF RV",ccfrv
              print "diagnostic",fits,template,int(teff),logg,feh,vsini,lsdshift

              plot_spectra_anl.main(fits,template,teff,logg,feh,vsini,lsdshift)
              rvlist = [objectname,fits,bjd,convHMS(ra),convDMS(dec),gmag,exptime,teff,logg,feh,snr,bcorr,telrv,vsini,lsdshift+bcorr,ccfrv]+list(rvlist[:,1])
              rvlist = writerv(rvlist)+"\n"
              out = open(fits+".rv","w")
              out.write(rvlist)
              out.close()


              plt.clf()
        
if __name__ == "__main__":

    folder = "/data/tfop/chiron_data/reduced/"
    main(folder)
