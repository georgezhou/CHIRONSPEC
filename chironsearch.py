import os,sys,pandas,string
from numpy import *
import argparse
database = pandas.read_csv("/data/tfop/chiron_analysis/database.csv")
import warnings
warnings.filterwarnings('ignore')



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


def mmd(x):
    stdev = nanstd(x)
    mask = abs(x-nanmedian(x)) < 5*stdev
    x = x[mask]
    return 1.35*nanmedian(abs(x-nanmedian(x))) / sqrt(float(len(x)))



def makereport(star):

    reportout = ""
    for i in range(len(star)):
        filepath = str(star["filepath"].iloc[i])
        filepath = os.path.basename(filepath)
        filepath = string.split(filepath,"_")[1]



        ordervels = []
        for order in arange(1,35):
            ordervels.append(star["order"+str(order)].iloc[i])

        ordervels = array(ordervels)

        errorbar = mmd(ordervels)


        reportout += star["objectname"].iloc[i]+" "
        reportout += filepath+" "
        reportout += str(star["bjd"].iloc[i])+" "
        reportout += "RV  "+str(round(star["lsdRV"].iloc[i],3))+" "+str(round(errorbar,3))+" "
        reportout += "Teff "+str(int(star["teff"].iloc[i]))+" "
        reportout += "logg "+str(round(star["logg"].iloc[i],2))+" "
        reportout += "feh "+str(round(star["feh"].iloc[i],2))+" "
        reportout += "vsini "+str(round(star["vsini"].iloc[i],1))+" "
        reportout += "SB2 Probability "+str(round(star["SB2_prob"].iloc[i],1))+" "
        reportout += star["filepath"].iloc[i]+"\n"
    
    #return reportout
    print reportout


def searchobjname(objectname):
    
    star = database[database["objectname"] == objectname]
    return star

def searchcoord(ra,dec):
    mask = abs(database["ra"]-ra) < 20/60./60.
    mask *= abs(database["dec"]-dec) < 20/60./60.

    star = database[mask]
    return star
    


if __name__ == "__main__":
    
    # Initiate the parser
    parser = argparse.ArgumentParser()

    # Add long and short argument
    parser.add_argument("--coord", "-c", help="coordinates \"HH:MM:SS DD:MM:SS\" or \"D.DDD +-D.DDD \" decimal deg. Please use quotation marks!")
    # Add long and short argument
    parser.add_argument("--objectname", "-o", help="Object name search. Exact match needed")

    # Read arguments from the command line
    args = parser.parse_args()


    if args.objectname != None:        
        star = searchobjname(args.objectname)
        print "Making report for ",args.objectname

        makereport(star)
        
    elif args.coord != None:

        ra,dec = string.split(args.coord)
        
        if is_number(ra) and is_number(dec):
            star = searchcoord(float(ra),float(dec))

        else:
            ra = convHMS(ra)
            dec = convDMS(dec)
            star = searchcoord(ra,dec)

        print "Making report for ",ra,dec
        makereport(star)
        
        

    
    else:
        print "please use this search properly gesh!",args



