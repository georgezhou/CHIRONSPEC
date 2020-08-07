import os,sys,glob,string,pandas
from numpy import *
import matplotlib.pyplot as plt
from scipy import optimize,interpolate
import pickle,emcee



def calc_bis(ccf,v0):

    ccf[:,1] /= max(ccf[:,1])
    vel = arange(min(ccf[:,0]),max(ccf[:,0]),0.001)
    ccf = interpolate.splrep(ccf[:,0],ccf[:,1],k=3)
    ccf = interpolate.splev(vel,ccf)
    mask_left = vel < v0
    mask_right = vel > v0

    def dolim(low_0,low_1,high_0,high_1):
        mask_low = (ccf > low_0) * (ccf < low_1)
        mask_high = (ccf > high_0) * (ccf < high_1)

        low = mean([mean(vel[mask_left*mask_low]),mean(vel[mask_right*mask_low])])
        high = mean([mean(vel[mask_left*mask_high]),mean(vel[mask_right*mask_high])])
        bis_out = (high-low) / (high_0-low_0)

        return bis_out
        

    bis = []
    bis.append(dolim(0.3,0.35,0.7,0.75))
    bis.append(dolim(0.35,0.4,0.75,0.8))
    bis.append(dolim(0.4,0.45,0.8,0.85))
    bis.append(dolim(0.45,0.5,0.85,0.9))
    biserr = std(bis)/sqrt(float(len(bis)))

    #plt.plot(vel,ccf)
    bis = []
    yval = arange(0.1,0.95,0.05)
    for y in yval:
        mask = (ccf > y-0.05) * (ccf < y+0.05)

        bis.append(mean([mean(vel[mask_left*mask]),mean(vel[mask_right*mask])]))
    bis = array(bis)
    plt.plot(bis-v0,yval)
    #plt.show()

    return dolim(0.3,0.4,0.7,0.8),biserr



def dofits(fits,v0):
    fits_path = os.path.dirname(fits)+"/"
    fitsname = os.path.basename(fits)
    ccf_list,ccf_master,vsini,shift = pickle.load(open(fits_path+"lsd_"+fitsname+".pkl","rb"))

    
    bis,biserr = calc_bis(ccf_master,v0)
    return bis,biserr

def dostar(star):

    database = "database.csv"
    database = pandas.read_csv(database)

    starname = star
    decker = "slicer"

    def mmd(x):
        return nanmedian(abs(x-nanmedian(x))) / sqrt(float(len(x)))

    star = database[database["objectname"] == starname]
    star = star[star["decker"] == decker]

    data = []
    for i in range(len(star)):
        ordervels = []
        for order in arange(1,35):
            ordervels.append(star["order"+str(order)].iloc[i])
        print ordervels
        rverr = mmd(ordervels)

        fits = star.iloc[i]["filepath"]
        bis,biserr = dofits(fits,star.iloc[i]["lsdRV"]-star.iloc[i]["bcorr"])
        
        data.append([star.iloc[i]["bjd"],star.iloc[i]["lsdRV"],rverr,bis,biserr])

    plt.show()

    data = array(data)


    plt.errorbar(data[:,1],data[:,3],xerr=data[:,2],yerr=data[:,4],linestyle="none")
    #plt.errorbar(data[:,0],data[:,3],yerr=data[:,4],linestyle="none")
    plt.show()

    
if __name__ == "__main__":
    dostar(sys.argv[1])
