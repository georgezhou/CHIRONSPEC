import os,sys,pickle
from numpy import *
import matplotlib.pyplot as plt
import matplotlib

from sklearn.ensemble import GradientBoostingRegressor,GradientBoostingClassifier
from sklearn.multioutput import MultiOutputRegressor


def binccf(ccf):
    binwidth = 1.
    vel = arange(-300,300,binwidth)
    ccfbin = zeros(len(vel))
    for i in range(len(vel)):
        mask = abs(ccf[:,0]-vel[i])<binwidth/2
        if len(ccf[:,1][mask]) > 0:
            ccfbin[i] = nanmean(ccf[:,1][mask])
        else:
            ccfbin[i] = 0

    mask = ccfbin != ccfbin
    ccfbin[mask] = 0

    return ccfbin

def checksb2(spectrum):
    ccf_list,ccf,vsini_init,rv = pickle.load(open(os.path.dirname(spectrum)+"/lsd_"+os.path.basename(spectrum)+".pkl","rb"))
    ccf[:,0] -= rv    
    ccf = binccf(ccf)
    
    sb2 = pickle.load(open("SPEC_ML/sb2.pkl","rb"))
    proba = sb2.predict_proba([ccf])
    print spectrum
    print proba
    
    #if proba[0][1] > 0.6:
        #plt.clf()
        #plt.plot(ccf)
        #plt.show()

    return proba[0][1]


def do_all():
    import pandas
    database = pandas.read_csv("database.csv")

    #database["SB2_prob"] = nan*zeros(len(database))
    
    for i in range(len(database)):

        
        fitsfile = database['filepath'].iloc[i]

        prob = checksb2(fitsfile)

        database["SB2_prob"].loc[i] = prob
        #@print database.loc[i]

    #print database

    database.to_csv("database_sb2.csv",index=False)
    
if __name__ == "__main__":
    
    do_all()
