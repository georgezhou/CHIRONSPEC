import os,sys,string,glob,pyfits
from numpy import *
import matplotlib.pyplot as plt
import pandas
from astropy.time import Time
from datetime import datetime,timedelta
import is_sb2


basedir = "/data/tfop/chiron_analysis/"
def create_database():
    columns = ["objectname"]

    database = pandas.DataFrame(["TEST"],columns=columns)
    database['decker'] = pandas.Series(["TEST"])#pandas.Series(zeros(len(database)))
    database['filepath'] = pandas.Series(["TEST"])#pandas.Series(zeros(len(database)))
    database['bjd'] = pandas.Series(zeros(len(database)))
    database['ra'] = pandas.Series(zeros(len(database)))
    database['dec'] = pandas.Series(zeros(len(database)))
    database['gmag'] = pandas.Series(zeros(len(database)))
    database['exptime'] = pandas.Series(zeros(len(database)))
    database['teff'] = pandas.Series(zeros(len(database)))
    database['logg'] = pandas.Series(zeros(len(database)))
    database['feh'] = pandas.Series(zeros(len(database)))
    database['snr'] = pandas.Series(zeros(len(database)))
    database['bcorr'] = pandas.Series(zeros(len(database)))
    database['telluricRV'] = pandas.Series(zeros(len(database)))
    database['vsini'] = pandas.Series(zeros(len(database)))
    database['lsdRV'] = pandas.Series(zeros(len(database)))
    database['ccfRV'] = pandas.Series(zeros(len(database)))
    database['lastupdate'] = pandas.Series(zeros(len(database)))
    for order in range(36):
        database['order'+str(order)] = pandas.Series(zeros(len(database)))

    database.to_csv(basedir+"database.csv",index=False)

def add_to_database(rvobs,database):
    decker = pyfits.getheader(string.replace(rvobs,".rv",""))["DECKER"]
    propid = pyfits.getheader(string.replace(rvobs,".rv",""))["PROPID"]
    rvobs = pandas.read_csv(rvobs,header=None)
    
    header = ["objectname","filepath","bjd","ra","dec","gmag","exptime","teff","logg","feh","snr","bcorr","telluricRV","vsini","lsdRV","ccfRV"]
    headerlen = len(header)
    for i in range(len(header),len(rvobs.iloc[0])):
        header.append("order"+str(i-headerlen))

    rvobs.columns = header
    rvobs['decker'] = decker
    rvobs['propid'] = propid
    rvobs['lastupdate'] = datetime.now()
    sb2prob = is_sb2.checksb2(rvobs['filepath'].iloc[0])
    rvobs['SB2_prob'] = sb2prob
    
    print rvobs

    #database = pandas.read_csv(basedir+"database.csv")

    ### check if observation exists
    if min(abs(database["bjd"]-float(rvobs["bjd"].iloc[0]))) < 10./(60*60*24): ### less than 1s:
        print "Observation exists"
        mask = abs(database["bjd"]-float(rvobs["bjd"].iloc[0])) < 10./(60*60*24)
        lastupdate = database["lastupdate"][mask].iloc[0]
        rvobs['lastupdate'] = lastupdate        
        indx = arange(len(database))
        database = database.drop(indx[mask])

    print "appending observation"
    database = database.append(rvobs.iloc[0],sort=False,ignore_index=True)
    return database

def add_folder(folder):
    rvlist = sort(glob.glob(folder+"*.rv"))
    database = pandas.read_csv(basedir+"database.csv")

    for rvobs in rvlist:
        print rvobs
        database = add_to_database(rvobs,database)

    database.to_csv(basedir+"database.csv",index=False)
        
if __name__ == "__main__":


    if sys.argv[1] == "--create_new":
        create_database()
    
    if sys.argv[1] == "--add_rv":
        rvobs = sys.argv[2]
        add_to_database(rvobs)

    if sys.argv[1] == "--add_folder":
        folder = sys.argv[2]
        add_folder(folder)
