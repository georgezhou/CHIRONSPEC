import os,sys,string,pandas
from numpy import *
import matplotlib.pyplot as plt
import astropy.time
from datetime import datetime,timedelta


def mmd(x):
    stdev = nanstd(x)
    mask = abs(x-nanmedian(x)) < 5*stdev
    x = x[mask]
    return 1.35*nanmedian(abs(x-nanmedian(x))) / sqrt(float(len(x)))



database = "database.csv"
database = pandas.read_csv(database)
maillist = pandas.read_csv("maillist.csv")


def makereport(checktime,propid):

    lastupdate = pandas.to_datetime(database["lastupdate"],errors="coerce")
    mask = lastupdate-checktime > timedelta(seconds=0.0)
    
    
    if not propid == None:
        mask *= database["propid"] == propid

    stars = unique(database["objectname"][mask])

    reportout = ""


    for starname in stars:
        reportout += "#######################################################\n"
        reportout += starname+"\n"
        star = database[database["objectname"] == starname]

        for i in range(len(star)):
            filepath = str(star["filepath"].iloc[i])
            filepath = os.path.basename(filepath)
            filepath = string.split(filepath,"_")[1]



            ordervels = []
            for order in arange(1,35):
                ordervels.append(star["order"+str(order)].iloc[i])

            ordervels = array(ordervels)

            errorbar = mmd(ordervels)


            reportout += filepath+" "
            reportout += str(star["bjd"].iloc[i])+" "
            reportout += "RV  "+str(round(star["lsdRV"].iloc[i],3))+" "+str(round(errorbar,3))+" "
            reportout += "Teff "+str(int(star["teff"].iloc[i]))+" "
            reportout += "logg "+str(round(star["logg"].iloc[i],2))+" "
            reportout += "feh "+str(round(star["feh"].iloc[i],2))+" "
            reportout += "vsini "+str(round(star["vsini"].iloc[i],1))+" "
            reportout += "SB2 Probability "+str(round(star["SB2_prob"].iloc[i],1))+"\n"


    timenow = string.replace(str(datetime.now())," ","T")
    out = open("/data/tfop/chiron_data/reports/"+timenow+"_PROPID_"+str(propid),"w")
    out.write(reportout)
    out.close()

    return "/data/tfop/chiron_data/reports/"+timenow+"_PROPID_"+str(propid)
    


def makereport_all(checktime):

    print "making reports for everything updated since",checktime
    lastupdate = pandas.to_datetime(database["lastupdate"],errors="coerce")
    mask = lastupdate-checktime > timedelta(seconds=0.0)
    
    unique_propids =  unique(database["propid"][mask])

    for propid in unique_propids:
        print "making report for",propid
        report = makereport(checktime,propid)

        ### find who to mail
        mask = maillist["propid"] == propid

        #for emails in maillist["email"][mask]:
        for i in range(len(maillist[mask])):
            emails = maillist["email"][mask].iloc[i]
            user = maillist["user"][mask].iloc[i]
            mail_cmd = "cat "+report+" | mail -s \"New Chiron data for "+str(propid)+" "+user+" dated " + str(checktime) + "\" "+emails       
            print mail_cmd
            os.system(mail_cmd)
            

if __name__ == "__main__":
    
    checktime = datetime.strptime("2020-05-21 10:00:00","%Y-%m-%d %H:%M:%S")
    print checktime
    makereport_all(checktime)
