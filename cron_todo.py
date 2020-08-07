import os,sys,glob,string,pandas
from datetime import datetime,timedelta
from numpy import *

os.chdir("/data/tfop/chiron_analysis/")


run_reduction = True

now = datetime.now()
print now


### check when the last reduction was
try:
    list_of_files = glob.glob('/data/tfop/chiron_data/reduced/*')
    latest_file = max(list_of_files, key=os.path.getctime)
    last_reduction =  datetime.fromtimestamp(os.path.getctime(latest_file))
    print "last reduction was",last_reduction
    
    if now - last_reduction < timedelta(minutes=30):
        run_reduction = False ### reduction script is probably already running
        print "a script is probably already running"

except ValueError:
    pass



### check if there are any new files waiting to be reduced
if run_reduction:

    import wget_newdata
    maillist = pandas.read_csv("maillist.csv")
    for user in unique(maillist["user"]):
        print "checking for new observations from ",user
        wget_newdata.main(user)

    
    list_of_files = glob.glob('/data/tfop/chiron_data/*')
    
    import reformat_files
    print "reformatting"
    reformat_files.reformat_all()

    # if len(glob.glob('/data/tfop/chiron_data/*')) == len(list_of_files):
    #     print "no new files generated, nothing to reduce"
    #     run_reduction = False

        
if run_reduction:

    if len(glob.glob('/data/tfop/chiron_data/reduced/*.fits')) == len(glob.glob('/data/tfop/chiron_data/reduced/*.rv')):
        print "no new files generated, nothing to reduce"
        run_reduction = False

    
        

if run_reduction:
    print "running reduction pipeline",now
    os.chdir("/data/tfop/chiron_analysis/")


    print "running analysis"
    import run_analysis
    run_analysis.main("/data/tfop/chiron_data/reduced/")

    import database
    database.add_folder("/data/tfop/chiron_data/reduced/")

    import makereport
    makereport.makereport_all(now)

