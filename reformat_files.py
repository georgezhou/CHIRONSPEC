import os,sys,glob,string,pyfits
from numpy import *
import matplotlib.pyplot as plt
import setjd


combineadjacent = True
outfolder = "/data/tfop/chiron_data/reduced/"

def get_object_list(fitslist):
    objectlist = []
    tharlist = []
    for fits in fitslist:
        if not os.path.exists(fits+".rv"):
            header = pyfits.getheader(fits)
            jd = float(header["JD"])
            objectname = header["OBJECT"]
            decker = header["DECKER"]
            if objectname in ["ThAr"]:
                tharlist.append([fits,jd,decker])
            else:
                objectlist.append([fits,jd,objectname,decker])

    return objectlist,tharlist

def combine_order(order_specs):
    medvals = nanmedian(order_specs,axis=1)
    normed_spec = order_specs.copy()/transpose([medvals])
    init_median = median(normed_spec,axis=0)

    
    errors = sqrt(order_specs)/order_specs
    for i in range(len(normed_spec[0])):

        ### if val deviates 2x from photon error, reject and replace with median
        mask = abs(normed_spec[:,i]-init_median[i]) > 3*errors[:,i]        
        order_specs[:,i][mask] = init_median[i]*medvals[mask]
        
    return sum(order_specs,axis=0)


def refit_wavelength(data):
    x = arange(len(data[:,0]))
    x.astype(float)
    w = polyfit(x,data[:,0],10)
    w = polyval(w,x)

    data[:,0] = w
    return data

def write_new_fits(headerlist,data):
    if len(headerlist) == 1:
        outheader = headerlist[0]
    if len(headerlist) > 1:
        outheader = headerlist[0]

        jdlist = []
        bjdlist = []
        hjdlist = []
        bcvlist = []
        exptimelist = []
        airmasslist = []

        for i in range(len(headerlist)):
            jdlist.append(float(headerlist[i]["JD"]))
            bjdlist.append(float(headerlist[i]["BJD"]))
            hjdlist.append(float(headerlist[i]["HJD"]))
            bcvlist.append(float(headerlist[i]["BCV"]))
            exptimelist.append(float(headerlist[i]["EXPTIME"]))
            airmasslist.append(float(headerlist[i]["AIRMASS"]))

        outheader.set("JD",mean(jdlist))
        outheader.set("BJD",mean(bjdlist))
        outheader.set("HJD",mean(hjdlist))
        outheader.set("BCV",mean(bcvlist))
        outheader.set("EXPTIME",sum(exptimelist))
        outheader.set("AIRMASS",sum(airmasslist))


    decker = outheader["decker"]
    flatfield = "flats/chi190211."+decker+"flat.fits"
    flatfield = pyfits.open(flatfield)[0].data

    shape_data = array(shape(data))
    shape_data[2] = shape_data[2] + 1
    newdata = zeros(shape_data)
    for order in range(len(data)):
        flat = flatfield[2,len(data)-order-1,:]
        flat /= max(flat)
        flat = flat[::-1]

        
        newdata[order,:,0] = data[order,:,0]
        newdata[order,:,1] = data[order,:,1]
        newdata[order,:,2] = data[order,:,1]/flat

        # print flat

        # plt.subplot(211)
        # plt.plot(newdata[order,:,0],newdata[order,:,1])
        # plt.plot(newdata[order,:,0],newdata[order,:,2])
        # plt.subplot(212)

        # plt.plot(newdata[order,:,0],flat)
        # plt.show()
    
    #sys.exit()
    data = newdata
        
    output_name = "CHIRON_"+string.replace(outheader["DATE-OBS"],":","-")+"_"+outheader["OBJECT"]+".fits"
    
    hduspec = pyfits.PrimaryHDU(data,header=outheader)
    hdulist = pyfits.HDUList([hduspec])
    hdulist.writeto(outfolder+output_name,clobber=True)

                               
    
def main(folder):
    os.system("mkdir "+folder+"rmed/")
    os.system("mv "+folder+"chi*.fits "+folder+"rmed/")
    os.system("mkdir "+outfolder)
    setjd.setjd_folder(folder)
    
    fitslist = sort(glob.glob(folder+"achi*.fits"))    
    objectlist,tharlist = get_object_list(fitslist)

    already_reduced = []

    i = 0
    while i < len(objectlist):

        if not i in already_reduced:

            ### Create a list of headers and data arrays for each spectrum
            ### also check here for multiple obs of same target

            headerlist = []
            wavelist = []
            datalist = []

            already_reduced.append(i)

            
            headerlist.append(pyfits.getheader(objectlist[i][0]))
            datalist.append(pyfits.getdata(objectlist[i][0]))
                              
            j = i+1
            while j < len(objectlist):
                if objectlist[j][2] == objectlist[i][2] and objectlist[i][3] == objectlist[j][3]:
                    if abs(objectlist[j][1] - objectlist[i][1]) < 0.2 and combineadjacent:
                        headerlist.append(pyfits.getheader(objectlist[j][0]))
                        datalist.append(pyfits.getdata(objectlist[j][0]))
                        already_reduced.append(j)

                j += 1



            thar_jd_diff = 0.5
            goodthar = None
            for j in range(len(tharlist)):
                if objectlist[i][1] - tharlist[j][1] > 0 and objectlist[i][1] - tharlist[j][1] < thar_jd_diff and objectlist[i][3] == tharlist[j][2]:
                    thar_jd_diff = objectlist[i][1] - tharlist[j][1]
                    goodthar = tharlist[j][0]

            if not goodthar == None:
                wavelist.append(pyfits.getdata(goodthar))
                        
            thar_jd_diff = 0.5
            goodthar = None
            for j in range(len(tharlist)):
                if  tharlist[j][1]-objectlist[i][1] > 0 and  tharlist[j][1] - objectlist[i][1] < thar_jd_diff and objectlist[i][3] == tharlist[j][2]:
                    thar_jd_diff = tharlist[j][1]-objectlist[i][1] 
                    goodthar = tharlist[j][0]
            if not goodthar == None:
                wavelist.append(pyfits.getdata(goodthar))
                        

            ### Combine the data arrays

            if len(datalist) == 1:
                master_data = datalist[0]

            else:
                master_data = datalist[0].copy()
                for order in range(len(master_data)):
                    order_array = []
                    for j in range(len(datalist)):
                        order_array.append(datalist[j][order][:,1])
                    order_array = array(order_array)
                    
                    master_data[order][:,1] = combine_order(order_array)

            ### combine the wave arrays
            if len(wavelist) <= 1:
                for order in range(len(master_data)):
                    master_data[order] = refit_wavelength(master_data[order])
                
                pass
            elif len(wavelist) == 2:
                for order in range(len(master_data)):
                    order_array = []
                    for j in range(len(wavelist)):
                        order_array.append(refit_wavelength(wavelist[j][order])[:,0])
                    order_array = array(order_array)
                    order_array = mean(order_array,axis=0)
                    master_data[order][:,0] = order_array

            else:
                print "Error: incorrect number of ThAr for exposure",objectlist[i][0]
                sys.exit(1)
            
            print master_data
            write_new_fits(headerlist,master_data)
            
        i += 1
    


def reformat_all():
    main_folder = "/data/tfop/chiron_data/"
    cwd = os.getcwd()+"/"
    print cwd
    os.chdir(main_folder)
    tgzlist = glob.glob("*.tgz")
    for tgz in tgzlist:
        os.chdir(main_folder)
        #unpacked_list = []
        tgz_unpacked = string.replace(tgz,".tgz","")
        if not os.path.exists(tgz_unpacked):
            print tgz
            print tgz_unpacked
    
            os.system("tar -xvf "+tgz)
            #unpacked_list.append(main_folder+tgz_unpacked+"/")

            os.chdir(cwd)
            main(main_folder+tgz_unpacked+"/")
        
if __name__ == "__main__":
    reformat_all()
