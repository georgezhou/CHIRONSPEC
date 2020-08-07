import os,sys,pickle,pyfits
import pandas
from numpy import *
import matplotlib.pyplot as plt
from scipy import interpolate,signal
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor

c = 299792.458
chiron_resolution = 80000
database = pandas.read_csv("database.csv")
speclib = "spectral_library/"


def normalise(spec,niter=1,sigma_low = -0.05,deg=1):

    try:

        ### normalise the spectrum to 1


        x = arange(len(spec))
        mask = spec == spec
        spec_iter = spec[mask]
        x_iter = x.copy()

        i = 0
        while i < niter:
            fit = polyfit(x_iter,spec_iter,deg)
            fit = polyval(fit,x_iter)

            mask = spec_iter - fit > sigma_low * std(spec_iter-fit)
            mask *= spec_iter - fit < 2. * std(spec_iter-fit)
            spec_iter = spec_iter[mask]
            x_iter = x_iter[mask]
            i += 1

        fit = polyfit(x_iter,spec_iter,deg)
        fit = polyval(fit,x)
        # plt.plot(spec)
        # plt.plot(fit)
        # plt.show()

        if min(fit) > 0:
            spec /= fit

        else:

            max_prefit = max(fit)
            spec -= fit
            spec += max_prefit

            mask = spec == spec
            spec_iter = spec[mask]
            x_iter = x[mask]

            fit = polyfit(x_iter,spec_iter,1)
            fit = polyval(fit,x_iter)

            mask = spec_iter - fit > sigma_low * std(spec_iter-fit)
            spec_iter = spec_iter[mask]
            x_iter = x_iter[mask]

            fit = polyfit(x_iter,spec_iter,1)
            fit = polyval(fit,x)


            spec /= fit

        #mask = spec > 1.1
        #spec[mask] = 1.


    except TypeError:
        # x = arange(len(spec))
        # mask = spec == spec
        # spec_iter = spec[mask]
        # x_iter = x.copy()[mask]
        # f = polyfit(x_iter,spec_iter,1)
        # f = polyval(f,x)
        # spec = spec / f

        pass
        
    return spec



def write_out_spec(spec):

    spectrum = pyfits.open(spec)
    spectrum = spectrum[0].data

    odata = []
    for order in arange(len(spectrum)):
        wave = spectrum[order,:,0]
        flux = spectrum[order,:,2]
        weights = spectrum[order,:,1]
        weights = signal.medfilt(weights, kernel_size=101)

        ### linearise
        indx = argsort(wave)
        wave_out = linspace(min(wave),max(wave),len(wave))
        vel_out = (wave_out-median(wave_out))/wave_out * c
        
        flux_out = interpolate.splrep(wave[indx],flux[indx],k=1)
        flux_out = interpolate.splev(wave_out,flux_out)
        weights_out = interpolate.splrep(wave[indx],weights[indx],k=1)
        weights_out = interpolate.splev(wave_out,weights_out)
                                        
        if chiron_resolution > 50000:
            gaussian_width = sqrt((c/44000.)**2 - (c/chiron_resolution)**2)
            gaussian_width /= 2.355

            gaussian = exp(-((vel_out)**2/(2*gaussian_width**2)))
            gaussian /= sum(gaussian)

            flux_out = convolve(flux_out, gaussian, mode="same")


        out = transpose(array([order*ones(len(wave_out)),wave_out,flux_out,weights_out]))

        ### find overlap
        if len(odata) > 1:
            mask_odata = array(odata)[:,1] > min(out[:,1])
            mask_odata *= array(odata)[:,1] < max(out[:,1])
            mask_out = out[:,1] > min(array(odata)[:,1])
            mask_out *= out[:,1] < max(array(odata)[:,1])

            if nanmedian(array(odata)[mask_odata][:,2]) > nanmedian(out[mask_out][:,2]):
                out = out[invert(mask_out)]
            else:
                odata = list(array(odata)[invert(mask_odata)])

        
        odata += list(out)
            

    odata = array(odata)
    odata[:,2] /= nanmedian(odata[:,2])
    odata[:,3] /= nanmax(odata[:,3])

    return odata

def binspec(spec):

    
    wave = arange(5000,5300,0.05)
    spec_interp = interpolate.splrep(spec[:,0],spec[:,1],k=1)
    #spec = interpolate.splev(wave,spec)

    binwidth=0.05
    specout = ones(len(wave))
    for i in range(len(wave)):
        mask = abs(spec[:,0]-wave[i])<binwidth/2
        if len(spec[:,1][mask]) > 0:
            specout[i] = nanmedian(spec[:,1][mask])
        else:
            specout[i] = interpolate.splev([wave[i]],spec_interp)[0]

    return specout

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
        

def roundnearest(val,r):
    val = val / r
    val = round(val)
    val = val*r
    return val

def make_spectrum(spectrum):



    ### Find the spectrum in the database:
    for indx in range(len(database)):
        if os.path.basename(database["filepath"].iloc[indx]) == os.path.basename(spectrum):
            break

    #rv = database["lsdRV"].iloc[indx]-database["bcorr"].iloc[indx]
    #vsini_init = database["vsini"].iloc[indx]
    snr = database["snr"].iloc[indx]
    teff_init = int(roundnearest(database["teff"].iloc[indx],250))
    template = speclib+"template_"+str(teff_init)+"_4.5_0.0.dat"

    template = loadtxt(template)

    
    ccf_list,ccf,vsini_init,rv = pickle.load(open(os.path.dirname(spectrum)+"/lsd_"+os.path.basename(spectrum)+".pkl","rb"))
    ccf[:,0] -= rv    
    ccf = binccf(ccf)
    # gaussian_width = c/44000. #sqrt((c/44000.)**2 - (c/chiron_resolution)**2)
    # gaussian_width /= 2.355

    # gaussian = exp(-((arange(-300,300,1.))**2/(2*gaussian_width**2)))
    # gaussian /= sum(gaussian)
    # ccf = convolve(ccf, gaussian, mode="same")
    # ccf /= max(ccf)

    spectrum = write_out_spec(spectrum)

    spec = []
    for order in unique(spectrum[:,0]):
        mask = spectrum[:,0] == order
        spec_i = spectrum[mask][:,1:]
        if nanmedian(spec_i[:,0]) > 4800 and nanmedian(spec_i[:,0]) < 5500:
            spec_i[:,1] = normalise(spec_i[:,1],deg=3,niter=3)
            spec_i[:,0] = spec_i[:,0]-(rv/c)*spec_i[:,0]
            spec += list(spec_i)

    spec = array(spec)
    indx = argsort(spec[:,0])
    spec = spec[indx]
    spec = binspec(spec)
    spec /= median(spec)

    output = list(spec)+list(ccf)
    output = array(output)
    mask = output == output
    output[invert(mask)] = 1.

    # plt.plot(output)
    # plt.show()

    return output,snr

def runML(spectrum):
    library = load("SPEC_ML/library.spc_corrected.npy")
    catalogue = pandas.read_csv("SPEC_ML/library.spc_corrected.csv")
    
    
    spec,snr = make_spectrum(spectrum)

    
    mor = pickle.load(open("SPEC_ML/mor.pkl","rb"))
    #gbr_teff,gbr_logg,gbr_mh,gbr_vsini = mor

    # teff = gbr_teff.predict([spec])
    # logg = gbr_logg.predict([spec])
    # feh = gbr_mh.predict([spec])
    # vsini = gbr_vsini.predict([spec])

    results = mor.predict([spec[:-600]])
    teff,logg,feh,vsini = results[0]
    teff = 10**teff
    

    #sb2 = pickle.load(open("/data/yjzhou0/SPEC_ML/sb2.pkl","rb"))
    #proba = sb2.predict_proba([spec[-600:]])
    proba = 0

    print teff,logg,feh,vsini,proba


    wave = arange(5000,5300,0.05)
    vel = arange(-300,300,1.)
    lsd = spec[-600:]
    spec = spec[:-600]

    plt.figure(figsize=(15,15))
    plt.subplots_adjust(left=0.05,right=0.98,top=0.9,bottom=0.07,hspace=0.3)
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(4,3)

    ax = plt.subplot(gs[0,:])
    plt.title(os.path.basename(spectrum))
    plt.plot(wave,spec,lw=3,alpha=0.5,color="k")
    plt.xlim(5000,5100)
    plt.ylabel("Flux",fontsize=15,weight="black")

    ax = plt.subplot(gs[1,:])
    plt.plot(wave,spec,lw=3,alpha=0.5,color="k")
    plt.xlim(5100,5200)
    plt.ylabel("Flux",fontsize=15,weight="black")

    ax = plt.subplot(gs[2,:])
    plt.plot(wave,spec,lw=3,alpha=0.5,color="k")
    plt.xlim(5200,5300)

    plt.xlabel("Wavelength $(\AA)$",fontsize=15,weight="black")
    plt.ylabel("Flux",fontsize=15,weight="black")

    
    ax = plt.subplot(gs[3,0])
    plt.plot(vel,lsd,lw=3,color="k",alpha=0.5)
    plt.xlim(min(vel),max(vel))
    plt.xlabel("Velocity (km/s)",fontsize=15,weight="black")
    plt.ylabel("LSD",fontsize=15,weight="black")

    # if proba[0][1] > 0.1 and proba[0][1] < 0.5:        
    #     plt.text(0.1,0.9,"Possible blend",color="r",fontsize=15, weight="black",ha="left",transform=ax.transAxes)
    # if proba[0][1] > 0.5:        
    #     plt.text(0.1,0.9,"Likely blend",color="r",fontsize=15, weight="black",ha="left",transform=ax.transAxes)

    
    ax = plt.subplot(gs[3,1])
    plt.scatter(catalogue["Teff_corrected"],catalogue["logg_corrected"],s=5,color="k",alpha=0.3)
    plt.scatter(teff,logg,s=100,color="r",marker="*",edgecolor="r")

    plt.xlim(8000,3500)
    plt.ylim(5,1.5)

    plt.xlabel("$T_\mathrm{eff}$",fontsize=18,weight="black")
    plt.ylabel("$\log g$",fontsize=18,weight="black")


    
    ax = plt.subplot(gs[3,2])
    plt.axis('off')

    plt.text(0.1,0.9,"$T_\mathrm{eff} = "+str(int(teff))+"$ K",color="k",fontsize=15,ha="left",transform=ax.transAxes)
    plt.text(0.1,0.7,"$\log g = "+str(round(logg,2))+"$",color="k",fontsize=15,ha="left",transform=ax.transAxes)
    plt.text(0.1,0.5,"$\mathrm{[M/H]} = "+str(round(feh,2))+"$",color="k",fontsize=15,ha="left",transform=ax.transAxes)
    plt.text(0.1,0.3,"$v\sin I_\star = "+str(round(vsini,2))+"$ km/s",color="k",fontsize=15,ha="left",transform=ax.transAxes)
    plt.text(0.1,0.1,"$S/N = "+str(int(snr))+"$",color="k",fontsize=15,ha="left",transform=ax.transAxes)


    

    
    
    if False:

        importance = gbr_teff.feature_importances_
        importance += gbr_logg.feature_importances_
        importance += gbr_mh.feature_importances_
        importance += gbr_vsini.feature_importances_

        mask = importance != importance
        mask += importance <= 0.
        importance[mask] = 10**-5 #nanmin(importance[invert(mask)])
                                            
        importance = log10(importance)

        print min(importance),max(importance)

        importance_spec = interpolate.splrep(wave,importance[:-600],k=3)
        wave_fine = arange(min(wave),max(wave),0.01)
        importance_spec = interpolate.splev(wave_fine,importance_spec)
        spec_fine = interpolate.splrep(wave,spec)
        spec_fine = interpolate.splev(wave_fine,spec_fine)
        
        importance_lsd = interpolate.splrep(vel,importance[-600:],k=3)
        vel_fine = arange(min(vel),max(vel),0.1)
        importance_lsd = interpolate.splev(vel_fine,importance_lsd)
        lsd_fine = interpolate.splrep(vel,lsd)
        lsd_fine = interpolate.splev(vel_fine,lsd_fine)
                                             

        ax = plt.subplot(gs[0,:])
        mask = wave_fine > 5000
        mask *= wave_fine < 5100
        plt.scatter(wave_fine[mask],spec_fine[mask],c=importance_spec[mask],cmap="hot",edgecolor="none",s=60,alpha=0.2,zorder=2)
        
        ax = plt.subplot(gs[1,:])
        mask = wave_fine > 5100
        mask *= wave_fine < 5200
        plt.scatter(wave_fine[mask],spec_fine[mask],c=importance_spec[mask],cmap="hot",edgecolor="none",s=60,alpha=0.2,zorder=2)
        
        ax = plt.subplot(gs[2,:])
        mask = wave_fine > 5200
        mask *= wave_fine < 5300
        plt.scatter(wave_fine[mask],spec_fine[mask],c=importance_spec[mask],cmap="hot",edgecolor="none",s=60,alpha=0.2,zorder=2)
        
        ax = plt.subplot(gs[3,0])
        plt.scatter(vel_fine,lsd_fine,c=importance_lsd,cmap="hot",edgecolor="none",s=60,alpha=0.2,zorder=2)

            


    if True:
        metric = ((array(catalogue["Teff_corrected"])-teff)/500)**2
        metric += ((array(catalogue["logg_corrected"])-logg)/0.5)**2
        metric += ((array(catalogue["feh_corrected"])-feh)/0.5)**2
        metric += (array(catalogue["vsini_corrected"])-vsini)**2


        indx = argmin(metric)
        ax = plt.subplot(gs[3,1])
        plt.scatter(catalogue["Teff_corrected"].iloc[indx],catalogue["logg_corrected"].iloc[indx],color="b",facecolor="none",s=100,edgecolor="orange",lw=3)
        print catalogue["Teff_corrected"].iloc[indx],catalogue["logg_corrected"].iloc[indx]



        
        model = library[indx]


        ax = plt.subplot(gs[0,:])
        plt.plot(wave,model[:-600]/median(model[:-600]),lw=1.5,alpha=0.8,color="r",zorder=3)

        ax = plt.subplot(gs[1,:])
        plt.plot(wave,model[:-600]/median(model[:-600]),lw=1.5,alpha=0.8,color="r",zorder=3)

        ax = plt.subplot(gs[2,:])
        plt.plot(wave,model[:-600]/median(model[:-600]),lw=1.5,alpha=0.8,color="r",zorder=3)
        
    plt.savefig(spectrum+".png")
    #plt.show()

    plt.clf()
    plt.close()


    return teff,logg,feh,vsini,proba

def runall(folder):
    import glob,pandas
    database = pandas.read_csv("database.csv")[:100]
    for i in range(len(database)):
        fits = database["filepath"].iloc[i]
        
        if os.path.exists(fits):
            print fits
            teff,logg,feh,vsini,proba = runML(fits)
            database[i,"teff"] = teff
            database[i,"logg"] = logg
            database[i,"feh"] = feh

    plt.scatter(database["teff"],database["logg"])
    plt.xlim(9000,3000)
    plt.ylim(5,3)
    plt.show()
        


if __name__ == "__main__":
    #runML(sys.argv[1])
    runall(sys.argv[1])
