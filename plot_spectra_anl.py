import os,sys,string,pickle
from numpy import *
import matplotlib.pyplot as plt
import matplotlib

matplotlib.pyplot.switch_backend('agg')
import matplotlib.gridspec as gridspec
import pyfits
from scipy import interpolate,optimize
import lsd,cc_rv,spectype
c = 3.*10**5



def normalise(spec,niter=1,sigma_low = 0.05,deg=1):

    ### normalise the spectrum to 1
    
    
    x = arange(len(spec))
    mask = spec == spec
    spec_iter = spec[mask]
    x_iter = x.copy()[mask]

    i = 0
    while i < niter:
        
        fit = polyfit(x_iter,spec_iter,deg)
        fit = polyval(fit,x_iter)

        mask = spec_iter - fit > sigma_low * nanstd(spec_iter-fit)
        spec_iter = spec_iter[mask]
        x_iter = x_iter[mask]
        i += 1

    fit = polyfit(x_iter,spec_iter,deg)
    fit = polyval(fit,x)

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

        mask = spec_iter - fit > sigma_low * nanstd(spec_iter-fit)
        spec_iter = spec_iter[mask]
        x_iter = x_iter[mask]

        fit = polyfit(x_iter,spec_iter,1)
        fit = polyval(fit,x)


        spec /= fit

    #mask = spec > 1.1
    #spec[mask] = 1.
    
    return spec



def normalise_with_template(spectrum,template,vsini=10,vshift=0):

    spectrum[:,1] = normalise(spectrum[:,1],deg=2)


    
    template_mask = template[:,0] > nanmin(spectrum[:,0])-10
    template_mask *= template[:,0] < nanmax(spectrum[:,0])+10
    template = template[template_mask]


    vaxis = template[:,0]
    vaxis = c*(vaxis-median(vaxis))/vaxis

    vprof = spectype.make_rot_prof(vaxis,vsini=vsini,epsilon=0.5,scale=1,vshift=vshift)
    vprof = spectype.vgauss(vaxis,vprof,macro=3.75/2.355,scale=1)
    vprof /= sum(vprof)
    template[:,1] = convolve(template[:,1],vprof,mode="same")


    
    ### Now fit template continuum to spectrum
    fit = interpolate.splrep(template[:,0],template[:,1],k=1)
    fit = interpolate.splev(spectrum[:,0],fit)

    diff = spectrum[:,1]/fit
    mask = diff-nanmedian(diff) < 1
    mask *= diff-nanmedian(diff) > -1

    zpt = nanmedian(spectrum[:,0][mask])
    

    fit = polyfit(spectrum[:,0][mask]-zpt,diff[mask],2)
    fit = polyval(fit,spectrum[:,0]-zpt)
            

    spectrum[:,1] /= fit

    return spectrum

def plotMgb(spectrum_array,template_master,ax,teff,logg,feh,vsini=12,vshift=0,decker="fiber"):
    
    spec = []
    for order in range(len(spectrum_array)):
        if max(spectrum_array[order][:,0]) > 5100 and min(spectrum_array[order][:,0]) < 5300:
            print order
            spec_normed = spectrum_array[order].copy()
            if decker == "slicer":
                spec_normed = spec_normed[300:-300]
            else:
                spec_normed = spec_normed
                
            #spec_normed[:,1] = normalise(spec_normed[:,1],niter=5,sigma_low = 0.05,deg=5)

            try:
                spec_normed = normalise_with_template(spec_normed,template_master.copy(),vsini=vsini,vshift=vshift)
            except TypeError:
                spec_normed[:,1] = normalise(spec_normed[:,1],niter=5,sigma_low = 0.05,deg=5)
                
            plt.plot(spec_normed[:,0],spec_normed[:,1],color="b",lw=2)

    template_mask = template_master[:,0] > 5000
    template_mask *= template_master[:,0] < 5400


    vaxis = template_master[:,0][template_mask]
    vaxis = c*(vaxis-median(vaxis))/vaxis

    if vsini < 1:
        vsini = 1

    vprof = spectype.make_rot_prof(vaxis,vsini=vsini,epsilon=0.5,scale=1,vshift=vshift)
    vprof = spectype.vgauss(vaxis,vprof,macro=3.75/2.355,scale=1)
    vprof /= sum(vprof)
    template_f = convolve(template_master[:,1][template_mask],vprof,mode="same")


    plt.plot(template_master[:,0][template_mask],template_f,color="r")

    plt.xlim(5150,5250)
    plt.ylim(0.02,1.2)

    [i.set_linewidth(2.) for i in ax.spines.itervalues()]
    plt.xlabel("Wavelength $(\AA)$",fontsize=15,weight="black")
    plt.text(0.99,0.05,"$v\sin i = "+str(round(vsini,1))+"$ km/s",transform=ax.transAxes,ha="right",fontsize=14)
    plt.text(0.99,0.2,"$T_\mathrm{eff} = "+str(int(teff))+"$ K",transform=ax.transAxes,ha="right",fontsize=14)
    plt.text(0.99,0.15,"$\log g = "+str(round(logg,2))+"$",transform=ax.transAxes,ha="right",fontsize=14)
    plt.text(0.99,0.1,"[Fe/H]$ = "+str(round(feh,2))+"$",transform=ax.transAxes,ha="right",fontsize=14)

    #if vsini > 12.5:
    #    vsini_val = sqrt(vsini**2-12.5**2)
    #    plt.text(0.99,0.05,"$v\sin i = "+str(round(vsini,1))+"$ km/s",transform=ax.transAxes,ha="right",fontsize=14)
    #else:
    #plt.text(0.99,0.05,"$v\sin i$ not measurable",transform=ax.transAxes,ha="right",fontsize=14)
        

def plot_ccf(spectrum_array,template,ax,vsini=10,vshift=0,bcv=0):

    vel_master = arange(-300,300,0.1)
    cc_master = []
    
    for order in range(len(spectrum_array)):
        if median(spectrum_array[order][:,0]) > 5000 and median(spectrum_array[order][:,0]) < 5500:

            try:
                spec_normed = normalise_with_template(spectrum_array[order],template.copy(),vsini=vsini,vshift=vshift)
                epos,drv,cc = cc_rv.cross_correlate_order(spec_normed[:,0],spec_normed[:,1],template)
                cc /= max(cc)

                cc = interpolate.splrep(drv,cc)
                cc = interpolate.splev(vel_master,cc)
                cc_master.append(cc)

            except TypeError:
                pass


    cc = mean(array(cc_master),axis=0)
    plt.plot(vel_master+bcv,cc,color="k",lw=2)

            
    [i.set_linewidth(2.) for i in ax.spines.itervalues()]

    vmin = vshift+bcv-10*vsini 
    vmax = vshift+bcv+10*vsini

    if vmin < min(vel_master):
        vmin = min(vel_master)
    if vmax > max(vel_master):
        vmax = max(vel_master)

    if vsini < 1:
        vmin = -300
        vmax = 300
    
    plt.xlim(vmin,vmax)
    plt.ylim(-0.3,1.1)

    plt.xlabel("Shift (km/s)",fontsize=15,weight="black")




def plot_lsd(lsd_master,ax,vsini=10,vshift=0,bcv=0):

    plt.plot(lsd_master[:,0]+bcv,lsd_master[:,1]/max(lsd_master[:,1]),color="k",lw=2)
            
    [i.set_linewidth(2.) for i in ax.spines.itervalues()]

    vmin = vshift+bcv-10*vsini 
    vmax = vshift+bcv+10*vsini

    if vmin < min(lsd_master[:,0]):
        vmin = min(lsd_master[:,0])
    if vmax > max(lsd_master[:,0]):
        vmax = max(lsd_master[:,0])
    if vsini < 1:
        vmin = -300
        vmax = 300
    
    plt.xlim(vmin,vmax)
    
    plt.ylim(-0.3,1.1)

    plt.xlabel("Shift (km/s)",fontsize=15,weight="black")


def plot_halpha(spectrum_array,ax):
    for order in range(len(spectrum_array)):
        if min(spectrum_array[order][:,0]) < 6562.8 and max(spectrum_array[order][:,0]) > 6562.8:
            
            plt.plot(spectrum_array[order][:,0],spectrum_array[order][:,1],color="k",lw=1)

            plt.xlim(min(spectrum_array[order][:,0]),max(spectrum_array[order][:,0]))
            break

    plt.axvline(x=6562.8,ymax=1,ymin=0.8,color="b",lw=2)
    
    plt.xlabel("Wavelength $(\AA)$",fontsize=15,weight="black")
    plt.ylabel("ADU",fontsize=15,weight="black")

def plot_NaD(spectrum_array,ax):
    for order in range(len(spectrum_array)):
        if min(spectrum_array[order][:,0]) < 5892 and max(spectrum_array[order][:,0]) > 5892:
            
            plt.plot(spectrum_array[order][:,0],spectrum_array[order][:,1],color="k",lw=1)

            plt.xlim(min(spectrum_array[order][:,0]),max(spectrum_array[order][:,0]))
            break

    plt.axvline(x=5889.950,ymax=1,ymin=0.8,color="b",lw=2)
    plt.axvline(x=5895.924,ymax=1,ymin=0.8,color="b",lw=2)
     
    plt.xlabel("Wavelength $(\AA)$",fontsize=15,weight="black")
    plt.ylabel("ADU",fontsize=15,weight="black")

    
def plot_Li(spectrum_array,ax):
    for order in range(len(spectrum_array)):
        if min(spectrum_array[order][:,0]) < 6708 and max(spectrum_array[order][:,0]) > 6708:
            
            plt.plot(spectrum_array[order][:,0],spectrum_array[order][:,1],color="k",lw=1)

            xmin = min(spectrum_array[order][:,0])
            xmax = max(spectrum_array[order][:,0])

            if xmin < 6708-30:
                xmin = 6708-30
            if xmax > 6708+30:
                xmax = 6708+30
            
            plt.xlim(xmin,xmax)
            break

    plt.axvline(x=6707.7,ymax=1,ymin=0.8,color="b",lw=2)
     
    plt.xlabel("Wavelength $(\AA)$",fontsize=15,weight="black")
    plt.ylabel("ADU",fontsize=15,weight="black")


def main(spectrum_file,template_file,teff,logg,feh,vsini,shift):
    ### format the spectrum file
    spectrum_hdulist = pyfits.open(spectrum_file)
    ra,dec,objectname,decker,bcv = spectrum_hdulist[0].header["RA"], spectrum_hdulist[0].header["DEC"], spectrum_hdulist[0].header["OBJECT"],spectrum_hdulist[0].header["DECKER"],float(spectrum_hdulist[0].header["BCV"])
    template_array = loadtxt(template_file)
    
    spectrum_name = os.path.basename(spectrum_file)
    spectrum_path = string.replace(spectrum_file,spectrum_name,"")

    lsd_list,lsd_master,dummy,dummy = pickle.load(open(spectrum_path+"lsd_"+spectrum_name+".pkl","rb"))

    spectrum_array = []
    for order in range(0,len(spectrum_hdulist[0].data)):
        wave = spectrum_hdulist[0].data[order,:,0][20:-20]
        flux = spectrum_hdulist[0].data[order,:,2][20:-20]
        mask = flux != flux
        flux[mask] = 1.
            
        spectrum_array.append(transpose(array([wave,flux])))


    ### CREATE GRID FOR SPECTRUM
    plt.figure(figsize=(25,15))
    plt.subplots_adjust(left=0.05,right=0.95,bottom=0.1,top=0.9)
    gs = gridspec.GridSpec(2, 4)


    ### Plot the mgb order    
    ax = plt.subplot(gs[0, 0:3])
    plt.title(objectname+" "+ra+" "+dec)
    plotMgb(spectrum_array,template_array,ax,teff,logg,feh,vsini=vsini,vshift=shift,decker=decker)

    ### Plot the ccf for a couple of orders
    ax = plt.subplot(gs[1, 0])
    try:
        plot_ccf(spectrum_array,template_array,ax,vsini=vsini,vshift=shift,bcv=bcv)
    except ValueError:
        pass
    plt.title("Non-rotating template correlation")


    ### Plot the lsd profile
    ax = plt.subplot(gs[1, 1])
    plot_lsd(lsd_master,ax,vsini=vsini,vshift=shift,bcv=bcv) 
    plt.title("LSD line profile")


    ### Plot the Halpha profile
    ax = plt.subplot(gs[0, 3])
    plot_halpha(spectrum_array,ax)
    plt.title("H-alpha")

    ### Plot the NaD profile
    ax = plt.subplot(gs[1, 2])
    plot_NaD(spectrum_array,ax)
    plt.title("Na-D")
    
    ### Plot the Li profile
    ax = plt.subplot(gs[1, 3])
    plot_Li(spectrum_array,ax)
    plt.title("Li")
    
    plt.savefig(spectrum_path+spectrum_name+".pdf")

if __name__ == "__main__":
    
    spectrum_file = "/data/tfop/chiron_data/reduced/CHIRON_2020-01-13T07-18-00.3_T0400595342.fits"
    template_file = "spectral_library/template_5750_4.5_0.0.dat"

    teff = 5727
    logg = 2.863235151436735
    feh = -0.2382474843289443
    vsini = 0.17080534296369684
    shift = 140.53363940746706

    main(spectrum_file,template_file,teff,logg,feh,vsini,shift)
    plt.show()
