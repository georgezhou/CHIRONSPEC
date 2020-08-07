import os,sys,string
from numpy import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.pyplot.switch_backend('agg')
import pyfits
from scipy import interpolate,optimize,signal,stats
from PyAstronomy import pyasl
import spectype
import lsd

c = 3*10**5.

def normalise(spec,niter=1,sigma_low = 0.05,deg=5):

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
        spec_iter = spec_iter[mask]
        x_iter = x_iter[mask]
        i += 1

    fit = polyfit(x_iter,spec_iter,deg)
    fit = polyval(fit,x)
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

    mask = spec > 1.1
    spec[mask] = 1.
    
    return spec


def apodize(x,ap=50):
    func = ones(len(x))
    xpos = arange(len(x))

    mask = xpos<ap
    func[mask] = 0.5*cos((xpos[mask]+ap)*2*pi/(ap*2))+0.5
    
    mask = max(xpos)-xpos < ap
    func[mask] = 0.5*cos((max(xpos)-ap-xpos[mask])*2*pi/(ap*2))+0.5

    return func*x

def fitquad(x,y,dp=50):
    pos = arange(len(x))
    xcenter = pos[argmax(y)]
    xmin = xcenter-dp
    xmax = xcenter+dp
    if xmin < 0:
        xmin = 0
    if xmax > len(pos):
        xmax = len(pos)


    x,y = x[xmin:xmax],y[xmin:xmax]
    x0 = polyfit(x,y,2)
    epos = -x0[1]/(2*x0[0])

    f = polyval(x0,x)
    # plt.plot(x,y)
    # plt.plot(x,f)
    # plt.show()

    return epos
    


def cross_correlate_order(wave,flux,template,ap=200):
    mask = template[:,0] > min(wave)-1000*median(wave)/c
    mask *= template[:,0] < max(wave)+1000*median(wave)/c

    flux = normalise(flux,deg=1)
    template_flux = normalise(template[:,1][mask],deg=0)
    #template_flux = template[:,1][mask]

    ### do some more normalisation against the template
    fit = interpolate.splrep(template[:,0][mask],template_flux)
    fit = interpolate.splev(wave,fit)

    diff = flux/fit
    diffmask = diff-median(diff) < 0.2
    diffmask *= diff-median(diff) > -0.6

    fit = polyfit(wave[diffmask],diff[diffmask],2)
    fit = polyval(fit,wave)
    flux /= fit

    flux = -1*(flux-1)
    template_flux = -1*(template_flux-1)


    #flux = apodize(flux,ap=ap)
    #template_flux = apodize(template_flux,ap=ap)

    drv,cc = pyasl.crosscorrRV(wave, flux, template[:,0][mask], template_flux, -300, +300, 0.01, mode='doppler', skipedge=0, edgeTapering=None)
    # plt.subplot(211)
    # plt.plot(wave,flux)
    # plt.plot(template[:,0][mask],template_flux)

    # plt.subplot(212)
    # plt.plot(drv,cc)
    # plt.show()

    try:
        #epos = pyasl.quadExtreme(drv, cc, mode='max', dp=(50, 50), exInd=None, fullOutput=False, fullPoint=False)[0]
        dp = len(cc[cc > 0.5*nanmax(cc)])/2
        print "len of cc peak find",dp
        epos = fitquad(drv,cc,dp=dp)
        print "rv",epos
    except:
        epos = nan
    

    return epos,drv,cc
def average_ccf(ccf_list,vel,vsini=10,shift=0):

    vgrid = arange(-300,300,0.01)
    good_ccf = []
    weights = []

    for ccf in ccf_list:

        if std(ccf) == std(ccf):

            mask_width = vsini*2
            if mask_width < 3:
                mask_width = 3

            mask = vel-shift > -1*mask_width
            mask *= vel-shift < mask_width

            out_mask1 = ccf[0]-shift > vsini
            out_mask2 = ccf[0]-shift < -1*vsini
            out_mask = out_mask1 + out_mask2        

            ccf /= nanmedian(ccf[mask])

            interp = interpolate.splrep(vel,ccf,k=1)
            interp = interpolate.splev(vgrid,interp)

            if median(interp) == median(interp):
                good_ccf.append(interp)
                weights.append(nanstd(ccf[out_mask]))
            
   
    good_ccf = array(good_ccf)
    weights = 1/array(weights)**2
    weights = weights / nansum(weights)

    if median(weights) != median(weights):
        weights = ones(len(good_ccf))

    ccf_final = zeros(len(good_ccf[0]))
    for i in range(len(good_ccf)):
        ccf_final += good_ccf[i]*weights[i]


    ccf_o = transpose(array([vgrid,ccf_final]))

    return ccf_o


def main(spectrum,template,vsini=10):
    spectrum_hdulist = pyfits.open(spectrum)
    
    template = loadtxt(template)

    cc_list = []

    #wavemin = 4500
    #wavemax = 6200

    wavemin = 5500
    wavemax = 6200

    rvlist = []
    for order in range(0,len(spectrum_hdulist[0].data)):
    #for order in [13]:
        wave = spectrum_hdulist[0].data[order,:,0][20:-20]
        flux = spectrum_hdulist[0].data[order,:,2][20:-20]

        if median(wave) > wavemin and median(wave) < wavemax:

            try:
                epos,drv,cc = cross_correlate_order(wave,flux,template)
                cc_list.append(cc)
            except TypeError:
                epos = nan
                pass



            rvlist.append([order,epos])

    rvlist = array(rvlist)

    cc_list = array(cc_list)
    ccf = average_ccf(cc_list,drv,vsini=10,shift=nanmedian(rvlist))
    #epos = pyasl.quadExtreme(ccf[:,0], ccf[:,1], mode='max', dp=(50, 50), exInd=None, fullOutput=False, fullPoint=False)[0]
    dp = len(cc[cc > 0.5*nanmax(cc)])/2
    print "len of cc peak find",dp
    epos = fitquad(ccf[:,0],ccf[:,1],dp=dp)
    
    #plt.plot(ccf[:,0],ccf[:,1])
    #plt.show()

    return rvlist,epos

def remove_stellar_template(template_array,spectrum_i,lsd):

    
    template_mask = template_array[:,0] > min(spectrum_i[:,0])-10
    template_mask *= template_array[:,0] < max(spectrum_i[:,0])+10
    template_i = template_array[template_mask]

    vel = c*(template_i[:,0]-median(template_i[:,0]))/template_i[:,0]
    lsd_interp = interpolate.splrep(lsd[:,0],lsd[:,1])
    lsd_interp = interpolate.splev(vel,lsd_interp,ext=1)
    lsd_interp /= sum(lsd_interp)

    template_i[:,1] = convolve(template_i[:,1],lsd_interp,mode="same")

    mask = template_i[:,1] == template_i[:,1]
    
    template_i_interp = interpolate.splrep(template_i[:,0][mask],template_i[:,1][mask])
    template_i_interp = interpolate.splev(spectrum_i[:,0],template_i_interp,ext=1)
    mask = template_i_interp != 0
    # mask *= spectrum_i[:,0] == spectrum_i[:,0]
    # mask *= spectrum_i[:,1] == spectrum_i[:,1]
    # mask *= abs(spectrum_i[:,1]) < inf


    offset = spectrum_i[:,1]/template_i_interp
    fit = polyfit(spectrum_i[:,0][mask],offset[mask],2)
    fit = polyval(fit,spectrum_i[:,0])
    template_i_interp *= fit

    

    spectrum_corrected = spectrum_i[:,1] / template_i_interp

    # plt.subplot(211)
    # plt.plot(spectrum_i[:,0],spectrum_i[:,1])
    # plt.plot(spectrum_i[:,0],template_i_interp)

    # plt.subplot(212)
    # plt.plot(spectrum_i[:,0],spectrum_corrected)
    # plt.show()
    
    return spectrum_corrected

def runlsd(spectrum,template):
    spectrum = spectrum[argsort(spectrum[:,0])]

    mask = template[:,0] > min(spectrum[:,0])-1000*median(spectrum[:,0])/c
    mask *= template[:,0] < max(spectrum[:,0])+1000*median(spectrum[:,0])/c
    template = template[mask]
    template[:,1] = normalise(template[:,1],deg=0)
    
    ### Now do lsd 10 times and create an average lsd profile\
    wave0 = min(spectrum[:,0])
    wave1 = max(spectrum[:,0])
    dwave = (wave1-wave0)/float(len(spectrum))

        
    spectrum_interp,template_interp,vdel=lsd.interpolate_spectra(spectrum,template)

    vel,ccf = lsd.lsd(template_interp,spectrum_interp,vdel)
    lsd_prof_list = []

    interp = interpolate.splrep(spectrum[:,0],spectrum[:,1],k=1)

    for i in range(10):
        wave_i = arange(random.normal(wave0,2.),random.normal(wave1,2.),dwave)
        spectrum_i = interpolate.splev(wave_i,interp)
        spectrum_i = transpose(array([wave_i,spectrum_i]))
        spectrum_interp,template_interp,vdel=lsd.interpolate_spectra(spectrum_i,template)

        spectrum_interp = apodize(1-spectrum_interp)
        spectrum_interp = 1-spectrum_interp
        template_interp = apodize(1-template_interp)
        template_interp = 1-template_interp

        vel_i,ccf_i = lsd.lsd(template_interp,spectrum_interp,vdel)

        ccf_i = interpolate.splrep(vel_i,ccf_i,k=1)
        ccf_i = interpolate.splev(vel,ccf_i)

        lsd_prof_list.append(ccf_i)

    ccf = median(array(lsd_prof_list),axis=0)

    ccf = transpose(array([vel,ccf]))

    #plt.plot(ccf[:,0],ccf[:,1])
    #plt.show()
    return ccf
    



def measure_telluric_rv(spectrum,stellar_template,vsini,vshift):
    spectrum_hdulist = pyfits.open(spectrum)
    stellar_template = loadtxt(stellar_template)
    template = loadtxt("transmission.dat")
    template[:,0] *= 10.
    epos = 0

    for order in range(1,len(spectrum_hdulist[0].data)):
        #s = read_spec.read_spec(spectrum,order=order)[20:-20]
        #wave,flux = s[:,0],s[:,1]
        wave = spectrum_hdulist[0].data[order,:,0][20:-20]
        flux = spectrum_hdulist[0].data[order,:,2][20:-20]
        mask = flux != flux
        flux[mask] = 1
        

        if min(wave) < 6270 and max(wave) > 6320:

            mask = wave > 6270
            mask *= wave < 6320


            vel = arange(-300,300,0.01)
            lsd_interp = spectype.make_rot_prof(vel,vsini=vsini,scale=1.,vshift=vshift)
            lsd_interp = spectype.vgauss(vel,lsd_interp,macro=5.66/2.355,scale=1.)
            lsd_interp = transpose(array([vel,lsd_interp]))
            #epos,drv,cc = cross_correlate_order(wave[mask],flux[mask],template,ap=10)
            #print epos

            flux = remove_stellar_template(stellar_template,transpose(array([wave,flux])),lsd_interp)


            ccf = runlsd(transpose(array([wave,flux])),template)
            ccfmask = abs(ccf[:,0]) < 60
            ccf = ccf[ccfmask]
            ccfmask = abs(ccf[:,0]) > 30
            fit = polyfit(ccf[:,0][ccfmask],ccf[:,1][ccfmask],1)
            fit = polyval(fit,ccf[:,0])
            ccf[:,1] -= fit

            # plt.plot(ccf[:,0],ccf[:,1])
            # plt.show()
            
            x0 = spectype.fitlsd_mcmc(ccf,vsini_init=12.0,nrun=500)
            epos = x0[2]


            
            #runlsd(transpose(array([wave,flux])),template)
            #sys.exit()
            
            #epos,drv,cc = cross_correlate_order(wave,flux,template,ap=10)

            # print "telluric rv",epos
            # plt.plot(drv,cc)
            # plt.axvline(x=epos)
            # plt.show()
            
            break

    return epos



if __name__ == "__main__":

    spectrum = "/data/tfop/chiron_data/reduced/CHIRON_2019-07-05T01-52-45.0_T0193641523.fits"
    template = "spectral_library/template_6000_4.0_0.0.dat"

    #rvlist,rv = main(spectrum,template)

    #rvlist,rv = main(spectrum,template)
    print measure_telluric_rv(spectrum,template,6.571710801257795, 44.384706922190816)

    

    
    
