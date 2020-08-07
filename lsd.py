import os,sys,string
from numpy import *
import matplotlib.pyplot as plt
import matplotlib

matplotlib.pyplot.switch_backend('agg')
from scipy import interpolate,optimize
import scipy.linalg as scilinalg
import pyfits
import pickle

c = 3.*10**5

    

def binspec(spec,binstep = 0.04):
    newwave = arange(min(spec[:,0]),max(spec[:,0]),binstep)
    newspec = ones(len(newwave))
    for i in range(len(newwave)):
        mask = abs(spec[:,0]-newwave[i]) < binstep/2.
        newspec[i] = mean(spec[:,1][mask])

    newspec = transpose(array([newwave,newspec]))
    return newspec

def normalise(spec,niter=1,sigma_low = 0.05,deg=5):

    ### normalise the spectrum to 1

    
    x = arange(len(spec))
    mask = spec == spec
    spec_iter = spec[mask].copy()
    x_iter = x[mask].copy()

    


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





def make_rot_prof(v,vsini=100,epsilon=0.5,scale=1,vshift=0):
    
    vrot = 2*(1-epsilon)*(1-((v-vshift)/vsini)**2)**(0.5) + 0.5*pi*epsilon*(1-((v-vshift)/vsini)**2)
    #vrot /= pi*vsini*(1-epsilon/3)


    mask = vrot != vrot
    vrot[mask] = 0

    vrot /= max(vrot)
    vrot *= scale
    
    return vrot


def vrot_minfunc(x0,vel,ccf):
    ### x0 = [vsini,epsilon,scale,vshift]
    #print x0

    good_x0 = True
    if x0[0] < 0 or x0[0] > 250:
        good_x0 = False
    if x0[1] < 0 or x0[1] > 1:
        good_x0 = False
    if x0[2] < 0 or x0[2] > 100.:
        good_x0 = False
    if abs(x0[3]) > 200:
        good_x0 = False

    if good_x0:
        vrot = make_rot_prof(vel,vsini=x0[0],epsilon=x0[1],scale=x0[2],vshift=x0[3])
        diff = sum((vrot-ccf)**2)
        return diff
    else:
        return nan

def fit_ccf(ccf,vel):

    maxheight = max(ccf)
    maxloc = vel[argmax(ccf)]

    x0 = []
    for i in range(10):
        x0_i = [random.normal(10,1.),random.normal(0.5,0.1),random.normal(maxheight,0.1*maxheight),random.normal(maxloc,5.)]
        x0_i = optimize.fmin(vrot_minfunc,x0_i,args=(vel,ccf),disp=0)
        
        x0.append(x0_i)

    x0 = array(x0)
    x0 = median(x0,axis=0)
        
    
    return x0





def convmtx(v, n):
    """Generates a convolution matrix
    
    Usage: X = convm(v,n)
    Given a vector v of length N, an N+n-1 by n convolution matrix is
    generated of the following form:
              |  v(0)  0      0     ...      0    |
              |  v(1) v(0)    0     ...      0    |
              |  v(2) v(1)   v(0)   ...      0    |
         X =  |   .    .      .              .    |
              |   .    .      .              .    |
              |   .    .      .              .    |
              |  v(N) v(N-1) v(N-2) ...  v(N-n+1) |
              |   0   v(N)   v(N-1) ...  v(N-n+2) |
              |   .    .      .              .    |
              |   .    .      .              .    |
              |   0    0      0     ...    v(N)   |
    And then it's trasposed to fit the MATLAB return value.     
    That is, v is assumed to be causal, and zero-valued after N.
    """
    N = len(v) + 2*n - 2
    xpad = concatenate([zeros(n-1), v[:], zeros(n-1)])
    X = zeros((len(v)+n-1, n))
    # Construct X column by column
    for i in xrange(n):
        X[:,i] = xpad[n-i-1:N-i]
    
    return X.transpose()

# def apodize(spectrum):
#     aplen = 0.2*len(spectrum)
#     spectrum = 1-spectrum
#     for i in range(len(spectrum)):
#         if i < aplen:
#             spectrum[i] *= float(i)/float(aplen)

#         if i > len(spectrum)-aplen:
#             spectrum[i] *= float(len(spectrum)-i)/float(aplen)

#     spectrum = 1-spectrum

#     # plt.plot(spectrum)
#     # plt.show()
#     return spectrum


def apodize(x,ap=0.2):
    func = ones(len(x))
    ap = int(ap*len(x))
    xpos = arange(len(x))

    mask = xpos<ap
    func[mask] = 0.5*cos((xpos[mask]+ap)*2*pi/(ap*2))+0.5
    
    mask = max(xpos)-xpos < ap
    func[mask] = 0.5*cos((max(xpos)-ap-xpos[mask])*2*pi/(ap*2))+0.5

    x = 1-x
    x = func*x
    x = 1-x

    return x


def interpolate_spectra(spectrum,template):
    ### interpolate both spectra to gridded velocity space

    oversample = 2.0


    ### convert wavelength axis to velocity
    if len(spectrum)%2==0:
        spectrum = spectrum[1:]
   
    linecentre = spectrum[int(round(float(len(spectrum))/2)),0]
    vel = -1*c*(linecentre-spectrum[:,0])/spectrum[:,0]
    vel_interp = linspace(min(vel),max(vel),len(vel)*oversample)
    vdel = (max(vel)-min(vel))/(float(len(vel))*oversample)

    #plt.plot(vel,spectrum[:,1])

    interp = interpolate.splrep(vel,spectrum[:,1],k=1)
    spectrum_interp = interpolate.splev(vel_interp,interp)

    vel = -1*c*(linecentre-template[:,0])/template[:,0]
    interp = interpolate.splrep(vel,template[:,1],k=1)
    template_interp = interpolate.splev(vel_interp,interp)


    mask = template_interp < 0
    mask += template_interp > 1
    mask += spectrum_interp < 0
    mask += spectrum_interp > 2.0
    mask = invert(mask)
    vel_interp,spectrum_interp,template_interp = vel_interp[mask],spectrum_interp[mask],template_interp[mask]

    # plt.plot(vel_interp,spectrum_interp)
    # plt.plot(vel_interp,template_interp)
    # plt.show()
    # sys.exit()
         


    return spectrum_interp,template_interp,vdel

def lsd(template,spectrum,vdel):
    len_ccf = 1200

    spectrum_pad = list(zeros(len_ccf/2-1))+list(1-spectrum)+list(zeros(len_ccf/2))
    spectrum_pad = array(spectrum_pad)

    #spectrum_pad = spectrum.copy()
    #template = template[len_ccf/2-1:-1*len_ccf/2]


    vel = arange(-1*len_ccf/2,len_ccf/2,1.)
    vel *= vdel

    M = transpose(convmtx(1-template,len_ccf))
    ccf = scilinalg.lstsq(M,spectrum_pad)[0]
    
    return vel,ccf


def lsd_all(spectrum_array,template_array):
    
    ccf_list = []

    wavemin = 4200
    wavemax = 6200
    


    for order in range(len(spectrum_array)):

        spectrum = spectrum_array[order]
        
        if median(spectrum[:,0]) > wavemin and median(spectrum[:,0]) < wavemax:
            indx = argsort(spectrum[:,0])
            spectrum = spectrum[indx]

            #plt.plot(spectrum[:,0],spectrum[:,1])
            if spectrum[1,0]-spectrum[0,0] < 0.04:
                spectrum = binspec(spectrum,binstep = 0.04)
            #plt.plot(spectrum[:,0],spectrum[:,1])
            #plt.show()
            spectrum[:,1] = normalise(spectrum[:,1],deg=1)
            spectrum = spectrum[spectrum[:,1]==spectrum[:,1]]
            mask = template_array[:,0] > min(spectrum[:,0])-50
            mask *= template_array[:,0] < max(spectrum[:,0])+50
            template = template_array[mask]

            
            ### Now fit template continuum to spectrum
            fit = interpolate.splrep(template[:,0],template[:,1])
            fit = interpolate.splev(spectrum[:,0],fit)

            diff = spectrum[:,1]/fit

            mask = diff-nanmedian(diff) < 0.2
            mask *= diff-nanmedian(diff) > -0.6

            zpt = nanmedian(spectrum[:,0][mask])
            
            fit = polyfit(spectrum[:,0][mask]-zpt,diff[mask],2)
            fitval = polyval(fit,spectrum[:,0]-zpt)

            fit = polyval(fit,spectrum[:,0]-zpt)
            spectrum[:,1] /= fit

            ### Now do lsd 10 times and create an average lsd profile\
            wave0 = min(spectrum[:,0])
            wave1 = max(spectrum[:,0])
            dwave = (wave1-wave0)/float(len(spectrum))
            spectrum_interp,template_interp,vdel=interpolate_spectra(spectrum,template)

            vel,ccf = lsd(template_interp,spectrum_interp,vdel)
            lsd_prof_list = []

            interp = interpolate.splrep(spectrum[:,0],spectrum[:,1],k=1)

            for i in range(10):
                wave_i = arange(random.normal(wave0,2.),random.normal(wave1,2.),dwave)
                spectrum_i = interpolate.splev(wave_i,interp)
                spectrum_i = transpose(array([wave_i,spectrum_i]))
                spectrum_interp,template_interp,vdel=interpolate_spectra(spectrum_i,template)

                spectrum_interp = apodize(spectrum_interp)
                template_interp = apodize(template_interp)
                vel_i,ccf_i = lsd(template_interp,spectrum_interp,vdel)
                
                ccf_i = interpolate.splrep(vel_i,ccf_i,k=1)
                ccf_i = interpolate.splev(vel,ccf_i)

                lsd_prof_list.append(ccf_i)

            ccf = median(array(lsd_prof_list),axis=0)

            # plt.subplot(211)
            # plt.plot(spectrum_interp)
            # plt.plot(template_interp)
            # plt.subplot(212)
            # plt.plot(vel,ccf)
            # plt.show()


            ccf_list.append([vel,ccf])


    return ccf_list


def average_ccf(ccf_list):

    hydrogen_lines = []

    vgrid = arange(-400,400,0.05)

    ### average once to get a general ccf
    temp_ccf = []

    for ccf in ccf_list:

        interp = interpolate.splrep(ccf[0],ccf[1],k=1)
        interp = interpolate.splev(vgrid,interp)
        temp_ccf.append(interp)
  
    temp_ccf = array(temp_ccf)
    temp_ccf = median(temp_ccf,axis=0)


    x0 = fit_ccf(temp_ccf,vgrid)

    vsini = x0[0]
    shift = x0[3]
    scale = x0[2]

    print "best fit vsini,shift,scale",vsini,shift,scale

    good_ccf = []
    weights = []

    order = 1
    for ccf in ccf_list:

        mask_width = vsini/2.
        if mask_width < 3:
            mask_width = 3
        

        mask = ccf[0]-shift > -1*mask_width
        mask *= ccf[0]-shift < mask_width

        out_mask1 = ccf[0]-shift > vsini
        out_mask2 = ccf[0]-shift < -1*vsini
        out_mask = out_mask1 + out_mask2        

        ccf[1] /= median(ccf[1][mask])

        interp = interpolate.splrep(ccf[0],ccf[1],k=1)
        interp = interpolate.splev(vgrid,interp)
        good_ccf.append(interp)
        weights.append(std(ccf[1][out_mask]))
            

        order += 1
   
    good_ccf = array(good_ccf)
    weights = 1/array(weights)**2
    weights = weights / sum(weights)

    ccf_final = zeros(len(good_ccf[0]))
    for i in range(len(good_ccf)):
        ccf_final += good_ccf[i]*weights[i]


    ccf_o = transpose(array([vgrid,ccf_final]))

    return ccf_o,vsini,shift


def run_spectrum(spectrum_file,template_file):
    
    ### format the spectrum file
    spectrum_hdulist = pyfits.open(spectrum_file)
    template_array = loadtxt(template_file)

    spectrum_array = []
    for order in range(0,len(spectrum_hdulist[0].data)):
        wave = spectrum_hdulist[0].data[order,:,0][20:-20]
        flux = spectrum_hdulist[0].data[order,:,2][20:-20]
        spectrum_array.append(transpose(array([wave,flux])))

    ccf_list = lsd_all(spectrum_array,template_array)
    ccf_master,vsini,shift = average_ccf(ccf_list)


    spectrum_name = os.path.basename(spectrum_file)
    spectrum_path = string.replace(spectrum_file,spectrum_name,"")

    pickle.dump([ccf_list,ccf_master,vsini,shift],open(spectrum_path+"lsd_"+spectrum_name+".pkl","wb"))

    
    return ccf_master,vsini,shift
                              

if __name__ == "__main__":

    spectrum_file = "/media/Onion/Data/Chiron/data/reduced/CHIRON_2019-04-22T02-03-30.7_T0123482865.fits"
    template_file = "/media/Onion/Data/spectral_library/lib/template_5500_4.0_0.0.dat"

    ccf_master,vsini,shift = run_spectrum(spectrum_file,template_file)
    print vsini,shift

    plt.plot(ccf_master[:,0],ccf_master[:,1])
    plt.show()
