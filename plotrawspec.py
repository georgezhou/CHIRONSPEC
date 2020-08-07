import os,sys,string,pyfits,pickle
from numpy import *
import matplotlib.pyplot as plt


if __name__ == "__main__":
    specpkl = sys.argv[1]

    filename, file_extension = os.path.splitext(specpkl)

    if file_extension == ".pkl":

        try:
            spectrum,background,tharspec = pickle.load(open(specpkl,"rb"))
        except ValueError:
            spectrum,background,tharspec,spectrum_noflat,background_noflat = pickle.load(open(specpkl,"rb"))

        for order in range(len(spectrum)):
            plt.subplot(311)
            plt.title(str(order))
            plt.plot(spectrum_noflat[order])
            plt.plot(background_noflat[order])

            plt.subplot(312)
            plt.plot(spectrum[order])
            plt.plot(background[order])

            plt.subplot(313)
            plt.plot(tharspec[order])

            plt.show()

    if file_extension == ".fits":
        
        spectrum = pyfits.open(specpkl)
        spectrum = spectrum[0].data


        for order in range(len(spectrum)):
            wave = spectrum[order,:,0]
            flux = spectrum[order,:,2]
            plt.title(str(order))
            plt.plot(wave,flux)
            plt.show()
