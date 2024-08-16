from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def flat_corr(num_of_x, num_of_y):
    for img in range(0, len(flat_array)):
        flat_image = 'flat_R.001/' + flat_array[img]
        image_data2 = fits.getdata(flat_image, ext=0)
        aver = np.average(image_data2)

        for i in range(0, num_of_x):
            for j in range(0, num_of_y):
                temp = float(image_data2[i][j])/aver
                image_data2[i][j] = temp

        name = 'flat_R.001/nFlat' + str(img+1) + ".fits"
        fits.writeto(name, image_data2, overwrite=True)

    norm_array = list()
    for n in range(1, 10):
        temp = "nFlat" + str(n) + ".fits"
        norm_array.append(temp)
    
    mflat(xax, yax, norm_array)


def mflat(num_of_x, num_of_y, n_array):
    nflat_fits = "flat_R.001/nFlat1.fits"
    mflat_data = fits.getdata(nflat_fits, ext=0)

    #for i in range(0, num_of_x):
    for i in range(0, 10):
        #print(i)
        #for j in range(0, num_of_y):
        for j in range(0, 10):
            median_arr = list()
            for img in range(0, len(n_array)):
                image = 'flat_R.001/' + n_array[img]
                image_data2 = fits.getdata(image, ext=0)
                median_arr.append(image_data2[i][j])
            med = np.median(median_arr)
            mflat_data[i][j] = med

    Mname = 'flat_R.001/masterFlat.fits'
    fits.writeto(Mname, mflat_data, overwrite=True)

    plt.figure()
    image_data5 = fits.getdata(Mname, ext=0)
    plt.imshow(image_data5, cmap='gray')
    plt.colorbar()
    plt.title("master flat-field")

    plt.figure()
    image_data3 = fits.getdata(nflat_fits, ext=0)
    plt.imshow(image_data3, cmap='gray')
    plt.colorbar()
    plt.title("normalized flat-field")
    
    my_fits2 = fits.open(nflat_fits)
    my_fits2[0].header.append(('COMMENT', 'Exercise 8.3'), end=True)
    fits.writeto(Mname, mflat_data, my_fits2[0].header, overwrite=True)

    plt.show()
    return Mname


def flatApplic(masterFlat):
    light_array = ['46P_3_R-0001.fit', '60383_C-0005_d.fit', 'IC_434_B_30s-002.fit',
                   'M66_R_1-003.fit', 'M74_R_120s-001.fit', 'NGC3628_R_1-002.fit',
                   'POINT_FLIP_360S_R_0002.fit']

    mflat_data = fits.getdata(masterFlat, ext=0)
    aver = np.average(mflat_data)

    #for i in range(0, len(light_array)):
    for i in range(0, 1):
        light = 'LIGHT/' + light_array[i+3]
        light_data = fits.getdata(light, ext=0)

        new = np.divide(light_data, aver)
        print(np.amin(new), np.amax(new))

        plt.imshow(new, cmap='gray', vmin=np.amin(new), vmax=22000*aver)
        plt.title("Application of master flat-field")
        plt.show()


if __name__ == "__main__":
    flat_array = list()
    for flat in range(10, 20):
        temp = "flat_R-00" + str(flat) + ".fit"
        flat_array.append(temp)

    image_data = fits.getdata('flat_R.001/flat_R-0010.fit', ext=0)
   
    xax = image_data.shape[0]
    yax = image_data.shape[1]

    flat_corr(xax, yax)
