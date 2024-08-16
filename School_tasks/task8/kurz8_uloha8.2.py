from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def clean(n):
    img = 'LIGHT/' + image_array[n]
    my_fits = fits.open(img)
    image_data = fits.getdata(img, ext=0)

    xax = my_fits[0].header["NAXIS1"]
    yax = my_fits[0].header["NAXIS2"]

    change = np.median(image_data)
    counter = 0
    for i in range(0, xax):
        for j in range(0, yax):
            if(image_data[i][j]) == 0 or (image_data[i][j]) == 65535:
                image_data[i][j] = change
                counter += 1

    print("Number of replaced pixels: " + str(counter))

    aver = np.average(image_data)
    plt.figure()
    plt.imshow(image_data, cmap='gray', vmin=np.amin(image_data), vmax=2*aver)
    plt.colorbar()
    plt.title(image_array[n])

    header_title = image_array[n] + "_new.fits"
    my_fits[0].header.append(('COMMENT', 'Exercise 8.2'), end=True)
    fits.writeto(header_title, image_data, my_fits[0].header, overwrite=True)


if __name__ == "__main__":
    image_array = ['46P_3_R-0001.fit', '60383_C-0005_d.fit', 'IC_434_B_30s-002.fit',
                'M66_R_1-003.fit', 'M74_R_120s-001.fit', 'NGC3628_R_1-002.fit',
                'POINT_FLIP_360S_R_0002.fit']

    clean(0)
    clean(3)
    clean(4)

    plt.show()
