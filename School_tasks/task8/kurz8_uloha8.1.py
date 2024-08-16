from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def visualize(n):
    img = 'LIGHT/' + image_array[n]
    image_data = fits.getdata(img, ext=0)

    mean = np.mean(image_data)
    sigma = np.std(image_data)

    fig = plt.figure(figsize=(8, 8))
    ax = fig.subplots(2, 2)

    h = 1
    for i in range(0, 2):
        for j in range(0, 2):
            ax[i, j].imshow(image_data, cmap='gray', vmin=np.amin(image_data), vmax=mean+h*sigma)
            h += 1
    plt.suptitle(str(image_array[n]) + "\nfor sigma 1 - 4")


if __name__ == "__main__":
    image_array = ['46P_3_R-0001.fit', '60383_C-0005_d.fit', 'IC_434_B_30s-002.fit',
                   'M66_R_1-003.fit', 'M74_R_120s-001.fit', 'NGC3628_R_1-002.fit',
                   'POINT_FLIP_360S_R_0002.fit']

    visualize(3)
    visualize(5)
    visualize(6)

    plt.show()
