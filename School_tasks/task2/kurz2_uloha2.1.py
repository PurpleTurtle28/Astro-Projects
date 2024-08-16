import matplotlib.pyplot as plt

def graf(my_file):
    name = my_file[5:]
    wavelength = list()
    intensity = list()

    with open(my_file, "r") as in_file:
        for line in in_file:
            w, i = line.split()
            wavelength.append(float(w))
            intensity.append(float(i))

    plt.plot(wavelength, intensity, c='g')
    plt.xlabel("wavelength [nm]")
    plt.ylabel("relative intensity [ADU]")

    plt.title(name)
    plt.show()

graf("data/M20161214_002620_PGRACAM-LP_wsum.txt")
