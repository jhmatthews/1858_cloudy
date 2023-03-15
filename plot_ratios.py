import numpy as np
import matplotlib.pyplot as plt
import jm_util
jm_util.set_mod_defaults()

o4blends = ["1397.20A", "1399.77A", "1401.16A", "1404.78A", "1407.38A"]
nvblends = ["1239.93A", "1240.40A"]
# do you want to include Mg II and OIV lines for NV and Si IV as Mauche does
INCLUDE_BLENDS = True
plot_without_blends_too = True


def get_cloudy_line_strengths(BB_TEMP, LOG_N, logxi, folder="cld_data/", suffix="",
                              line_names=["c4", "n5", "he2", "si4", "o4", "o4sum", "mg2sum", "si4_raw", "n5_raw"]):
    """Get cloudy line strengths from a given file for a range of IPs

    Args:
            BB_TEMP (float, optional): black body temperature in eV. 
            LOG_N (float, optional): log electron density. 
            logxi (array-like): logarithms of ionization parameters
            folder (str, optional): data folder name. Defaults to "cld_data/".
            suffix (str, optional): CNO or not? ==cno if so, blank otherwise
            line_names (list, optional): list of line names. Defaults to ["c4", "n5","he2", "si4"].

    Returns:
            cld_data (dict): dictionary of cloudy data
    """
    cld_data = dict()

    # create arrays to hold line luminosities
    for l in line_names:
        cld_data[l] = np.zeros_like(logxi)

    for i, logxi in enumerate(logxis):
        try:
            root = root = "{}/BB_{:.1f}_1858_n{:.1f}_xi{:.3f}{}".format(
                folder, BB_TEMP, LOG_N, logxi, suffix)
            f = open(root + ".lines", "r")
            for line in f:
                data = line.split()
                if data[1] == "Blnd" and data[2] == "1397.00A":
                    cld_data["si4"][i] = 10.0**float(data[4])
                if data[1] == "Blnd" and data[2] == "1549.00A":
                    cld_data["c4"][i] = 10.0**float(data[4])
                if data[1] == "Blnd" and data[2] == "1240.00A":
                    cld_data["n5"][i] = 10.0**float(data[4])
                if data[1] == "He" and data[2] == "2" and data[3] == "1640.43A":
                    cld_data["he2"][i] = 10.0**float(data[5])
                if data[1] == "Blnd" and data[2] == "1402.00A":
                    cld_data["o4"][i] = 10.0**float(data[4])

                # do you want to include Mg II and O IV blended lines as contributions?
                if INCLUDE_BLENDS:
                    if data[1] == "Mg" and data[2] == "2" and data[3] in nvblends:
                        cld_data["mg2sum"][i] += 10.0**float(data[5])
                    if data[1] == "O" and data[2] == "4" and data[3] in o4blends:
                        cld_data["o4sum"][i] += 10.0**float(data[5])

        except FileNotFoundError:
            print("no file for log xi of ", logxi)

        # store lines without the blends
        cld_data["n5_raw"][i] = cld_data["n5"][i]
        cld_data["si4_raw"][i] = cld_data["si4"][i]

        # include blended lines a la Mauche
        cld_data["si4"][i] += cld_data["o4sum"][i]
        cld_data["n5"][i] += cld_data["mg2sum"][i]

        # store them as log quantities
        for l in line_names:
            cld_data[l][i] = np.log10(cld_data[l][i])

        print(root, cld_data["o4"][i], cld_data["o4sum"]
              [i], cld_data["n5"][i], cld_data["c4"][i])
        # if
        f.close()

    return (cld_data)


def make_figure(logxis, folder="cld_data/", BB_TEMP=50.0, LOG_N=9.0):
    """Make the figures for comparison with Mauche+ 1997

    Args:
            logxis (_type_): log ionization parameters
            folder (str, optional): folder that contains the cloudy data. Defaults to "cld_data/".
            BB_TEMP (float, optional): black body temperature in eV. Defaults to 50.0.
            LOG_N (float, optional): log electron density. Defaults to 9.0.
    """

    # get data
    cld_data = get_cloudy_line_strengths(BB_TEMP, LOG_N, logxis, folder=folder)
    cld_data_cno = get_cloudy_line_strengths(
        BB_TEMP, LOG_N, logxis, folder=folder, suffix="cno")

    # array used to plot a 1:1 line
    logratios = np.linspace(-4, 4, 100)

    # read in the Mauche data for comparison
    lxi_mauche, si4ratio_mauche, n5ratio_mauche, he2ratio_mauche = np.genfromtxt(
        "mauche_data/mauche_data_{:.0f}eV.dat".format(BB_TEMP), unpack=True)

    # kwargs for plotting
    mauche_kwargs = {"marker": "x", "linestyle": "--",
                     "color": "C1", "label": "Mauche"}
    cno_kwargs = {"marker": None, "linestyle": "-",
                  "color": "C2", "label": "Cloudy CNO", "lw": 3}
    raw_kwargs = {"marker": None, "linestyle": "-", "alpha": 0.3,
                  "color": "C2", "label": "Cloudy CNO (no blends)", "lw": 4}
    std_kwargs = {"marker": None, "linestyle": "-",
                  "color": "C0", "label": "Cloudy", "lw": 3}
    save_kwargs = {"dpi": 300}

    sed_string = "bb{:.0f}".format(BB_TEMP)
    title_string = "Blackbody, {:.0f} eV".format(BB_TEMP)

    # make some figures
    plt.figure()
    plt.title(title_string)
    plt.plot(logxis, cld_data["si4"]-cld_data["c4"], **std_kwargs)
    plt.plot(logxis, cld_data_cno["si4"]-cld_data_cno["c4"], **cno_kwargs)

    if plot_without_blends_too:
        plt.plot(logxis, cld_data_cno["si4_raw"]-cld_data_cno["c4"], **raw_kwargs)
    plt.plot(lxi_mauche, si4ratio_mauche, **mauche_kwargs)
    plt.xlabel(r"$\log \xi$", fontsize=18)
    plt.ylabel(r"$\log ({\rm Si IV / CIV})$", fontsize=18)
    plt.legend()
    plt.tight_layout()
    # plt.xlim(-3,2)
    plt.ylim(-3,2)
    plt.savefig("mauche_cld_comp_si4c4_{}.png".format(
        sed_string), **save_kwargs)

    plt.figure()
    plt.title(title_string)
    plt.plot(cld_data["si4"]-cld_data["c4"],
             cld_data["n5"]-cld_data["c4"], **std_kwargs)
    plt.plot(cld_data_cno["si4"]-cld_data_cno["c4"],
             cld_data_cno["n5"]-cld_data_cno["c4"], **cno_kwargs)
    if plot_without_blends_too:
        plt.plot(cld_data_cno["si4_raw"]-cld_data_cno["c4"],
             cld_data_cno["n5_raw"]-cld_data_cno["c4"], **raw_kwargs)
    plt.plot(si4ratio_mauche, n5ratio_mauche, **mauche_kwargs)
    plt.ylabel(r"$\log ({\rm N V / CIV})$", fontsize=18)
    plt.xlabel(r"$\log ({\rm Si IV / CIV})$", fontsize=18)
    plt.plot(logratios, logratios, ls="--", c="k")
    plt.scatter(1.1, 1.8, c="C3", marker="*", s=200)
    plt.legend()
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.tight_layout()
    plt.savefig("mauche_cld_comp_ratio1_{}.png".format(
        sed_string), **save_kwargs)

    plt.figure()
    plt.title(title_string)
    plt.plot(cld_data["si4"]-cld_data["c4"],
             cld_data["he2"]-cld_data["c4"], **std_kwargs)
    plt.plot(cld_data_cno["si4"]-cld_data_cno["c4"],
             cld_data_cno["he2"]-cld_data_cno["c4"], **cno_kwargs)
    if plot_without_blends_too:
        plt.plot(cld_data_cno["si4_raw"]-cld_data_cno["c4"],
             cld_data_cno["he2"]-cld_data_cno["c4"], **raw_kwargs)
    plt.plot(si4ratio_mauche, he2ratio_mauche, **mauche_kwargs)
    plt.plot(logratios, logratios, ls="--", c="k")
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.scatter(1.1, 1.1, c="C3", marker="*", s=200)
    plt.ylabel(r"$\log ({\rm He II / CIV})$", fontsize=18)
    plt.xlabel(r"$\log ({\rm Si IV / CIV})$", fontsize=18)
    plt.legend()
    plt.tight_layout()
    plt.savefig("mauche_cld_comp_ratio2_{}.png".format(
        sed_string), **save_kwargs)


if __name__ == "__main__":
    logxis = np.arange(-3, 2, 0.125)
    make_figure(logxis, BB_TEMP=10.0)
    make_figure(logxis, BB_TEMP=30.0)
    make_figure(logxis, BB_TEMP=50.0)
