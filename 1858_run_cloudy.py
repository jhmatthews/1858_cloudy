import os
import numpy as np
import convert_ip
import argparse
from tqdm import tqdm

# some constants
# CAUTION: COLUMN DENSITY IS HARDWIRED AT THE MOMENT
BOLTZMANN = 1.38062e-16
EV2ERGS = 1.602192e-12
COLUMN = 23
OUTPUT_FOLDER = "cld_output"

# Different strings to control abundances
CNO_STRING = """
abundances H II region
element scale factor helium 2
element scale factor carbon 0.5
element scale factor nitrogen 7"""

CNO_STRING_EQ = """
abundances H II region
element scale factor carbon 0.018
element scale factor nitrogen 7.5
element scale factor oxygen 0.037"""


def run_cloudy(fname, cloudy_path="/Users/matthewsj/software/c17.01/source/", folder="cld_data/"):
    isys = os.system(
        "{}/cloudy.exe < {}/{}.in > {}/{}.out".format(cloudy_path, folder, fname, folder, fname))
    return (isys)


def get_param_string(root, log_xi, log_U, sed="bb", temp=None,
                     abundance="solar", LOG_N=9, folder="cld_data"):
    '''
    Create a string to write to a cloudy input file. 

    Parameters:
            root 	str
                            root filename to use for naming

            log_xi  float 
                            log of xi ionization parameter (only used for title)

            log_U 	float 
                            log of U ionization parameter

            temp float 
                            temperature in K (linear)

    Returns:
            my_string str
            string containing text to use for input file
    '''
    if abundance == "solar":
        abundance_string = "abundances solar"
    elif abundance == "cno":
        abundance_string = CNO_STRING
    elif abundance == "cno_eq":
        abundance_string = CNO_STRING_EQ

    if sed == "bb":
        sed_string = "Blackbody {}".format(temp)
    elif sed == "brems":
        sed_string = "brems {}".format(temp)

    my_string = '''title {} model, log density {:.1f}, log xi {:.1f}
# =========
# commands controlling continuum =========
{}
ionization parameter {:.1f}
# =========
# commands for density & abundances
# =========
hden {:.1f}{}
# =========
# commands controlling geometry  
# =========
stop column density 19
# =========
# other commands for details
# =========
iterate convergence 
# these are to try to speed things up
init "c84.ini"
no molecules
no level2 lines
no fine opacities
atom h-like levels small
atom he-like levels small
COSMIC RAY BACKGROUND
element limit off -8
# =========
# commands controlling output    =========
# =========
normalize to "H  1" 1215.67A 100 
# print line faint 1  # get rid of lots of faint lines c
# print off
print lines column
save lines, array "{}/{}.lines" last
save continuum "{}/{}.cont" 
#
# ========================================
'''.format(sed, LOG_N, log_xi, sed_string, log_U, LOG_N, abundance_string, folder, root, folder, root)
    return my_string


def initialise_cloudy_sim(root, log_xi, log_U, LOG_N=9, sed="bb", temp=580225.906078, abundance="solar", folder="cld_data/"):
    """Initialise a cloudy simulation 

    Args:
        root (_type_): root filename to use for naming
        log_xi (_type_): log of xi ionization parameter (only used for title)
        log_U (_type_): log of U ionization parameter
        LOG_N (int, optional): log electron density. Defaults to 9.
        sed (str, optional): SED type to use: must be bb or brems.. Defaults to "bb".
        temp (float, optional): Temperature in K (linear). Defaults to 580225.906078.
        abundance (str, optional): abundances to use. . Defaults to "solar".
        folder (str, optional): folder to store data. Defaults to "cld_data/".
    """

    # get the string to write to file
    params = get_param_string(
        root, log_xi, log_U, LOG_N=LOG_N, temp=temp, sed=sed, abundance=abundance, folder=folder)

    # chuck the text in an input file
    fname = "{}/{}.in".format(folder, root)
    f = open(fname, "w")
    f.write(params)
    f.close()


def run_logxi_grid(logxis=np.arange(-3, 2, 0.5), sed="bb", TEMP=50.0, LOG_N=9, abundance="solar", folder="cld_data/"):
    """Run a grid of blackbody cloudy simulations 

    Args:
        logxis (array-like, optional): array of floats containing log10 of ionization parameter xi. Defaults to np.arange(-3, 2, 0.5).
        sed (str, optional): SED type to use: must be bb or brems. Defaults to "bb".
        TEMP (float, optional): Blackbody temperature in ELECTRON VOLTS (linear). Defaults to 50.0.
        LOG_N (int, optional): log of density. Defaults to 9.
        abundance (str, optional): abundances to use. Defaults to "solar".
        folder (str, optional): folder to save the data. Defaults to "cld_data/".
    """
    for i, logxi in enumerate(tqdm(logxis)):
        #print (i, logxi)

        # we need to calculate U, the IP used by Cloudy.
        xi = 10.0 ** logxi
        if sed == "bb":
            factor = convert_ip.get_conversion_factor_blackbody(temp_eV=TEMP)
            root = "BB_{:.1f}_1858_n{:.1f}_xi{:.3f}".format(TEMP, LOG_N, logxi)
            # root filename to use for input and output files
        elif sed == "brems":
            factor = convert_ip.get_conversion_factor_brems(
                temp_eV=TEMP * 1000.0)
            root = "brems_{:.1f}_1858_n{:.1f}_xi{:.3f}".format(
                TEMP, LOG_N, logxi)
        U = xi * factor
        logU = np.log10(U)

        if "cno" in abundance:
            root += abundance

        #print (BB_TEMP * EV2ERGS / BOLTZMANN)
        if sed == "bb":
            temp = TEMP * EV2ERGS / BOLTZMANN
        elif sed == "brems":
            temp = np.log10(TEMP * 1000.0 * EV2ERGS / BOLTZMANN)
        # Initialise simulation by creating input file
        initialise_cloudy_sim(root, logxi, logU, LOG_N=LOG_N, temp=temp,
                              abundance=abundance, folder=folder, sed=sed)

        # run the sim!
        run_cloudy(root, folder=folder)


def parse_args():
    """parse command line arguments

    Returns:
        args: arguments to use in program
    """
    parser = argparse.ArgumentParser(
        prog='1858_run_cloudy',
        description='Run a cloudy grid for the 1858 abundances project')

    parser.add_argument('-t', '--temp', nargs=1, type=float,
                        help="Temperature in eV, or keV if brems",
                        default=50.0)

    parser.add_argument('-sed', '--sed', nargs=1, type=str,
                        help="SED (must be bb or brems)",
                        default="bb", choices={"bb", "brems"})

    parser.add_argument('-a', '--abun', nargs=1, type=str,
                        help="which abundances to use. cno_eq means CNO equilibrium.",
                        default="solar", choices={"solar", "cno", "cno_eq"})

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_args()
    TEMP = args.temp
    abundance = args.abun
    sed = args.sed

    # density hardwired to 1e9
    run_logxi_grid(logxis=np.arange(-3, 2, 0.125), TEMP=TEMP,
                   LOG_N=9, abundance=abundance, sed=sed)
