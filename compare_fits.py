"""
This script compares FITS files in order to see at which point they differ.
The file 'set_values.cfg' is the config file for this application and has examples
of how to format values appropriately. If any bugs are found or any alterations
are needed, please contact Jesse Averbukh at javerbukh@stsci.edu
"""

import ConfigParser
import argparse
from astropy.io import fits
import numpy as np
import sys


def get_config_file_names(config,args):
    """
    Reads in the config file specified in the command line arguement used to
    run this script
    """
    config.read(args.chosen_config)
    all_file_pairs = []

    for section in config.sections():
        output_file = config.get(section, 'output_file')
        expected_output_file = config.get(section, 'expected_output_file')

        file_pair = output_file, expected_output_file
        all_file_pairs.append(file_pair)

    return all_file_pairs

def tests_for_size_four(hdulist_out, hdulist_expected, elem):
    if np.allclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06):
        print ("Both files have similar values for {}".format(elem))
    else:
        equivalence_arr = np.isclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06)
        where_false = np.where(equivalence_arr == False)
        print ("The files are unequal at the following pixel locations in {}:".format(elem))
        for x in range(0,len(where_false[0])):
            if x>10:
                print("\tetc...ending forloop for sake of brevity")
                break
            print ("\tIntegration: {}, Group: {}, X,Y: ({},{})".format(where_false[0][x],
            where_false[1][x], where_false[2][x], where_false[3][x]))
            print ("\t\tOutput file value at this location = {}, Expected file = {}".format(hdulist_out[elem].data[where_false[0][x]][where_false[1][x]][where_false[2][x]][where_false[3][x]],
                hdulist_expected[elem].data[where_false[0][x]][where_false[1][x]][where_false[2][x]][where_false[3][x]]))

def tests_for_size_three(hdulist_out, hdulist_expected, elem):
    if np.allclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06):
        print ("Both files have similar values for {}".format(elem))
    else:
        equivalence_arr = np.isclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06)
        where_false = np.where(equivalence_arr == False)
        print ("The files are unequal at the following pixel locations in {}:".format(elem))
        for x in range(0,len(where_false[0])):
            if x>10:
                print("\tetc...ending forloop for sake of brevity")
                break
            print ("\tGroup: {}, X,Y: ({},{})".format(where_false[0][x],
            where_false[1][x], where_false[2][x]))
            print ("\t\tOutput file value at this location = {}, Expected file = {}".format(hdulist_out[elem].data[where_false[0][x]][where_false[1][x]][where_false[2][x]],
                hdulist_expected[elem].data[where_false[0][x]][where_false[1][x]][where_false[2][x]]))

def tests_for_size_two(hdulist_out, hdulist_expected, elem):
    if np.allclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06):
        print ("Both files have similar values for {}".format(elem))
    else:
        equivalence_arr = np.isclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06)
        where_false = np.where(equivalence_arr == False)
        print ("The files are unequal at the following pixel locations in {}:".format(elem))
        for x in range(0,len(where_false[0])):
            if x>10:
                print("\tetc...ending forloop for sake of brevity")
                break
            print ("\tX,Y: ({},{})".format(where_false[0][x],
            where_false[1][x]))
            print ("\t\tOutput file value at this location = {}, Expected file = {}".format(hdulist_out[elem].data[where_false[0][x]][where_false[1][x]],
                hdulist_expected[elem].data[where_false[0][x]][where_false[1][x]]))

def tests_for_size_one(hdulist_out, hdulist_expected, elem):
    if np.allclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06):
        print ("Both files have similar values for {}".format(elem))
    else:
        equivalence_arr = np.isclose(hdulist_out[elem].data,hdulist_expected[elem].data,atol=1e-06)
        where_false = np.where(equivalence_arr == False)
        print ("The files are unequal at the following pixel locations in {}:".format(elem))
        for x in range(0,len(where_false[0])):
            if x>10:
                print("\tetc...ending forloop for sake of brevity")
                break
            print ("\tIndex: {}".format(where_false[0][x]))
            print ("\t\tOutput file value at this location = {}, Expected file = {}".format(hdulist_out[elem].data[where_false[0][x]],
                hdulist_expected[elem].data[where_false[0][x]]))

def tests_for_BinTableHDU(hdulist_out, hdulist_expected, elem):
    rows_not_present = []

    for out_row in range(0,len(hdulist_out[elem].data)):
        checker = False
        checker_val = 0
        for expected_row in range(0,len(hdulist_expected[elem].data)):
            if set(hdulist_out[elem].data[out_row]).issubset(set(hdulist_expected[elem].data[expected_row])):
                checker = True
        if not checker:
            rows_not_present.append(hdulist_out[elem].data[out_row])

    if rows_not_present:
        print ("\tThe following rows are not identical:")
        for row in rows_not_present:
            print ("\t\t{}".format(row))
    else:
        print ("\tRows in table are identical")

def even_deeper_tests(output_file, expected_output_file, all_elems_out):
    """
    Checks to see where exactly the files differ on the scale of pixel values within
    all of the elements of hdulist
    """
    hdulist_out = fits.open(output_file)
    hdulist_expected = fits.open(expected_output_file)

    for elem in all_elems_out:
        if hdulist_out[elem].data is None:
            print ("{} has no dimensions to check".format(elem))
            pass
        elif (isinstance(hdulist_out[elem], fits.hdu.table._TableBaseHDU)):
            print ("Checking values for {} inside its BinTable".format(elem))
            tests_for_BinTableHDU(hdulist_out, hdulist_expected, elem)
        elif not (isinstance(hdulist_out[elem], (fits.hdu.image.ImageHDU, fits.hdu.image.PrimaryHDU))):
            print ("{} is not of the appropriate type (ImageHDU or PrimaryHDU)".format(elem))
            pass
        elif hdulist_out[elem].shape != hdulist_expected[elem].shape:
            print ("The two files do not have the same dimensions for {}".format(elem))
            pass
        elif len(hdulist_out[elem].data.shape) == 4:
            tests_for_size_four(hdulist_out, hdulist_expected, elem)
        elif len(hdulist_out[elem].data.shape) == 3:
            tests_for_size_three(hdulist_out, hdulist_expected, elem)
        elif len(hdulist_out[elem].data.shape) == 2:
            tests_for_size_two(hdulist_out, hdulist_expected, elem)
        elif len(hdulist_out[elem].data.shape) == 1:
            tests_for_size_one(hdulist_out, hdulist_expected, elem)
        else:
            print ("Was not able to check {}".format(elem))
    return True

def deeper_tests(output_file, expected_output_file):
    """
    Checks to see if files have the same headers and if they have the same SCI
    dimensions, in that order
    """
    hdulist_out = fits.open(output_file)
    hdulist_expected = fits.open(expected_output_file)
    all_elems_out = []
    for elem_out in range(0,len(hdulist_out)):
        all_elems_out.append(''.join(str(hdulist_out[elem_out].name)))
    print ("{} headers: {}".format(output_file, all_elems_out))

    all_elems_expected = []
    for elem_expected in range(0,len(hdulist_expected)):
        all_elems_expected.append(''.join(str(hdulist_expected[elem_expected].name)))
    print ("{} headers: {}".format(expected_output_file, all_elems_expected))

    if all_elems_out == all_elems_expected:
        print ("Both files have the same headers, now checking dimensions and pixels for each header...\n")
        even_deeper_tests(output_file, expected_output_file, all_elems_out)
        return True
    else:
        print ("Files do not have the same headers")
        return False
    #
    # print ("{} SCI dimensions: {}".format(output_file, hdulist_out['SCI'].shape))
    # print ("{} SCI dimensions: {}".format(expected_output_file, hdulist_expected['SCI'].shape))
    #
    # if hdulist_out['SCI'].shape == hdulist_expected['SCI'].shape:
    #     print ("Both files have the same SCI dimensions, now checking individual pixels...\n")
    # else:
    #     print("Files do not have the same SCI dimensions")
    #     return False

    even_deeper_tests(output_file, expected_output_file, all_elems_out)
    return True

def comparison_tests(all_file_pairs):
    """
    Takes two files for the config file and checks to see if they are identical.
    If not, further tests are run in deeper_tests()
    """
    for (output_file, expected_output_file) in all_file_pairs:
        print ("\nRunning tests for the following files:\nOutput file: {}\nExpected output file: {}\n".format(output_file, expected_output_file))
        hdulist_out = fits.open(output_file)
        hdulist_expected = fits.open(expected_output_file)
        fd = fits.FITSDiff(hdulist_out,hdulist_expected)
        if fd.identical:
            print ("FITSDiff is True, files {} and {} are identical".format(output_file, expected_output_file))
        else:
            print ("FITSDiff returned False, running deeper tests now...\n")
            deeper_tests(output_file, expected_output_file)
        print ("------------------------------------------------------------------------------------------------\n")


parser = argparse.ArgumentParser()
parser.add_argument("chosen_config", nargs = "?", default = "None", help= "The Config file with the files that are to be compared, or file 1 of the files to be compared")
parser.add_argument("chosen_file", nargs='?', default='None', help = "(Optional) file 2 of files to be compared in command line")
args = parser.parse_args()

if len(sys.argv) == 2 and not args.chosen_config == "None":
    config = ConfigParser.ConfigParser()
    all_file_pairs = get_config_file_names(config,args)
    comparison_tests(all_file_pairs)
elif len(sys.argv) == 3 and not args.chosen_file == "None" and not args.chosen_config == "None":
    file_pair = [(args.chosen_config, args.chosen_file)]
    comparison_tests(file_pair)
else:
    print ("Incorrect number of arguements, format is either:")
    print ("\tpython compare_fits.py [chosen_config]")
    print ("\tpython compare_fits.py [file_to_compare1] [file_to_compare2]")
