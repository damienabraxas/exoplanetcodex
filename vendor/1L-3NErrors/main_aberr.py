#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""NLTE corrections calculator

This code takes stellar parameters and line parameters for a specific line
and calculates the NLTE corrections, based on results by Amarsi+2022. 
A neural network is used for interpolation. For more information about
how to prepare the input files and run the code, see the README documentation. 

Visit Github site:
https://github.com/sliljegren/NLTE_corrections

Reach out directly at:
Sofie Liljegren <sofie.liljegren@astro.su.se>

"""

import sys
import numpy as np
import pandas as pd

from function_aberr import *

'''


'''



# def main():
#
# input arguments defined in the terminal
input_args = sys.argv

#defalt 2 input arguments
if len(input_args) == 3:
    input_fname = input_args[1]
    output_fname = input_args[2]
    extrapolate = False

# special case 3 input argument, to turn on extrapolation
elif len(input_args) == 4:
    input_fname = input_args[1]
    output_fname = input_args[2]
    extrapolate = input_args[3]

# wrong number of arguments
else:
    # error message returned
    print(
        '''
        The number of inputs are wrong. Please see README for
        information how to correctly run this script. 
        ---'''
    )
    # program terminated
    sys.exit()

#headers for the data frame
header2 = ['teff/K',
 'lg(g/cms^-2)',
 'A(Fe)',
 'vmic/kms^-1',
 'Elo/eV',
 'Eup',
 'lggf',
 'Species']

# reads the line parameter from the defined input file
input_parameters = pd.read_csv(input_fname, names = header2, comment = '#')

# remove extra spaces in the species input
input_parameters['Species'] = input_parameters['Species'].str.strip()

# define an output data frame
outp = input_parameters.copy()

# add default aberr to the data frame if its not calculated by the script
outp['aberr'] = np.ones(len(input_parameters.index)) * -999

# add extrapolation and in_grid column to output df
outp['extrapolate'] = [extrapolate for i in range(len(input_parameters.index))]
outp['in_grid'] = [False for i in range(len(input_parameters.index))]


# checks if the values in the input are within the boundaries of the original model grid
in_tgrid = (outp['teff/K'] >= 5000) & (outp['teff/K'] <= 6500)
in_lgggrid = (outp['lg(g/cms^-2)'] >= 4) & (outp['lg(g/cms^-2)'] <= 4.5)
in_vmicgrid = (outp['vmic/kms^-1'] >= 0) & (outp['vmic/kms^-1'] <= 3)
in_Afegrid = (outp['A(Fe)'] >= 4.5) & (outp['A(Fe)'] <= 7.5)

indi_e = in_tgrid & in_lgggrid & in_vmicgrid & in_Afegrid

# add True to in_grid column in output df if this is the case
outp.loc[indi_e, 'in_grid'] = True

# if extrapolate, calculates aberr for values outside the original grid
if extrapolate:
    indi_e = outp['extrapolate']

# checks for Fe1 and elo < 2eV
if any(input_parameters['Species'] == 'Fe1') and any(input_parameters['Elo/eV'] < 2.0):
    # index for these parameters
    indi11 = (input_parameters['Species']  == 'Fe1') & (input_parameters['Elo/eV'] < 2.0) & indi_e
    # calculated aberr
    aberr11 = function_fe1_lt02(input_parameters.drop(columns = ['Species'])[indi11])
    # added to the output df
    outp.loc[indi11,'aberr'] = aberr11

# checks for Fe1 and elo > 2eV
if any(input_parameters['Species'] == 'Fe1') and any(input_parameters['Elo/eV'] > 2.0):
    # index for these parameters
    indi12 = (input_parameters['Species']  == 'Fe1') & (input_parameters['Elo/eV'] > 2.0)  & indi_e
    # calculated aberr
    aberr12 = function_fe1_gt02(input_parameters.drop(columns = ['Species'])[indi12])
    # added to the output df
    outp.loc[indi12, 'aberr'] = aberr12

# checks for Fe2
if any(input_parameters['Species'] == 'Fe2'):
    # index for these parameters
    indi2 = (input_parameters['Species']  == 'Fe2') & indi_e
    # calculated aberr
    aberr2 = function_fe2(input_parameters.drop(columns = ['Species'])[indi2])
    # added to the output df
    outp.loc[indi2, 'aberr'] = aberr2


# output df written to defined output file name and which contains: 
# index, teff/K,lg(g/cms^-2), A(Fe), vmic/kms^-1, Elo/eV, Eup, lggf, Species, aberr, extrapolate, in_grid
outp.to_csv(output_fname, index_label = '# index', float_format='%16.8f')

print('---')
print('aberr has been saved to %s'%(output_fname))
print('---')

#
# if __name__ == '__main__':
#     main()