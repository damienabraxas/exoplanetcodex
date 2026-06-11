#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""NLTE corrections functions

Functions to calculate the NLTE corrections aberr, used in the main_aberr.py program.
For more information about how to prepare input files and run the main program, or
how to use the functions provided in this script see the README files.

Visit Github site:
https://github.com/sliljegren/NLTE_corrections

Reach out directly at:
Sofie Liljegren <sofie.liljegren@astro.su.se>

"""

import pickle
import sklearn

# Load the Fe1 and Fe2 model
[scaler_fe1_lt02, fe1_ml_model_lt02] = pickle.load(open('fe1_model_lt02.p', 'rb'))
[scaler_fe1_gt02, fe1_ml_model_gt02] = pickle.load(open('fe1_model_gt02.p', 'rb'))

[scaler_fe2, fe2_ml_model] = pickle.load(open('fe2_model.p', 'rb'))

def function_fe1_lt02(input_parameters):
    '''
    This function applies a neural network machine learning algorithm, that has been trained on Fe1 data with elo < 0.2, on lines with parameters defined in the input_parameter array to return their NLTE corrections. The metodology is outlined in the README and in the paper, Amarsi+2022.
    
    Parameters
    ----------
    input_parameter:    array
                        Df of parameters of lines - columns are 'teff/K','lg(g/cms^-2)','A(Fe)', 'vmic/kms^-1', 'Elo/eV', 'Eup', 'lggf'
    
    Returns
    ------
    abcorr:             array
                        The NLTE corrections of the lines, defined by the input parameters.
    '''
    # re-scales the data
    input_parameters = scaler_fe1_lt02.transform(input_parameters)
    
    # uses model to predict NLTE corrections
    abcorr = fe1_ml_model_lt02.predict(input_parameters)
    return abcorr

def function_fe1_gt02(input_parameters):
    '''
    This function applies a neural network machine learning algorithm, that has been trained on Fe1 data with elo > 0.2, on lines with parameters defined in the input_parameter array to return their NLTE corrections. The metodology is outlined in the README and in the paper, Amarsi+2022.
    
    Parameters
    ----------
    input_parameter:    pandas df
                        Df of parameters of lines - columns are 'teff/K','lg(g/cms^-2)','A(Fe)', 'vmic/kms^-1', 'Elo/eV', 'Eup', 'lggf'
    
    Returns
    ------
    abcorr:             array
                        The NLTE corrections of the lines, defined by the input parameters.
    '''
    # re-scales the data
    input_parameters = scaler_fe1_gt02.transform(input_parameters)
    
    # uses model to predict NLTE corrections
    abcorr = fe1_ml_model_gt02.predict(input_parameters)
    return abcorr




def function_fe2(input_parameters):
    '''
    This function applies a neural network machine learning algorithm, that has been trained on Fe2 data, on lines with parameters defined in the input_parameter array to return their NLTE corrections. The metodology is outlined in the README and in the paper, Amarsi+2022.
    
    Parameters
    ----------
    input_parameter:    pandas df
                        Df of parameters of lines - columns are 'teff/K','lg(g/cms^-2)','A(Fe)', 'vmic/kms^-1', 'Elo/eV', 'Eup', 'lggf'
    
    Returns
    ------
    abcorr:             array
                        The NLTE corrections of the lines, defined by the input parameters.
    '''
    # re-scales the data
    input_parameters = scaler_fe2.transform(input_parameters)
    
    # uses model to predict NLTE corrections
    abcorr = fe2_ml_model.predict(input_parameters)
    return abcorr