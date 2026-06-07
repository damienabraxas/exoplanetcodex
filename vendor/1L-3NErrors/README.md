# 1L-3NErrors

Using machine learning neural network trained on Fe data, calculate Fe1 and Fe2
1D LTE - 3D non-LTE abundance errors - for more information see [Amarsi+2022](https://arxiv.org/abs/2209.13449).


 ---

1. Dependencies
2. Running the script
3. Using only `function_aberr.py`

 ---

## 1. Dependencies

To run the script we suggest setting up a virtual environment with python and install the version of packages used to develop the machine learning models. For more information about virtual environments see [the documentation](https://docs.python.org/3/tutorial/venv.html). 

This script uses the following dependencies:

    pandas==1.1.3
    numpy==1.19.3
    scikit-learn==1.0.2

Before using these need to be installed. 

## 2. Running the script

### 2.1 Prepare the input file

The input file contains information about the lines for which the 1D LTE - 3D
non-LTE abundance errors `aberr` are calculated and has the following syntax

    #teff/K, lg(g/cms^-2), A(Fe;3N), vmic/kms^-1,   Elo/eV,   Eup/eV,    lggf, Spec.
    5051,            4.0,      4.5,         0.0, 2.998045, 5.586893, -2.563,   Fe1
    ...

Each row in the input file is the stellar parameter and then the line parameters
for a specific spectral line.

Please note that the abundance errors are functions of the 3D non-LTE iron
abundance, A(Fe;3N).  To correct 1D LTE abundances the user should adopt an
initial guess for the 3D non-LTE abundance and iterate.

It is possible to provide both Fe1, Fe2, and other lines in the same file if the
species are correctly defined. The output species other than Fe1 or Fe2 will be
-999. 

For an example input file, see `test_data.csv`. 

### 2.2 Running the script

This is a `python3` script, that takes (obligatorily) 2 inputs, the
`input_file.csv` input filename, and `output_file.csv` output filename as

    python3 main_aberr.py input_file.csv ouput_file.csv


### 2.3 Extrapolation

The default setting of the code is to not allow extrapolation, i.e. if the given
stellar parameters are outside the calculated model grid, the code will not
calculate a aberr value for that input. It is possible to turn it on by adding a
third input as 

    python3 main_aberr.py input_file.csv ouput_file.csv True

The validity of such calculations has not been tested by the authors, and we
cannot supply errors as such. 

**Using extrapolation is therefore at the discretion of the user, and the
authors cannot guarantee that the results are reasonable.**


### 2.4 The output

The output will have a similar syntax to the input file, with a few added output
columns:

    #index, teff/K, lg(g/cms^-2), A(Fe;3N), vmic/kms^-1,   Elo/eV,   Eup/eV,   lggf, Spec., aberr, extrapolate, in_grid
     0, 5051,          4.0,      4.5,	        0.0,      3.0, 5.586893, -2.563,   Fe1, -0.136,False,       True
	...

- The first column is an index column - `int`
- The aberr column contains the calculated NLTE corrections - `float`
- If the input values are outside the stellar parameter grid, of the species
    is wrong, this value is a default -999 value
- The extrapolation column states if extrapolation is allowed - `bolean`
- The in_grid column states if the stellar parameters of the line are within the
  model grid values - `bolean`


## 3. Using only `function_aberr.py`

It is possible to use the functions provided in `function_aberr.py` directly in
a python code as well by adding 

    from function_aberr import *

to the top of your python script.

`function_aberr.py` contains 3 different functions that calculate the aberr:

    aberr = function_fe1_lt02(input_parameters) #-> calculates aberr for Fe1 with elo < 2eV
	
	aberr = function_fe1_gt02(input_parameters) #-> calculates aberr for Fe1 with elo > 2eV
	
	aberr = function_fe2(input_parameters) 		#-> calculates aberr for Fe2

The output is an numpy array with the aberr results. 

Two caveats
- User has to make sure to use the correct function for Fe1 and Fe2 lines. 
- User has to make sure to use the correct function for Fe1 less than 2 Elo/eV or Fe1 greater than 2 Elo/eV. 

The input is pandas df with the following columns

	    teff/K  lg(g/cms^-2)  A(Fe)  vmic/kms^-1    Elo/eV       Eup     lggf
	0   5051             4.0    4.5          0.0  3.236720  5.825065 -1.76299
	1   5051             4.0    4.5          0.0  3.641640  6.222374 -1.51400
	...

Notice that the functions can also take numpy arrays, as

	array([[5051, 4.0, 4.5, 0.0, 3.23672, 5.825065, -1.76299],
	       [5051, 4.0, 4.5, 0.0, 3.64164, 6.222373999999999, -1.514],
		   ...])
    

which will yield the same aberr output as the data frame, however, with a
warning.
