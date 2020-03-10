#!/usr/bin/env python
#=====================================================================================
#
#                      The UQ Toolkit (UQTk) version 3.1.0
#                          Copyright (2020) NTESS
#                        https://www.sandia.gov/UQToolkit/
#                        https://github.com/sandialabs/UQTk
#
#     Copyright 2020 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
#     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
#     retains certain rights in this software.
#
#     This file is part of The UQ Toolkit (UQTk)
#
#     UQTk is open source software: you can redistribute it and/or modify
#     it under the terms of BSD 3-Clause License
#
#     UQTk is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     BSD 3 Clause License for more details.
#
#     You should have received a copy of the BSD 3 Clause License
#     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.
#
#     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
#     Sandia National Laboratories, Livermore, CA, USA
#=====================================================================================
# This file contains many of the tools need for the polynomial tutorial,
# including the class for the model

from __future__ import print_function # To make print() in Python 2 behave like in Python 3

import math
import sys
import getopt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


try:
    import numpy as np
except ImportError:
    print("Numpy module could not be found")


try:

    import PyUQTk.inference as uqtkinf
except ImportError:
    print("PyUQTk inference module not found")

################################################################################
def fit_dram(model):
    """
    Fit the model parameters with Bayesian inference using the DRAM sampler from UQTk

    Input:
        model: Model object, examples can be found below
               All relativent fitting infomation is stored in the model_info attribute of the object
               need:    model_info["method"]:either am, or dram
                        model_info["nsteps"]: number of MCMC steps to take
                        model_info["cvini"]: initial covariance matrix as a diagonal matrix with parameter sampling ranges divided by scale factor as entries.
                        model_info["spllo"]: Lower bound of prior
                        model_info["splhi"]:upper bound of prior
                        model_info["gamma"]: factor to multiply proposed jump size with in the chain past the burn-in phase
                                             (Reduce this factor to get a higher acceptance rate.)
                        model_info["cini"]:initial guess for MCMC chain


    Output:
        A list with samples and other information. See the dram function for more details.
    """

    # define MCMC parameters

    model_info = model.model_info
    method = model_info["method"]
    nsteps = model_info["nsteps"]
    cvini = model_info["cvini"]
    spllo = model_info["spllo"]
    splhi = model_info["splhi"]
    gamma = model_info["gamma"]
    cini = model_info["cini"]

    #make dictonary to pass into UQTk function
    opts = {'method':method,'nsteps':nsteps,'nburn':int(nsteps/10),'nadapt':int(nsteps/1000),'nfinal':10000000,
           'inicov':cvini,'coveps':1.e-10,'burnsc':5,'gamma':gamma,'ndr':2,'drscale':[5,4,3],
           'spllo':spllo,'splhi':splhi}

    sol = uqtkinf.mcmc.dram(opts,cini,model.postAWSM,{})
    return sol

################################################################################
class model(object):
    """Making the basic model class. This one is the general form, will not be used for any particular model
    Currently does not do anything, but could be useful if we add more models"""
    def __init__(self):
        pass

################################################################################
class model_letter(model):
    """This is the basic structure for the letter models (A, B, C, D, ...)
    Contains all the basic functions
    Cannot be used for a model object. The prediction method needs to be specificed in a more specific model

    To declare a class like this, you need an xml file like the example, input.xml
    """

    #Attributes:
    #blacnk dictionay that will be filled in init
    model_info = {}

    def __init__(self, case_info_tree, my_model, input_data, error_level):
        '''
        Sets up the model object

        Gets all model-specific information out of the xml tree and puts it in the
        proper data strutures for later use by the MCMC and postAWSM routines.

        Note: format is such that it can accommmodates a general set of model setting
        parameters that are specified through the xml input file.

        input:
            case_info_tree: the root of the case info, example: ceria-zinkevich
            my_model: which model from the case, example: model_A
            input_data: the raw data, in delta, p, T format
        '''
        super(model_letter, self).__init__()
        print("Making object for ", my_model)

        # Get model information (including information needed to compute log posterior)
        model_info_tree = case_info_tree.find("model_types/" + my_model)
        # Get names and initial values of parameters, as well as the low and high sample limits
        # param_names_dict: dictonary of parameter names, with parameter names as keys and indexes as values
        params_elem = model_info_tree.find("infer_params")
        if params_elem is None:
            print("ERROR: get_model_info: Could not get parameter element from input file")
            sys.exit(1)
        param_names_dict = {}

        #get prior ranges and intital starting point for all parameters
        cini_list = []
        spllo_list = []
        splhi_list = []
        idx = 0
        for param in params_elem.findall("param"):
            param_names_dict[param.get("name")] = idx
            cini_list.append(float(param.get("cini")))
            spllo_list.append(float(param.get("spllo")))
            splhi_list.append(float(param.get("splhi")))
            idx += 1
        if idx == 0:
            print("WARNING: get_model_info: no model parameters to infer found in input file")

        #turn prior ranges into numpy arrays, (makes it easier to do certain functions later)
        cini = np.array(cini_list)
        spllo = np.array(spllo_list)
        splhi = np.array(splhi_list)

        # Form dictionary with all relevant information for the model
        self.model_info={'param_names_dict':param_names_dict, 'cini':cini, 'spllo':spllo, 'splhi':splhi, 'model_name': my_model}

        #get mcmc inflromation
        mcmc_info_tree = model_info_tree.find("mcmc_settings")
        # Number of steps to take
        self.model_info["nsteps"]  = int(mcmc_info_tree.get("n_steps"))

        # Factor to divide prior range by to set initial covariance matrix of proposal distribution
        scale_fac = float(mcmc_info_tree.get("scale_fac"))

        # Compute an initial covariance matrix as a diagonal matrix with parameter sampling ranges
        # divided by scale factor as entries.
        self.model_info["cvini"] = np.diag((splhi-spllo)/scale_fac)
        self.model_info["method"]  = mcmc_info_tree.get("method")

        #factor to multiply proposed jump size with in the chain past the burn-in phase (Reduce this factor to get a higher acceptance rate.)
        self.model_info["gamma"]  = float(mcmc_info_tree.get("gamma"))

        # Assemble dictionary with information for log posterior computation by adding
        # the input data file to the model_info dictionary.
        self.model_info['data'] = input_data

        #get output file name
        self.model_info["output_file"] = model_info_tree.find("mcmc_settings").get("output_file")

        # Post processing settings
        post_settings_tree = model_info_tree.find("post_proc_settings")
        self.model_info["n_skip"] = int(post_settings_tree.get("skip"))   # number of samples to skip at beginning of chain
        self.model_info["stride"] = int(post_settings_tree.get("stride")) # stride for reading in samples

        #number of iterations of MCMC
        self.model_info["number"] = 0

        #save the error error_level
        self.model_info["error_level"] = error_level

    ####################################################
    def postAWSM(self, fit_params, dict):
        """
        Compute the posterior value for the current parameter sample. Since we are working
        with a uniform prior, the posterior probability is directly proportional to the likelihood, so
        we are computing the likelihood here only. For numerical stability, the log posterior is
        computed. The noise model assumes uncorrelated Gaussian noise on the z data only. The standard
        deviation is inferred as a hyper parameter.

        Input:
            fit_params: current sample of the parameter set
            dict: currently unused dictionary, only here because dram calls the
            function with a dictonary of input values, If you need to add additional
            values to compute the posterior, include here
            All other relativent infomation is stored in the model_info attribute of the object
            need:   param_names_dict: dictonary of parameter names and indexes
                    model_name: name of this model
                    fit_params: the parameters that are currently being tested
                    spllo: array of lower end of prior
                    splhi: array of upper end of prior
                    number: number of MCMC steps so far

        Output:
            [ln_likelihood, ln_prior]: a list with the natural log of the
                likelihood and prior  for the current sample
        """

        # Extract data from model_info
        param_names_dict = self.model_info['param_names_dict']
        model_name = self.model_info["model_name"]

        #sigma is the error level specified in the xml file
        sig = self.model_info["error_level"]
        # Noise model parameters
        # sig_index = param_names_dict['sigma']
        # sig = fit_params[sig_index]  # assumed standard deviation of noise model

        #Read in bounds
        spllo = self.model_info["spllo"]
        splhi = self.model_info["splhi"]

        #consider when out of bounds
        #when out of range return ln_like = -max because its the log
        for i in range(len(spllo)):
            if(fit_params[i] < spllo[i]):
                return [-sys.maxsize, -sys.maxsize]
            elif(fit_params[i] > splhi[i]):
                return [-sys.maxsize, -sys.maxsize]


        #Calculate prior

        #find area of prior range
        area = 1;
        for i in range(len(spllo)):
            area = area * (splhi[i] - spllo[i])

        ln_prior = math.log(1/area)

        #compute the difference between observations and forward model using model specific routine
        self.model_info["fit_params"] = fit_params
        y_diff = self.get_y_diff()

        n_data = len(y_diff)

        # Compute log likelihood in three terms
        ln_like  = -np.dot(y_diff,y_diff)/(2.0*sig*sig)
        ln_like += -n_data*math.log(sig)
        ln_like += -0.5*n_data*math.log(2.0*math.pi)

        #Show round number every 100000 samples for checking progress
        self.model_info["number"] += 1
        if (self.model_info["number"] % 100000 ==0):
            print(self.model_info["number"])

        #return a list of log likelihood and log prior
        return [ln_like, ln_prior]

    ####################################################################################
    def get_y_diff(self):
        """
        Compute the difference between the model predictions and the observational data
        for letter models

        Input:
            All needed input will be taken from the model_info dictionary
            need    data = x,y input data
                    param_names_dict = dictionay of parameter names and index values
                    fit_params = current fit parameters to test
        Output:
            The difference [y_observed - y_model] (1D numpy array)
        """
        # Extract data from self.model_info
        data = self.model_info['data']
        param_names_dict = self.model_info['param_names_dict']
        fit_params = self.model_info["fit_params"]


        pred = self.prediction()

        data_pred = self.make_pred_array(pred)

        diff = self.get_diff(data, data_pred)

        return diff


    ################################################################################
    def prediction(self):
        """
        This is a place holder. will be implemented in more specific model
        """
        print("Not implemented here\n Need more specific model")
        return 0

    ################################################################################
    def compute_derivative(self,x_values, params):
        """
        This is a place holder. will be implemented in more specific model
        """
        print("Not implemented here\n Need more specific model")
        return 0

    ################################################################################
    def make_pred_array(self, y_pred):
        """makes the predicted array of [x, y_pred]"""

        data = self.model_info["data"]
        data_pred = np.zeros_like(data)
        data_pred[:,1] = y_pred
        data_pred[:,0] = data[:,0]
        return data_pred

    ##############################################################################
    def get_diff(self,data, data_pred):
        """ Compute difference between the predicted and observed y values"""
        y_diff = data[:,1] - data_pred[:,1]

        return y_diff

####################################################################################
class model_A(model_letter):
    """Model_A, parameters a, b are active here"""
    def __init__(self, case_info_tree, input_data, error_level):
        super(model_A, self).__init__(case_info_tree, "model_A", input_data, error_level)

    ################################################################################
    def prediction(self):
        """
        Compute y for model_A: a + bx

        Input:
            All  inputs are read in from the dictonary model_info
                need:   data: numpy array with observations in x y
                        fit_params: array with current sample of the parameter set
                        param_names_dict: dictionary of parameter indices
        Output:
            y_pred: numpy array with model predictions
        """
        data = self.model_info["data"]
        param_names_dict = self.model_info['param_names_dict']
        fit_params = self.model_info["fit_params"]

        # Extract parameters
        a   = fit_params[param_names_dict["a"]]
        b  = fit_params[param_names_dict["b"]]


        # Compute prediction of z with fitting function using current parameter sample
        y_pred  = a + b * data[:,0]

        return y_pred

    ################################################################################
    def compute_derivative(self,x_values, params):
        """
        Calculates the deivative at the x_values

        inputs:
            x_values: array of x_values to calculate the dervaitve at
            params: current parameters to find derivative at
            Other input read in from the dictonary model_info
            need    param_names_dict: dictionary of parameter indices
        """

        param_names_dict = self.model_info['param_names_dict']

        # Extract parameters
        a   = params[param_names_dict["a"]]
        b  = params[param_names_dict["b"]]

        #this derivative will be constant, just return an array of same value in right size.
        return np.full(len(x_values), b)



####################################################################################
class model_B(model_letter):
    """Model_B, parameters a, b, c are active here"""
    def __init__(self, case_info_tree, input_data, error_level):
        super(model_B, self).__init__(case_info_tree, "model_B", input_data, error_level)

    ################################################################################
    def prediction(self):
        """
        Compute y for model_B: a + bx + cx^2

        Input:
            All  inputs are read in from the dictonary model_info
                need:   data: numpy array with observations in x y
                        fit_params: array with current sample of the parameter set
                        param_names_dict: dictionary of parameter indices
        Output:
            y_pred: numpy array with model predictions
        """
        data = self.model_info["data"]
        param_names_dict = self.model_info['param_names_dict']
        fit_params = self.model_info["fit_params"]

        # Extract parameters
        a   = fit_params[param_names_dict["a"]]
        b  = fit_params[param_names_dict["b"]]
        c  = fit_params[param_names_dict["c"]]


        # Compute prediction of z with fitting function using current parameter sample
        y_pred  = a + b * data[:,0] + c * data[:,0]**2

        return y_pred

    ################################################################################
    def compute_derivative(self,x_values, params):
        """
        Calculates the deivative at the x_values

        inputs:
            x_values: array of x_values to calculate the dervaitve at
            params: current parameters to find derivative at
            Other input read in from the dictonary model_info
            need    param_names_dict: dictionary of parameter indices
        """

        param_names_dict = self.model_info['param_names_dict']

        # Extract parameters
        a   = params[param_names_dict["a"]]
        b  = params[param_names_dict["b"]]
        c  = params[param_names_dict["c"]]

        return np.full(len(x_values), b) + 2 * c * x_values

####################################################################################
class model_C(model_letter):
    """Model_C, parameters a, b, c, d are active here"""
    def __init__(self, case_info_tree, input_data, error_level):
        super(model_C, self).__init__(case_info_tree, "model_C", input_data, error_level)

    ################################################################################
    def prediction(self):
        """
        Compute y for model_C: a + bx + cx^2 + dx^3

        Input:
            All  inputs are read in from the dictonary model_info
                need:   data: numpy array with observations in x y
                        fit_params: array with current sample of the parameter set
                        param_names_dict: dictionary of parameter indices
        Output:
            y_pred: numpy array with model predictions
        """
        data = self.model_info["data"]
        param_names_dict = self.model_info['param_names_dict']
        fit_params = self.model_info["fit_params"]

        # Extract parameters
        a   = fit_params[param_names_dict["a"]]
        b  = fit_params[param_names_dict["b"]]
        c  = fit_params[param_names_dict["c"]]
        d  = fit_params[param_names_dict["d"]]


        # Compute prediction of z with fitting function using current parameter sample
        y_pred  = a + b * data[:,0] + c * data[:,0]**2  + d * data[:,0]**3

        return y_pred

    ################################################################################
    def compute_derivative(self,x_values, params):
        """
        Calculates the deivative at the x_values

        inputs:
            x_values: array of x_values to calculate the dervaitve at
            params: current parameters to find derivative at
            Other input read in from the dictonary model_info
            need    param_names_dict: dictionary of parameter indices
        """
        param_names_dict = self.model_info['param_names_dict']

        # Extract parameters
        a   = params[param_names_dict["a"]]
        b  = params[param_names_dict["b"]]
        c  = params[param_names_dict["c"]]
        d  = params[param_names_dict["d"]]

        return np.full(len(x_values), b) + 2 * c * x_values + 3 * d * x_values ** 2
####################################################################################
class model_D(model_letter):
    """Model_D, parameters a, b, c, d, e are active here"""
    def __init__(self, case_info_tree, input_data, error_level):
        super(model_D, self).__init__(case_info_tree, "model_D", input_data, error_level)

    ################################################################################
    def prediction(self):
        """
        Compute y for model_D: a + bx + cx^2 + dx^3 + ex^4

        Input:
            All  inputs are read in from the dictonary model_info
                need:   data: numpy array with observations in x y
                        fit_params: array with current sample of the parameter set
                        param_names_dict: dictionary of parameter indices
        Output:
            y_pred: numpy array with model predictions
        """
        data = self.model_info["data"]
        param_names_dict = self.model_info['param_names_dict']
        fit_params = self.model_info["fit_params"]

        # Extract parameters
        a   = fit_params[param_names_dict["a"]]
        b  = fit_params[param_names_dict["b"]]
        c  = fit_params[param_names_dict["c"]]
        d  = fit_params[param_names_dict["d"]]
        e  = fit_params[param_names_dict["e"]]


        # Compute prediction of z with fitting function using current parameter sample
        y_pred  = a + b * data[:,0] + c * data[:,0]**2  + d * data[:,0]**3 + e * data[:,0]**4

        return y_pred

    ################################################################################
    def compute_derivative(self,x_values, params):
        """
        Calculates the deivative at the x_values

        inputs:
            x_values: array of x_values to calculate the dervaitve at
            params: current parameters to find derivative at
            Other input read in from the dictonary model_info
            need    param_names_dict: dictionary of parameter indices
        """
        param_names_dict = self.model_info['param_names_dict']

        # Extract parameters
        a   = params[param_names_dict["a"]]
        b  = params[param_names_dict["b"]]
        c  = params[param_names_dict["c"]]
        d  = params[param_names_dict["d"]]
        e  = params[param_names_dict["e"]]

        return np.full(len(x_values), b) + 2 * c * x_values + 3 * d * x_values ** 2 + 4 * e * x_values ** 3

####################################################################################
class model_E(model_letter):
    """Model_E, parameters a, b, c, d, e, f are active here"""
    def __init__(self, case_info_tree, input_data, error_level):
        super(model_E, self).__init__(case_info_tree, "model_E", input_data, error_level)

    ################################################################################
    def prediction(self):
        """
        Compute y for model_E: a + bx + cx^2 + dx^3 + ex^4 + fx^5

        Input:
            All  inputs are read in from the dictonary model_info
                need:   data: numpy array with observations in x y
                        fit_params: array with current sample of the parameter set
                        param_names_dict: dictionary of parameter indices
        Output:
            y_pred: numpy array with model predictions
        """
        data = self.model_info["data"]
        param_names_dict = self.model_info['param_names_dict']
        fit_params = self.model_info["fit_params"]

        # Extract parameters
        a   = fit_params[param_names_dict["a"]]
        b  = fit_params[param_names_dict["b"]]
        c  = fit_params[param_names_dict["c"]]
        d  = fit_params[param_names_dict["d"]]
        e  = fit_params[param_names_dict["e"]]
        f  = fit_params[param_names_dict["f"]]


        # Compute prediction of z with fitting function using current parameter sample
        y_pred  = a + b * data[:,0] + c * data[:,0]**2  + d * data[:,0]**3 + \
                e * data[:,0]**4 + f * data[:,0]**5

        return y_pred

    ################################################################################
    def compute_derivative(self,x_values, params):
        """
        Calculates the deivative at the x_values

        inputs:
            x_values: array of x_values to calculate the dervaitve at
            params: current parameters to find derivative at
            Other input read in from the dictonary model_info
            need    param_names_dict: dictionary of parameter indices
        """
        param_names_dict = self.model_info['param_names_dict']

        # Extract parameters
        a   = params[param_names_dict["a"]]
        b  = params[param_names_dict["b"]]
        c  = params[param_names_dict["c"]]
        d  = params[param_names_dict["d"]]
        e  = params[param_names_dict["e"]]
        f  = params[param_names_dict["f"]]

        return np.full(len(x_values), b) + 2 * c * x_values + 3 * d * x_values ** 2 + 4 * e * x_values ** 3 + 5 * f * x_values ** 4

####################################################################################
###Try to find a better way to implement this part
def make_model_object(my_model, case_info_tree, input_data, error_level):
    """ creates the object for the correct model
    inputs:
        my_model: name of the model to create object for
        case_info_tree: the starting part for the xml tree
        input_data: the input_data that was made in get_data.py
    output:
        a model object of teh correct type
    """
    if my_model == "model_A":
        model = model_A(case_info_tree, input_data, error_level)
    elif my_model == "model_B":
        model = model_B(case_info_tree, input_data, error_level)
    elif my_model == "model_C":
        model = model_C(case_info_tree, input_data, error_level)
    elif my_model == "model_D":
        model = model_D(case_info_tree, input_data, error_level)
    elif my_model == "model_E":
        model = model_E(case_info_tree, input_data, error_level)
    elif my_model == "model_F":
        model = model_F(case_info_tree, input_data, error_level)
    else:
        print(my_model," has not been implemented yet.")
        sys.exit(0)

    return model

################################################################################

def meta_info_correct_posterior(meta_info):
    """
    change the meta_info list from [acceptance, log likelihood, log prior]
    to [acceptance_prob, log posterior]
    This is to pass into writeSamples, so it can pass into inference/postproc.py nicely
    """

    nsteps = len(meta_info)
    new_meta_info = np.zeros((nsteps,2)) # columns for acceptance prob and log posterior

    for i in range(nsteps):
        new_meta_info[i] = [meta_info[i][0],meta_info[i][1] + meta_info[i][2]]

    return new_meta_info

################################################################################
def evaluate_function(coeff, x_values):
    """
    Evaulates a polynoial with the given coefficents at the given x_values

    inputs:
        coeff: list of coefficents, in order of increasing power, i.e. coeff[0] for x^0, coeff[1] for x^1, ...
        x_values: values to evaluate function at

    outputs:
        list of y_values, same size on x_values
    """
    y_values = np.zeros(np.shape(x_values))
    for i in range(len(x_values)):
        x = x_values[i]
        for j in range(len(coeff)):
            y_values[i] += coeff[j] * x ** j
    return y_values

################################################################################
def evaluate_derivative(coeff, x_values):
    """
    Evaulates the derivative of a polynoial with the given coefficents at the given x_values

    inputs:
        coeff: list of coefficents, in order of increasing power, i.e. coeff[0] for x^0, coeff[1] for x^1, ...
        x_values: values to evaluate function at

    outputs:
        list of deriavtives, same size on x_values
    """
    y_values = np.zeros(np.shape(x_values))
    for i in range(len(x_values)):
        x = x_values[i]
        for j in range(1, len(coeff)):
            y_values[i] += j * coeff[j] * x ** (j-1)
    return y_values


################################################################################
def writeSamples(file_name,samples,map_state,meta_info,var_names):
    """
    Write out an ASCII file with the MCMC samples, in the format that
    can be read by the UQTk postproc.py postprocessing script

    Input:
        file_name : name of the output file_name
        samples   : sample values to write out (nsteps x number of chain variables)
        map_state : [map parameter values, posterior probability of the MAP values]
        meta_info : [acceptance probabilities, log posterior] for each sample
        var_names : List of variable names in chain
    Output:
        ASCII file with the MCMC samples, in the format that
        can be read by the UQTk postproc.py postprocessing script
    """
    # Write out file line by line.
    # First column is MCMC step number
    # Next to last column is acceptance probability
    # Last column is the log posterior probability

    # Open the file
    fobj = open(file_name,'w')

    ######Can probably clean this up a little bit

    # First line contains the labels of each variable.
    fobj.write('#steps ')
    for label in var_names:
        fobj.write(label)
        fobj.write(' ')
    fobj.write("acc_prob ")
    fobj.write("log_post ")
    fobj.write("\n")

    # Print sample values
    for row in range(samples.shape[0]):
        fobj.write('{:5d}'.format(row))
        fobj.write(' ')
        for value in samples[row,:]:
            fobj.write('{: .15e}'.format(value))
            fobj.write(' ')
        for value in meta_info[row,:]:
            fobj.write('{: .15e}'.format(value))
            fobj.write(' ')
        fobj.write("\n") #Finish the line

    # Print MAP values on last line, with a posterior value of -1.0, to flag the Python
    # postprocessing script that this is the MAP line
    row += 1
    fobj.write('{:5d}'.format(row))
    fobj.write(' ')
    for value in map_state[0]:
        fobj.write('{: .15e}'.format(value))
        fobj.write(' ')
    fobj.write("-1.0 ")
    fobj.write('{: .15e}'.format(map_state[1]))
    fobj.write(' ')

    fobj.close()

    return
