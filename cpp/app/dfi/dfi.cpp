/* =====================================================================================

                      The UQ Toolkit (UQTk) version @UQTKVERSION@
                          Copyright (@UQTKYEAR@) NTESS
                        https://www.sandia.gov/UQToolkit/
                        https://github.com/sandialabs/UQTk

     Copyright @UQTKYEAR@ National Technology & Engineering Solutions of Sandia, LLC (NTESS).
     Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government
     retains certain rights in this software.

     This file is part of The UQ Toolkit (UQTk)

     UQTk is open source software: you can redistribute it and/or modify
     it under the terms of BSD 3-Clause License

     UQTk is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     BSD 3 Clause License for more details.

     You should have received a copy of the BSD 3 Clause License
     along with UQTk. If not, see https://choosealicense.com/licenses/bsd-3-clause/.

     Questions? Contact the UQTk Developers at <uqtk-developers@software.sandia.gov>
     Sandia National Laboratories, Livermore, CA, USA
===================================================================================== */

// Include std header files
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <string>
#include <list>

// Include UQTk header files
#include "tools.h"
#include "mcmc.h"
#include "amcmc.h"
#include "arrayio.h"
#include "PCSet.h"

// Include other header files
#include "utils.h"
#include "model.h"
#include "likelihood.h"

// Specify default values
#define PRIOR_FILE "prior.dat"                                     // Default prior file
#define DATA_FILES "data.{d}.dat"                                  // Default regex for data files
#define MINDEX_FILES "mindex.{d}.{n}.dat"                          // Default PCE multi-indices file
#define PCCF_FILES "pccf.{d}.{n}.dat"                              // Default PCE coefficients files
#define PUSHFORWARD_MINDEX_FILES "pushforward.mindex.{d}.{n}.dat"  // Default pushforward PCE multi-indices file
#define PUSHFORWARD_PCCF_FILES "pushforward.pccf.{d}.{n}.dat"      // Default pushforward PCE coefficients files
#define PUSHFORWARD_OUTPUT_FILES "pushforward.{d}.{n}.dat"         // Default pushforward output file
#define PUSHFORWARD_MAP_OUTPUT_FILES "pushforward.map.{d}.{n}.dat" // Default pushforward output file
#define CHAIN_FILE "chain.dat"                                     // Default chain output file
#define SYNTHETIC_DATA_FILES "synthetic.data.{d}.dat"              // Default synthetic data sets input and output file
#define CHAIN_INIT_FILE "chain.init.dat"                           // Default chain initial values file
#define NB_OF_MCMC 100000                                          // Default number of MCMC steps
#define BETA 1                                                     // Default scaling factor for the standard deviation of the synthetic data sets
#define NB_OF_SYNTHETIC_DATA_SETS 100                              // Default number of synthetic data instances per data set
#define OUTPUT_FREQUENCY 0.1                                       // Default output frequency
#define PROPOSAL_JUMP_SIZE 0.01                                    // Default proposal jump size
#define BURNIN 0.1                                                 // Default MCMC chain burnin
#define EVERY 1                                                    // Default MCMC chain subsamling rate

/**
 * Print options to screen.
 */
void help()
{
    std::cout << "Bayesian calibration for data summaries." << std::endl;
    std::cout << std::endl;
    std::cout << "This program performs calibration of a model represented by a set of" << std::endl;
    std::cout << "polynomial surrogate models using Bayesian inference. Files that contain" << std::endl;
    std::cout << "the experimental data, polynomial expansion coefficients and multi-indices" << std::endl;
    std::cout << "and the prior specifications must be provided." << std::endl;
    std::cout << std::endl;
    std::cout << "Usage:" << std::endl;
    std::cout << "dfi.x [-d <data_files>] [-c <pccf_files>] [-i <mindex_files>] [-p <prior_file>] " << std::endl;
    std::cout << "      [-n <n_mcmc>] [-k <n_synth>] [-b <betas>] [-j <weights>] [-r <seed>] " << std::endl;
    std::cout << "      [-g <jump_size>] [-o] [-u <out_freq>] [-m <burnin>] [-e <every>] [-z] [-q]" << std::endl;
    std::cout << "      [-a <chain>] [-t <chain_init>] [-t <chain_init>] [-f <pushforward_out>]" << std::endl;
    std::cout << "      [-x <pushforward_out>] [-v <pushforward_map_out>]" << std::endl;
    std::cout << "      [-w <pushforward_mindex_files>] [-l] [-s <synthetic_data_files>] [-h]" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -d  Experimental data file name(s) (default = " << DATA_FILES << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "      This file should contain two columns, a first column with the measurement" << std::endl;
    std::cout << "      values and a second column with the measurement errors." << std::endl;
    std::cout << "      Multiple numbered data files can be provided by replacing the '{d}' in " << std::endl;
    std::cout << "      the file name by the appropriate data set number. " << std::endl;
    std::cout << std::endl;
    std::cout << "  -c  Polynomial surrogate coefficient file name(s) (default = " << PCCF_FILES << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "      This file should contain a single column with the values for the polynomial" << std::endl;
    std::cout << "      coefficients. Each file should correspond to a particular measurement station. " << std::endl;
    std::cout << "      Multiple files can be provided by replacing the '{d}' and '{n}' in the " << std::endl;
    std::cout << "      given file name by the appropriate data set number and measurement station number." << std::endl;
    std::cout << std::endl;
    std::cout << "  -i  Polynomial surrogate multi-index file name(s) (default = " << MINDEX_FILES << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "      This file should contain the polynomial chaos multi-index set. Every row contains" << std::endl;
    std::cout << "      the order of the basis polynomial in every dimension. Each file should correspond to" << std::endl;
    std::cout << "      a particular measurement station. Multiple files can be provided by replacing the" << std::endl;
    std::cout << "      '{d}' and '{n}' in the given file name by the appropriate data set number and" << std::endl;
    std::cout << "      measurement station number." << std::endl;
    std::cout << std::endl;
    std::cout << "      Note: Currently, only Legendre uniform polynomial chaos expansions are supported." << std::endl;
    std::cout << std::endl;
    std::cout << "  -p  Name of file with prior specifications (default = " << PRIOR_FILE << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "      This file contains one row for each parameter. Only uniform and Gaussian priors" << std::endl;
    std::cout << "      are supported at this point. Specify uniform priors as 'uniform a b' where 'a' is" << std::endl;
    std::cout << "      the lower bound and 'b' is the upper bound for the given paramer. Specify Gaussian" << std::endl;
    std::cout << "      priors as 'gaussian mu sigma' where 'mu' is the mean and 'sigma' is the standard deviation." << std::endl;
    std::cout << std::endl;
    std::cout << "  -n  Number of MCMC iterations (default = " << NB_OF_MCMC << ")" << std::endl;
    std::cout << "  -k  Number of synthetic data sets (single value or comma-separated list, default = " << NB_OF_SYNTHETIC_DATA_SETS << ")" << std::endl;
    std::cout << "  -b  Scaling factor for synthetic data set variance (single value or comma-separated list, default = " << BETA << ")" << std::endl;
    std::cout << "  -j  Weights for each data set when using weighted likelihood (comma-separated list)" << std::endl;
    std::cout << std::endl;
    std::cout << "      Note: The special argument 'data_size' uses weights that scale the likelihood" << std::endl;
    std::cout << "      according to the data set size" << std::endl;
    std::cout << std::endl;
    std::cout << "  -r  Random seed used to generate synthetic data (default = random)" << std::endl;
    std::cout << "  -g  Proposal jump size in MCMC algorithm (default = " << PROPOSAL_JUMP_SIZE << ")" << std::endl;
    std::cout << "  -o  Perform L-BFGS to find good starting values for chain when specified" << std::endl;
    std::cout << "  -u  Output frequency for MCMC chain" << std::endl;
    std::cout << "  -m  Burn-in for computing pushforward posterior (default = " << BURNIN << ")" << std::endl;
    std::cout << "  -e  Subsampling rate for computing pushforward posterior (default = " << EVERY << ")" << std::endl;
    std::cout << "  -z  Compute the standard deviation of the pushforward posterior predictions when specified" << std::endl;
    std::cout << "  -q  Compute the expectation of the Fisher information matrix over the posterior when specified" << std::endl;
    std::cout << std::endl;
    std::cout << "  -a  MCMC chain output file name (default = " << CHAIN_FILE << ")" << std::endl;
    std::cout << "  -t  Name of file with initial values for the chain (default = " << CHAIN_INIT_FILE << ")" << std::endl;
    std::cout << "  -f  Pushforward output file name(s) (default = " << PUSHFORWARD_OUTPUT_FILES << ")" << std::endl;
    std::cout << "  -x  Pushforward MAP output file name(s) (default = " << PUSHFORWARD_MAP_OUTPUT_FILES << ")" << std::endl;
    std::cout << "  -v  Pushforward polynomial surrogate coefficient file name(s) (default = " << PUSHFORWARD_PCCF_FILES << ")" << std::endl;
    std::cout << "  -w  Pushforward polynomial surrogate multi-index file name(s) (default = " << PUSHFORWARD_MINDEX_FILES << ")" << std::endl;
    std::cout << "  -l  Read MCMC chain from file and compute pushforward posteriors when specified" << std::endl;
    std::cout << "  -s  Synthetic data file name(s) (default = " << SYNTHETIC_DATA_FILES << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "  -h  Print this help message" << std::endl;
}

/**
 * Print success message to screen.
 */
int success()
{
    // Print success message
    std::cout << "╭────────────────────────────────────────────────────╮" << std::endl;
    std::cout << "│                       success!                     │" << std::endl;
    std::cout << "╰────────────────────────────────────────────────────╯" << std::endl;

    // Return value
    return 0;
}

/**
 * Print failure message to screen.
 */
int failure()
{
    // Print success message
    std::cout << "╭────────────────────────────────────────────────────╮" << std::endl;
    std::cout << "│                       failure!                     │" << std::endl;
    std::cout << "╰────────────────────────────────────────────────────╯" << std::endl;

    // Return value
    return 1;
}

/**
 * Runs the main program.
 */
int main(int argc, char *argv[])
{
    // Print startup message
    std::cout << "╭────────────────────────────────────────────────────╮" << std::endl;
    std::cout << "│                        dfi                         │" << std::endl;
    std::cout << "╰────────────────────────────────────────────────────╯" << std::endl;

    //
    // Set default values
    //

    std::string prior_file = PRIOR_FILE;                                     // Name of file that contains lower and upper bounds
    std::string data_set_files = DATA_FILES;                                 // Name of file that contains the measurements
    std::string mindex_files = MINDEX_FILES;                                 // Name of file that contains the multi-indices
    std::string pccf_files = PCCF_FILES;                                     // Name of file that contains the PCE coefficients
    std::string betas = std::string();                                       // Synthetic data set variacne scaling factors
    std::string chain_file = CHAIN_FILE;                                     // Name of the MCMC chain output file
    std::string chain_init_file = CHAIN_INIT_FILE;                           // Name of the file that contains the initial values of the chain
    std::string pushforward_mindex_files = PUSHFORWARD_MINDEX_FILES;         // Name of file that contains the pushforward multi-indices
    std::string pushforward_pccf_files = PUSHFORWARD_PCCF_FILES;             // Name of file that contains the pushforward PCE coefficients
    std::string pushforward_output_files = PUSHFORWARD_OUTPUT_FILES;         // Name of the pushforward output file
    std::string pushforward_map_output_files = PUSHFORWARD_MAP_OUTPUT_FILES; // Name of the pushforward map output file
    std::string synthetic_data_files = SYNTHETIC_DATA_FILES;                 // Name of the synthetic data file
    std::string nb_of_synthetic_data_sets = std::string();                   // Number of data instances to use
    std::string weights = std::string();                                     // Likelihood weights for each data set
    int nb_of_mcmc = NB_OF_MCMC;                                             // Number of MCMC samples
    int seed = time(NULL);                                                   // Random seed to inialize MCMC chain
    double burnin = BURNIN;                                                  // Chain burnin used to compute pushforward posterior
    int every = EVERY;                                                       // Chain downsampling ratio used to compute pushforward posterior
    double proposal_jump_size = PROPOSAL_JUMP_SIZE;                          // Proposal jump size in MCMC chain
    double output_frequency = OUTPUT_FREQUENCY;                              // Output frequency
    bool optimize = false;                                                   // When true, run optimizer before running MCMC
    bool read_synthetic_data = false;                                        // When true, read synthetic data sets from file
    bool read_chain_init = false;                                            // When true, read initial values of chain from file
    bool read_chain = false;                                                 // When true, load chain from file and compute pushforward
    bool compute_fisher_information = false;                                 // When true, compute the expectation of the Fisher information over the posterior
    bool compute_pushforward_variance = false;                               // When true, compute the variance of the pushforward posterior
    Info info;                                                               // Info object to store problem-specific parameters

    //
    // Parse options
    //

    // Loop over all input arguments
    while (true)
    {
        switch (getopt(argc, argv, "a:b:c:d:e:f:g:hi:j:k:lm:n:op:qr:s:t:u:v:w:x:z"))
        {
        case 'a':
            chain_file = optarg;
            continue;
        case 'b':
            betas = optarg;
            continue;
        case 'c':
            pccf_files = optarg;
            continue;
        case 'd':
            data_set_files = optarg;
            continue;
        case 'e':
            every = std::atoi(optarg);
            continue;
        case 'f':
            pushforward_output_files = optarg;
            continue;
        case 'g':
            proposal_jump_size = std::atof(optarg);
            continue;
        case 'h':
            help();
            return success();
        case 'i':
            mindex_files = optarg;
            continue;
        case 'j':
            weights = optarg;
            continue;
        case 'k':
            nb_of_synthetic_data_sets = optarg;
            continue;
        case 'l':
            read_chain = true;
            continue;
        case 'm':
            burnin = std::atof(optarg);
            continue;
        case 'n':
            nb_of_mcmc = std::atoi(optarg);
            continue;
        case 'o':
            optimize = true;
            continue;
        case 'p':
            prior_file = optarg;
            continue;
        case 'q':
            compute_fisher_information = true;
            continue;
        case 'r':
            seed = std::atoi(optarg);
            continue;
        case 's':
            synthetic_data_files = optarg;
            read_synthetic_data = true;
            continue;
        case 't':
            chain_init_file = optarg;
            read_chain_init = true;
            continue;
        case 'u':
            output_frequency = std::atof(optarg);
            continue;
        case 'v':
            pushforward_pccf_files = optarg;
            continue;
        case 'w':
            pushforward_mindex_files = optarg;
            continue;
        case 'x':
            pushforward_map_output_files = optarg;
            continue;
        case 'z':
            compute_pushforward_variance = true;
            continue;
        case -1:
            break;
        case '?':
            help();
            return failure();
        }
        break;
    }

    //
    // Check input parameters
    //

    // TODO ?

    //
    // Prior specification
    //

    // Read the file that contains the prior type and prior parameters
    if (not file_exists(prior_file))
    {
        std::cout << "Prior file " << prior_file << " not found!" << std::endl;
        return failure();
    }
    std::ifstream f(prior_file.c_str());
    std::string prior_type;
    double a, b;
    info.nb_of_parameters = 0;
    while (f >> prior_type >> a >> b)
    {
        Array1D<double> prior;
        if (prior_type == "uniform")
        {
            prior.PushBack(0);
        }
        else if (prior_type == "gaussian")
        {
            prior.PushBack(1);
        }
        else
        {
            std::cout << "Unknown prior type " << prior_type << "!";
            return failure();
        }
        prior.PushBack(a); // First prior parameter
        prior.PushBack(b); // Second prior parameter
        info.prior.PushBack(prior);
        info.nb_of_parameters += 1; // Update number of parameters
    }

    // Print all parameters
    std::cout << "Found prior:" << std::endl;
    for (int i = 0; i < info.nb_of_parameters; i++)
    {
        if (info.prior(i)(0) == 0)
        {
            std::cout << "  => p" << std::left << std::setw(4) << i + 1 << ":"
                      << " min = " << std::setw(10) << info.prior(i)(1) << ", max = " << std::setw(10) << info.prior(i)(2) << std::endl;
        }
        if (info.prior(i)(0) == 1)
        {
            std::cout << "  => p" << std::left << std::setw(4) << i + 1 << ":"
                      << " mu = " << std::setw(10) << info.prior(i)(1) << ", sigma = " << std::setw(10) << info.prior(i)(2) << std::endl;
        }
    }

    //
    // Experimental data sets
    //

    // Load experimental data sets
    info.nb_of_data_sets = 0;
    while (true)
    {
        std::string data_set_file = data_set_files;
        size_t pos = data_set_files.find("{d}", 0);
        if (pos != std::string::npos)
            data_set_file.replace(pos, 3, std::to_string(info.nb_of_data_sets + 1));
        if (file_exists(data_set_file))
        {
            Array2D<double> data;
            read_datafileVS(data, data_set_file.c_str());
            info.data_sets.PushBack(data);
            info.nb_of_data_sets += 1;
            info.nb_of_measurement_stations.PushBack(data.XSize());
        }
        else
        {
            break;
        }
        if (pos == std::string::npos)
        {
            break;
        }
    }

    // Check if we have at least 1 data set
    if (info.nb_of_data_sets < 1)
    {
        std::cout << "No data sets found!" << std::endl;
        return failure();
    }

    // Print all data points
    std::cout << "Found " << info.nb_of_data_sets << " data set" << ((info.nb_of_data_sets > 1) ? "s:" : ":") << std::endl;
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        std::cout << "  => Data set number " << d + 1 << "/" << info.nb_of_data_sets << ":" << std::endl;
        for (int i = 0; i < info.nb_of_measurement_stations(d); i++)
        {
            std::cout << "       -> y" << std::left << std::setw(4) << i + 1 << "="
                      << " " << std::setw(10) << info.data_sets(d)(i, 0) << "s" << std::left << std::setw(4) << i + 1 << "="
                      << " " << std::setw(10) << info.data_sets(d)(i, 1) << std::endl;
        }
    }

    //
    // PCE surrogate models
    //

    // Resize pccfs and mindices
    info.pccfs.Resize(info.nb_of_data_sets);
    info.mindices.Resize(info.nb_of_data_sets);

    // Loop over each data set
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        std::cout << "Reading " << info.nb_of_measurement_stations(d) << " PCE models for data set number " << d + 1 << "/" << info.nb_of_data_sets << "..." << std::endl;

        // Resize pccfs and mindices
        info.pccfs(d).Resize(info.nb_of_measurement_stations(d));
        info.mindices(d).Resize(info.nb_of_measurement_stations(d));

        // Loop over each measurement station and load PCE coefficients and multi-indices
        for (int n = 0; n < info.nb_of_measurement_stations(d); n++)
        {

            // Change {d} and {n} in pccf_files by the data set and measurement station number
            std::string pccf_file = pccf_files;
            size_t pos = pccf_file.find("{d}", 0);
            if (pos != std::string::npos)
                pccf_file.replace(pos, 3, std::to_string(d + 1));
            pos = pccf_file.find("{n}", 0);
            if (pos != std::string::npos)
                pccf_file.replace(pos, 3, std::to_string(n + 1));

            // Read and store the PCE coefficients
            if (file_exists(pccf_file))
            {
                read_datafileVS(info.pccfs(d)(n), pccf_file.c_str());
            }
            else
            {
                std::cout << "Could not find file " << pccf_file << std::endl;
                return failure();
            }

            // Change {d} and {n} in mindex_files by the data set and measurement station number
            std::string mindex_file = mindex_files;
            pos = mindex_file.find("{d}", 0);
            if (pos != std::string::npos)
                mindex_file.replace(pos, 3, std::to_string(d + 1));
            pos = mindex_file.find("{n}", 0);
            if (pos != std::string::npos)
                mindex_file.replace(pos, 3, std::to_string(n + 1));

            // Read and store the PCE coefficients
            if (file_exists(mindex_file))
            {
                read_datafileVS(info.mindices(d)(n), mindex_file.c_str());
            }
            else
            {
                std::cout << "Could not find file " << mindex_file << std::endl;
                return failure();
            }

            // Print PCE summary
            std::cout << "  => PCE number " << std::left << n + 1 << ":" << std::endl;
            std::cout << "       -> PCE coefficients read from file " << pccf_file << " (" << info.pccfs(d)(n).XSize() << " coefficients)" << std::endl;
            std::cout << "       -> Multi-indices read from file " << mindex_file << " (" << info.mindices(d)(n).XSize() << " x " << info.mindices(d)(n).YSize() << ")" << std::endl;
        }
    }

    // Consistency check: loop over each data set and print first PCE coefficient and measurement mean
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        std::cout << "Checking consistency of PCE models and data set " << d + 1 << "/" << info.nb_of_data_sets << "..." << std::endl;
        for (int i = 0; i < info.nb_of_measurement_stations(d); i++)
        {
            std::cout << "  => PCE number " << std::left << std::setw(5) << i + 1 << ": first PCE coefficient = " << std::setw(10) << info.pccfs(d)(i)(0) << "; measurement mean = " << std::setw(10) << info.data_sets(d)(i, 0) << std::endl;
        }
    }

    // Compose and save PCEs
    std::cout << "Constructing PCE surrogate models..." << std::endl;
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        Array1D<PCSet *> pces;
        for (int i = 0; i < info.nb_of_measurement_stations(d); i++)
        {
            PCSet *pce = new PCSet("NISPnoq", info.mindices(d)(i), "LU");
            pces.PushBack(pce);
        }
        info.pces.PushBack(pces);
        std::cout << "  => " << info.nb_of_measurement_stations(d) << " PCE surrogate models created for data set " << d + 1 << "/" << info.nb_of_data_sets << std::endl;
    }

    //
    // Synthetic data sets
    //

    // Parse standard deviation scaling factors
    if (betas.empty())
    {
        info.betas.Resize(info.nb_of_data_sets, BETA);
    }
    else
    {
        std::string delimiter = ",";
        size_t pos = 0;
        while (pos != std::string::npos)
        {
            pos = betas.find(delimiter);
            info.betas.PushBack(std::stod(betas.substr(0, pos)));
            betas.erase(0, pos + delimiter.length());
        }
    }

    // Print scaling factors
    std::cout << "Found " << info.betas.XSize() << " variance scaling factor" << ((info.nb_of_data_sets > 1) ? "s:" : ":") << std::endl;
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        std::cout << "  => Data set " << std::left << std::setw(5) << d + 1 << ": variance scaling factor = " << std::setw(10) << info.betas(d) << std::endl;
    }

    // Synthetic data sets
    if (read_synthetic_data) // Read synthetic data sets from file
    {
        // Resize synthetic data sets in info
        info.synthetic_data_sets.Resize(info.nb_of_data_sets);
        info.nb_of_synthetic_data_sets.Resize(info.nb_of_data_sets);

        // Loop over all data sets
        for (int d = 0; d < info.nb_of_data_sets; d++)
        {
            std::cout << "Reading synthetic data set for data set number " << d + 1 << "/" << info.nb_of_data_sets << "..." << std::endl;

            // Change {d} in synthetic_data_files by the data set number
            std::string synthetic_data_file = synthetic_data_files;
            size_t pos = synthetic_data_files.find("{d}", 0);
            if (pos != std::string::npos)
                synthetic_data_file.replace(pos, 3, std::to_string(d + 1));

            // Read and store the data set
            if (file_exists(synthetic_data_file))
            {
                read_datafileVS(info.synthetic_data_sets(d), synthetic_data_file.c_str());
                info.nb_of_synthetic_data_sets(d) = info.synthetic_data_sets(d).YSize();

                // Print synthetic data set info
                std::cout << "  => Synthetic data read from file " << synthetic_data_file << " (size " << info.synthetic_data_sets(d).XSize() << "x" << info.synthetic_data_sets(d).YSize() << ")" << std::endl;
            }
            else
            {
                std::cout << "Could not find file " << synthetic_data_file << std::endl;
                return failure();
            }
        }
    }
    else // Generate synthetic data sets
    {

        // Parse number of synthetic data sets
        if (nb_of_synthetic_data_sets.empty())
        {
            info.nb_of_synthetic_data_sets.Resize(info.nb_of_data_sets, NB_OF_SYNTHETIC_DATA_SETS);
        }
        else
        {
            std::string delimiter = ",";
            size_t pos = 0;
            int k = 0;
            while (pos != std::string::npos)
            {
                pos = nb_of_synthetic_data_sets.find(delimiter);
                k = std::stod(nb_of_synthetic_data_sets.substr(0, pos));
                if (k < 1)
                {
                    std::cout << "Number of synthetic data sets must be positive, got " << k << "!" << std::endl;
                    return failure();
                }
                info.nb_of_synthetic_data_sets.PushBack(k);
                nb_of_synthetic_data_sets.erase(0, pos + delimiter.length());
            }

            // Fill nb_of_synthetic_data_sets if needed
            int a = info.nb_of_synthetic_data_sets.XSize();
            int b = info.nb_of_data_sets;
            if (a == 1 && b > a)
            {
                for (int i = 0; i < b - a; i++)
                {
                    info.nb_of_synthetic_data_sets.PushBack(info.nb_of_synthetic_data_sets(0));
                }
            }
        }

        // Create a random number generator and reset the seed
        std::default_random_engine engine;
        engine.seed(seed);

        info.synthetic_data_sets.Resize(info.nb_of_data_sets);

        // Generate synthetic data
        for (int d = 0; d < info.nb_of_data_sets; d++)
        {
            std::cout << "Generating synthetic data set for data set number " << d + 1 << "/" << info.nb_of_data_sets << "..." << std::endl;
            info.synthetic_data_sets(d).Resize(info.nb_of_measurement_stations(d), info.nb_of_synthetic_data_sets(d));
            for (int n = 0; n < info.nb_of_measurement_stations(d); n++)
            {
                // Create normal distribution
                std::normal_distribution<double> distribution(info.data_sets(d)(n, 0), std::sqrt(info.betas(d)) * info.data_sets(d)(n, 1));

                // Generate the synthetic data
                for (int k = 0; k < info.nb_of_synthetic_data_sets(d); k++)
                {
                    info.synthetic_data_sets(d)(n, k) = distribution(engine);
                }
            }

            // Change {d} in synthetic_data_files by the data set number
            std::string synthetic_data_file = synthetic_data_files;
            size_t pos = synthetic_data_file.find("{d}", 0);
            if (pos != std::string::npos)
                synthetic_data_file.replace(pos, 3, std::to_string(d + 1));

            // Write data to file
            write_datafile(info.synthetic_data_sets(d), synthetic_data_file.c_str());

            // Print synthetic data set info
            std::cout << "  => Generated synthetic data set of size " << info.synthetic_data_sets(d).XSize() << " x " << info.synthetic_data_sets(d).YSize() << std::endl;
            std::cout << "  => Written synthetic data set to file " << synthetic_data_file << std::endl;
        }
    }

    //
    // Likelihood weights
    //

    // Parse likelihood weights
    if (weights.empty())
    {
        info.weights.Resize(info.nb_of_data_sets, 1);
    }
    else
    {
        if (weights == "data_size")
        {
            double denominator = 0;
            for (int d = 0; d < info.nb_of_data_sets; d++)
            {
                denominator += 1.0 / info.nb_of_measurement_stations(d);
            }
            for (int d = 0; d < info.nb_of_data_sets; d++)
            {
                info.weights.PushBack((info.nb_of_data_sets / ((double)info.nb_of_measurement_stations(d))) / denominator);
            }
        }
        else
        {
            std::string delimiter = ",";
            size_t pos = 0;
            while (pos != std::string::npos)
            {
                pos = weights.find(delimiter);
                info.weights.PushBack(std::stod(weights.substr(0, pos)));
                weights.erase(0, pos + delimiter.length());
            }
        }

        std::cout << "Using weighted likelihood with weights:" << std::endl;
        for (int d = 0; d < info.nb_of_data_sets; d++)
        {
            std::cout << "  => Data set number " << std::left << d + 1 << ":" << std::setw(10) << info.weights(d) << std::endl;
        }
    }

    //
    // Pushforward PCE surrogate models
    //

    // Pushforward PCE multi-indices and coefficients
    Array1D<Array1D<Array1D<double>>> pushforward_pccfs;
    Array1D<Array1D<Array2D<int>>> pushforward_mindices;

    // Loop over each data set
    for (int d = 0; d < info.nb_of_data_sets; d++)
    {
        std::cout << "Reading pushforward PCE models for data set number " << d + 1 << "/" << info.nb_of_data_sets << "..." << std::endl;
        Array1D<Array1D<double>> pccfs;
        Array1D<Array2D<int>> mindices;

        // Loop over each measurement station and load PCE coefficients and multi-indices
        int n = 0;
        while (true)
        {
            // Change {d} and {n} in pccf_files by the data set and measurement station number
            std::string pushforward_pccf_file = pushforward_pccf_files;
            size_t pos = pushforward_pccf_file.find("{d}", 0);
            if (pos != std::string::npos)
                pushforward_pccf_file.replace(pos, 3, std::to_string(d + 1));
            pos = pushforward_pccf_file.find("{n}", 0);
            if (pos != std::string::npos)
                pushforward_pccf_file.replace(pos, 3, std::to_string(n + 1));

            // Read and store the PCE coefficients
            if (file_exists(pushforward_pccf_file))
            {
                Array1D<double> pccf;
                read_datafileVS(pccf, pushforward_pccf_file.c_str());
                pccfs.PushBack(pccf);
            }
            else
            {
                break;
            }

            // Change {d} and {n} in mindex_files by the data set and measurement station number
            std::string pushforward_mindex_file = pushforward_mindex_files;
            pos = pushforward_mindex_file.find("{d}", 0);
            if (pos != std::string::npos)
                pushforward_mindex_file.replace(pos, 3, std::to_string(d + 1));
            pos = pushforward_mindex_file.find("{n}", 0);
            if (pos != std::string::npos)
                pushforward_mindex_file.replace(pos, 3, std::to_string(n + 1));

            // Load mindex file
            Array2D<int> mindex;
            read_datafileVS(mindex, pushforward_mindex_file.c_str());
            mindices.PushBack(mindex);

            // Print PCE summary
            std::cout << "  => PCE number " << std::left << n + 1 << ":" << std::endl;
            std::cout << "       -> PCE coefficients read from file " << pushforward_pccf_file << " (" << pccfs(n).XSize() << " coefficients)" << std::endl;
            std::cout << "       -> Multi-indices read from file " << pushforward_mindex_file << " (" << mindices(n).XSize() << " x " << mindices(n).YSize() << ")" << std::endl;

            // Update counter
            n += 1;
        }

        if (pccfs.XSize() == 0)
        {
            std::cout << "  => No pushforward PCE surrogate model found for data set " << d + 1 << std::endl;
        }

        // Add PCE coefficients and multi-indices to list of pushforward PCE models
        pushforward_pccfs.PushBack(pccfs);
        pushforward_mindices.PushBack(mindices);

        // Avoid double counting when no {d} in pushforward
        size_t pos = pushforward_pccf_files.find("{d}", 0);
        if (pos == std::string::npos)
            break;
    }

    //
    // Set MCMC chain burnin and subsampling factor
    //
    info.burnin = (0 < burnin && burnin < 1) ? burnin * nb_of_mcmc : burnin;
    info.every = every;

    //
    // Set online quntities to monitor
    //
    info.compute_fisher_information = compute_fisher_information;
    info.compute_pushforward_variance = compute_pushforward_variance;

    //
    // Run MCMC
    //

    Array2D<double> samplesT;
    if (read_chain)
    {
        std::cout << "Reading chain from file " << chain_file << "..." << std::endl;
        Array2D<double> chain;
        read_datafileVS(chain, chain_file.c_str());
        std::cout << "  => Done" << std::endl;

        // Subsample the chain
        std::cout << "Subsampling MCMC chain..." << std::endl;
        int nrows = (chain.XSize() - burnin)/every + 2;
        int ncols = chain.YSize() - 3;
        samplesT.Resize(nrows, ncols);
        for (int row = 0; row < nrows - 1; row ++) {
            for (int col = 0; col < ncols; col ++) {
                samplesT(row, col) = chain(burnin + row*every, col + 1);
            }
        }

        // Append the MAP
        int max = chain(0, ncols - 1);
        int argmax = 0;
        for (int row = 0; row < nrows; row++) {
            double val = chain(row, ncols - 1);
            if (val > max) {
                max = val;
                argmax = row;
            }
        }
        for (int col = 0; col < ncols; col ++) {
            samplesT(nrows - 1, col) = chain(argmax, col + 1);
        }
        std::cout << "  => Done" << std::endl;
    }
    else
    {

        // Load the MCMC chain
        AMCMC chain(eval_log_posterior, &info);

        // Set the dimension of the chain to the number of parameters
        chain.setChainDim(info.nb_of_parameters);

        // Set the random seed of the chain to the given value
        chain.setSeed(seed);

        // Initialize the covariance of the chain to a default value
        double covariance_value = 0.01;
        Array1D<double> chain_covariance(info.nb_of_parameters, covariance_value);
        chain.initChainPropCovDiag(chain_covariance);

        // Set the proposal jump size
        chain.initAMGamma(proposal_jump_size);

        // Set the chain adaptivity parameters (hardcoded for now)
        int adaptstart = 10000;
        int adaptstep = 10;
        int adaptend = nb_of_mcmc;
        chain.initAdaptSteps(adaptstart, adaptstep, adaptend);

        double eps_cov = 1e-8;
        chain.initEpsCov(eps_cov);

        // Specify output options
        int freq = std::max(1, (int)(nb_of_mcmc * output_frequency));
        chain.setOutputInfo("txt", chain_file, freq, freq);

        // Check chain length
        if (nb_of_mcmc < 1)
        {
            std::cout << "Number of MCMC iterations must be positive, got " << nb_of_mcmc << "!" << std::endl;
            return failure();
        }

        // Print chain details
        std::cout << "Setting up MCMC chain of length " << nb_of_mcmc << "..." << std::endl;
        std::cout << "  => Set random seed of the chain to " << seed << std::endl;
        std::cout << "  => Set chain covariance to " << covariance_value << std::endl;
        std::cout << "  => Set proposal jump size to " << proposal_jump_size << std::endl;
        std::cout << "  => Set adaptivity parameters to:" << std::endl;
        std::cout << "       -> Adapt start = " << adaptstart << std::endl;
        std::cout << "       -> Adapt step = " << adaptstep << std::endl;
        std::cout << "       -> Adapt end = " << adaptend << std::endl;
        std::cout << "       -> Eps cov = " << eps_cov << std::endl;
        std::cout << "  => Set chain print frequency to " << freq << std::endl;
        std::cout << "  => Set burnin to " << info.burnin << std::endl;
        std::cout << "  => Set subsampling rate to " << info.every << std::endl;
        if (compute_pushforward_variance)
        {
            std::cout << "  => Computing variance of pushforward posteriors" << std::endl;
        }
        if (compute_fisher_information)
        {
            std::cout << "  => Computing exected Fisher information matrix in posterior samples" << std::endl;
        }

        // Array to hold initial value of the chain
        Array1D<double> chain_init;
        if (read_chain_init)
        {
            read_datafileVS(chain_init, chain_init_file.c_str());
            if (chain_init.XSize() != info.nb_of_parameters)
            {
                std::cout << "Error while reading initial values for MCMC chain, expected " << info.nb_of_parameters << " got " << chain_init.XSize() << "!";
                return 1;
            }
        }
        else
        {
            for (int i = 0; i < info.nb_of_parameters; i++)
            {
                int prior_type = (int)info.prior(i)(0);
                double a = info.prior(i)(1);
                double b = info.prior(i)(2);
                if (prior_type == 0) // Uniform prior
                {
                    chain_init.PushBack((a + b) / 2);
                }
                else if (prior_type == 1) // Gaussian prior
                {
                    chain_init.PushBack(a);
                }
            }
        }

        // Print chain initial values
        std::cout << "  => Intialized chain to" << std::endl;
        for (int i = 0; i < info.nb_of_parameters; i++)
        {
            std::cout << "       -> p" << std::left << std::setw(4) << i + 1 << ": " << std::setw(10) << chain_init(i) << std::endl;
        }

        // Run optimizer
        if (optimize)
        {
            // Print status
            std::cout << "Status of chain before running optimizer:" << std::endl;
            print_status(chain_init, info);

            // Run optimizer
            std::cout << "Running optimizer..." << std::endl;
            chain.runOptim(chain_init);
        }

        // Reset standard deviation parameters
        if (compute_pushforward_variance)
        {
            for (int d = 0; d < info.nb_of_data_sets; d++)
            {
                Array1D<double> m(info.nb_of_measurement_stations(d), 0);
                info.m.PushBack(m);
                Array1D<double> s(info.nb_of_measurement_stations(d), 0);
                info.s.PushBack(s);
            }
        }

        // Reset Fisher information parameters
        if (compute_fisher_information)
        {
            info.grad_log_likelihood.Resize(info.nb_of_parameters, 0);
        }

        // Reset parameters for online updating of statistics
        info.c1 = 0;
        info.c2 = 0;

        // Print status
        std::cout << "Status of chain before running AMCMC:" << std::endl;
        print_status(chain_init, info);

        // Run the MCMC chain
        std::cout << "Starting AMCMC..." << std::endl;
        chain.runChain(nb_of_mcmc, chain_init);

        // Get MAP parameters and save them into the info object
        chain.getMode(chain_init);

        // Compare model outputs to measurements
        std::cout << "Status of chain after running AMCMC:" << std::endl;
        print_status(chain_init, info);

        // Print the standard deviation of the pushforward posteriors
        if (compute_pushforward_variance)
        {
            std::cout << "Standard deviations of pushforward posteriors:" << std::endl;
            for (int d = 0; d < info.nb_of_data_sets; d++)
            {
                std::cout << "  => Data set number " << d + 1 << "/" << info.nb_of_data_sets << std::endl;
                for (int n = 0; n < info.nb_of_measurement_stations(d); n++)
                {
                    std::cout << "       -> x" << std::left << std::setw(4) << n + 1 << ": standard deviation = " << std::setw(20) << std::sqrt(info.s(d)(n) / (info.c2 - 1)) << std::endl;
                }
            }
        }

        // Print the expected Fisher information matrix over the posterior, if required
        if (compute_fisher_information)
        {
            // WARNING: only multivariate Gaussian priors implemented for now...
            std::cout << "Likelihood-informed subspace:" << std::endl;
            for (int i = 0; i < info.nb_of_parameters; i++)
            {
                std::cout << "  => p" << std::left << std::setw(4) << i + 1 << ": " << std::setw(10) << info.prior(i)(2) * info.prior(i)(2) * info.grad_log_likelihood(i) << std::endl;
            }
        }

        // Get the samples in the chain assuming given burnin and downsampling rate
        Array2D<double> samples;
        chain.getSamples(info.burnin, info.every, samples);
        transpose(samples, samplesT);

        // Append the MAP
        samplesT.insertRow(chain_init, samplesT.XSize());

    }

    //
    // Compute the pushforward
    //

    // Evaluate the pushforward PCEs
    std::cout << "Computing pushforward posteriors:" << std::endl;
    int nb_of_pushforwards = pushforward_mindices.XSize();
    int nb_of_mcmc_iters = samplesT.XSize();
    for (int d = 0; d < nb_of_pushforwards; d++)
    {
        std::cout << "  => Computing pushforward posteriors for data set " << d + 1 << "/" << nb_of_pushforwards << "... " << std::endl;
        int nb_of_pushforward_pces = pushforward_mindices(d).XSize();
        if (nb_of_pushforward_pces > 0)
        {
            Array2D<double> pushforward(nb_of_mcmc_iters - 1, 1);
            Array2D<double> pushforward_map(1, 1);
            Array1D<double> y;
            for (int n = 0; n < nb_of_pushforward_pces; n++)
            {
                std::cout << "     => Measurement station " << n + 1 << "/" << nb_of_pushforward_pces << "... " << std::endl;
                // TODO: replace by call to model?
                PCSet pce("NISPnoq", pushforward_mindices(d)(n), "LU");
                int batchsize = 100000;
                for (int start = 0; start < nb_of_mcmc_iters; start += batchsize)
                {
                    int stop = start + batchsize;
                    if (stop > nb_of_mcmc_iters)
                    {
                        stop = nb_of_mcmc_iters;
                    }
                    int nrows = stop - start;
                    std::cout << "        --> Batching samples " << start + 1 << " to " << start + nrows << "... " << std::endl;
                    Array2D<double> samplesPlaceholder(nrows, info.nb_of_parameters);
                    for (int row = 0; row < nrows; row++)
                    {
                        for (int col = 0; col < info.nb_of_parameters; col++)
                        {
                            samplesPlaceholder(row, col) = samplesT(start + row, col);
                        }
                    }
                    pce.EvalPCAtCustPoints(y, samplesPlaceholder, pushforward_pccfs(d)(n));
                    for (int row = 0; row < nrows; row++)
                    {
                        if (start + row < nb_of_mcmc_iters - 1) {
                            pushforward(start + row, 0) = y(row);
                        }
                        else {
                            pushforward_map(0, 0) = y(row);
                        }
                    }
                }

                // Change {d} and {n} in pushforward_output_files by the data set number and measurement station number
                std::string pushforward_output_file = pushforward_output_files;
                size_t pos = pushforward_output_file.find("{d}", 0);
                if (pos != std::string::npos)
                    pushforward_output_file.replace(pos, 3, std::to_string(d + 1));
                pos = pushforward_output_file.find("{n}", 0);
                if (pos != std::string::npos)
                    pushforward_output_file.replace(pos, 3, std::to_string(n + 1));

                // Write pushforward to file
                write_datafile(pushforward, pushforward_output_file.c_str());
                std::cout << "       -> Pushforward written to file " << pushforward_output_file << std::endl;

                // Change {d} and {n} in pushforward_map_output_files by the data set number and measurement station number
                std::string pushforward_map_output_file = pushforward_map_output_files;
                pos = pushforward_map_output_file.find("{d}", 0);
                if (pos != std::string::npos)
                    pushforward_map_output_file.replace(pos, 3, std::to_string(d + 1));
                pos = pushforward_map_output_file.find("{n}", 0);
                if (pos != std::string::npos)
                    pushforward_map_output_file.replace(pos, 3, std::to_string(n + 1));

                // Write pushforward map to file
                write_datafile(pushforward_map, pushforward_map_output_file.c_str());
                std::cout << "       -> Pushforward map written to file " << pushforward_map_output_file << std::endl;
            }
        }
        else
        {
            std::cout << "       -> No pushforward posterior PCE found" << std::endl;
        }
    }

    // Print success message
    return success();
}
