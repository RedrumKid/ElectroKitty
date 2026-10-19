#include <iostream>
#include <cmath>
#include <vector>
#include <new>
#include <algorithm>
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "electokitty_helper_file.cpp"

using namespace std;
using namespace std::chrono;
namespace py = pybind11;

/*
The C++ implementation of the Python simulator
Binded with Pybind11 to Python. It is compiled seperately to the rest of the library.
The Pybind11 bind is done at the end, the rest of the code is the simulator
*/

class Electrokitty_simulator: MINPAC{
public:

/*
The simulator class and variable declarations. Using the vectro library, which is standard with most compilers.
The vectors are the same as with Python, except mechanism list, this one is chopped up, since C++ does not like arrays of mixed types
*/
        double F; // Faraday constant
        double R; // gas constant
        double PI; // Pi, the constant

        vector<double> t; // time vector from the main Python object, needs constant dt
        vector<double> E_generated; // the potential program, supplied by Python code
        vector<double> current; // the current vector to writte into

        vector<vector<double>> concentration_profile; // the concentration profile vector to copy into, is passed as part of the solution
        vector<vector<double>> surface_profile; // the adsorbed species profile

        vector<double> cell_const; // cell constants temperature, resistance, capacitance, area (T, Ru, Cdl, A)

        vector<double> diffusion_const; // diffusion constants for the dissolved species
        vector<double> isotherm; // isotherm constants for the adsorbed species

        vector<double> spectators; // which species are spectators, NOT IMPLEMENTED
        vector<double> spatial_information; // contains how to set up the concentration profile (first point percentage, number of points, viscosity, rotation rate)

        vector<vector<double>> species_information; // initial condition for all species
        vector<vector<double>> kin; // kinetic constants

        vector<vector<double>> cons; // other contants are stored here

        vector<double> x_dir; // the spatial points locations
        int number_of_diss_spec; // number of dissolved species
        int number_of_surf_conf; // number of adsorbed species
        vector<double> E_Corr; // vector to store the actual potential, modified due to resistance
        //mechanism list
        vector<vector<string>> spec; // lisitng species names
        vector<vector<vector<vector<int>>>> index; // index of which species are in which reaction
        vector<vector<int>> types; // which reaction are irrevesible
        vector<vector<int>> r_ind; // indicies of reactions
        vector<double> num_el; // electrons for each electrochemical reaction

        // not used
        vector<int> tells; // not used for fitting
        int gammaposition; 
        Fit_Params fit_params;

        string kinetic_model = "BV"; // which kinetic model to use

// ######## funkcije
        Electrokitty_simulator(){
                /*Initialization, just to set up the constants*/
                F = 96485.3321;
                R = 8.314;
                PI = 3.141592653589793238462643383279502884197;
        }

        void set_params(
        vector<double> scell_const,
        vector<double> sdiffusion_const,
        vector<double> sisotherm,
        vector<vector<double>> sspectators,
        vector<double> sspatial_information,
        vector<vector<double>> sspecies_information,
        vector<vector<double>> skin,
        vector<vector<string>> sspec,
        vector<vector<vector<vector<int>>>> sindex,
        vector<vector<int>> stypes,
        vector<vector<int>> sr_ind,
        vector<double> snum_el,
        string kin_model){
                /*
                function that takes the constants and the information given by the parser and stores it in the object
                
                this information is given by the Python simulator class
                */
                cell_const = scell_const;
                diffusion_const = sdiffusion_const;
                isotherm = sisotherm;
                for(int i = 0; i < 2; i++){
                        for (int j = 0; j<sspectators[i].size(); j++){
                                spectators.push_back(sspectators[i][j]);
                        }
                }
                spatial_information = sspatial_information;
                species_information = sspecies_information;
                kin = skin;
                spec = sspec;
                index = sindex;
                types = stypes;
                r_ind = sr_ind;
                num_el = snum_el;
                kinetic_model = kin_model;
        }

        void set_sim_prog(vector<double> Time, vector<double> E){
                /*
                function that takes the time and potential arrays and stores it in the object

                Input:
                - vector<double> Time: the time vector containing time points in s
                -  vector<double> E: potential for each time point in V
                */
                t = Time;
                E_generated = E;
        }

        vector<double> simulate(){
                /*
                fuction to call for the simulator to run a single cimulation and output the current for each time point

                Output:
                - vector<double> current
                */
                current.clear(); // clear the memory for each vector, otherwise the simulator will keep appending
                E_Corr.clear();
                surface_profile.clear();
                concentration_profile.clear();

                cons.push_back(cell_const); // store the constants
                cons.push_back(diffusion_const);
                cons.push_back(isotherm);
                cons.push_back(spectators);

                simulator_main_loop(
                        spec, index, types, r_ind, num_el, kin, cons,spatial_information, t, species_information,
                        E_generated, 0); // call the main loop
                cons.clear();
                return current;
        }

        vector<vector<double>> give_surf_profile(){
                /*
                a function so that Python can obtain the adsorbed species profile
                */
                return surface_profile;
        }

        vector<vector<double>> give_concentration_profile(){
                /*
                function so that Python can get the adsorbed species profile
                */
                return concentration_profile;
        }

        vector<double> give_current(){
                /*
                function so Python can get the current
                */
                return current;
        }

        vector<double> give_E_Corr(){
                /*
                function so that Python can access the true potentials
                */
                return E_Corr;
        }

        void simulator_main_loop(
                vector<vector<string>> sspec,
                vector<vector<vector<vector<int>>>> iindex,
                vector<vector<int>> ttypes,
                vector<vector<int>> rr_ind,
                vector<double> nnum_el,
                vector<vector<double>> kinetic_cons,
                vector<vector<double>> constants,
                // odstran če se da - se ne da
                vector<double> spatial_info,
                //
                vector<double> time,
                vector<vector<double>> species_info,
                vector<double> potential_program,
                int eqilibration = 0
                )
                {
                        /*
                        the main loop of the simulator going over all time points. 
                        takes the inputed parameters to setup the simulation by perparing helping objects and performs the needed
                        precalc steps, such as creating the derivative coefficients
                        */
                vector<double> isotherm_cons_burner; // vector to store the isotherm constants after renormaliting them and adding 0s for dissolved species
                vector<vector<double>> ads_cons; // rate constants for chemical reactions at the surface
                vector<vector<double>> bulk_cons; // rate constants for reactions if the diffusion layer
                vector<vector<double>> EC_cons; // rate constants for electrochemical reactions
                vector<vector<double>> bound1; // storing coefficients for the boundary at the electrode
                vector<double> bound2; // concentrations at the other boundary, equal to the intial condition
                vector<double> xp; // storing solution to the DAE set here for use 
                vector<double> delta_e; // stores dE/dt for points above the first time point
                vector<vector<vector<double>>> a; // stores coefficient for calculating the derivatives in the diffusion layer
                double dt, velocity_c, i_kin_prev; // time difference, fluid velocity, previous faradaic current
                int n; // size of the solution array and xp
                int iflag = 1; // check of solver
                int info; // info from solver

                // copying
                t = time;
                spec = sspec; 
                index = iindex;
                types = ttypes;
                r_ind = rr_ind;
                num_el = nnum_el;
                number_of_diss_spec = int(spec[1].size()); // get the number of each type of species
                number_of_surf_conf = int(spec[0].size());

                // copy the constants
                cell_const = constants[0];
                diffusion_const = constants[1];
                isotherm = constants[2];
                isotherm_cons_burner = constants[2];
                spectators = constants[3];

                species_information = species_info; // intial condition

                vector<double> null(number_of_diss_spec); // create 0 values for dissolved species, so that the isotherm correction is 1
                for (int i = 0; i<null.size(); i++){
                        null[i] = 0.;
                }
                
                if (number_of_surf_conf > 0){ // check if there are adsorbed species
                        double max_surf_conc;
                        max_surf_conc = *max_element(species_info[0].begin(), species_info[0].end());
                        for (int i = 0; i<isotherm_cons_burner.size(); i++){
                                isotherm_cons_burner[i]/=max_surf_conc; // divide the isotherm constant by the maximum surface concentration 
                        }
                }
                for(int i = 0; i < number_of_diss_spec; i++){
                        isotherm_cons_burner.push_back(0.); // add 0s to the isothem array to avoid errors
                }
                
                // create rate constant vectors to use in the simulation, this is simply reoragnising things
                ads_cons = create_constant_list(r_ind[0], kinetic_cons);
                bulk_cons = create_constant_list(r_ind[1], kinetic_cons);
                EC_cons = create_constant_list(r_ind[2], kinetic_cons);
                
                ads_cons = get_kinetic_constants(ads_cons, types[0]); // update for irreversibility
                bulk_cons = get_kinetic_constants(bulk_cons, types[1]);
                
                vector<vector<vector<double>>> various_constants = {
                        ads_cons,
                        bulk_cons,
                        EC_cons}; // store them in a central vector
                
                dt = time[1] - time[0]; // calculate the time step difference

                if(number_of_diss_spec > 0){
                        x_dir = space_ranges(time[time.size()-1], *max_element(diffusion_const.begin(), diffusion_const.end()),
                        spatial_info[0], int(spatial_info[1])); // calculate the diffusion layer points
                }else{
                        x_dir = space_ranges(time[time.size()-1], 1.,
                        spatial_info[0], int(spatial_info[1])); // if there are no dissolved species create dummy points
                }

                velocity_c = -0.51/sqrt(spatial_info[2])*pow(2*PI*spatial_info[3], 1.5); // calculate the fluid velocity for an RDE

                a = calc_main_coeficients(x_dir, dt, diffusion_const, int(x_dir.size())-2, velocity_c); // get the derivative coefficients

                bound1 = calc_boundary_condition(x_dir, 0, diffusion_const, 3, velocity_c); // get the coefficients for the boundary at the electrode

                n = number_of_surf_conf+number_of_diss_spec*(int(x_dir.size())) + 2; // the size of the solution array 

                // allocate memory for the solver
                double wa[100000];
                int lw = 100000;
                double x[1000];
                double f[1000];
                
                // create the first guess, by assigning the solution array the initial condition
                for (int i = 0; i<number_of_surf_conf; i++){// first slice
                        x[i] = species_info[0][i];  
                }

                for (int i = 0; i < x_dir.size(); i++){ // looping over the diffusion layer
                        for (int j = 0; j<number_of_diss_spec; j++){
                                x[number_of_surf_conf+number_of_diss_spec*i+j] = species_info[1][j];
                        }
                        
                }

                x[n-2] = potential_program[0]; // add the potential
                x[n-1] = 0.; // add capacitive current

                for (int i = 0; i<number_of_diss_spec; i++){
                        bound2.push_back(species_info[1][i]);// create the boundry away from the electrode
                }

                for (int i = 1; i < time.size(); i++){
                        delta_e.push_back((potential_program[i]-potential_program[i-1])/dt); // calculate dE/dt
                }

                // if equilibration is set first a solution is found for the initial condition
                if (eqilibration == 0){
                        params.set_params(int(spatial_info[1]), dt, number_of_surf_conf, number_of_diss_spec, 
                        bound1, bound2, a, null, various_constants, index, isotherm_cons_burner, 
                        spectators, 1., cell_const, kinetic_model);
                        params.set_ec_params(cell_const[0], num_el, types[2]);
                }else{
                        params.set_params(int(spatial_info[1]), dt, number_of_surf_conf, number_of_diss_spec, 
                        bound1, bound2, a, null, various_constants, index, isotherm_cons_burner, 
                        spectators, 0., cell_const, kinetic_model);
                        params.set_ec_params(cell_const[0], num_el, types[2]);
                }
                
                xp = params.copy_array_to_vec(x, n); // create a copy of the initial condition to use in the solver
                
                vector<double> c_bound; // concentrations at the electrode boundary, includes adsorbed species

                // first point is due to initial condition
                // xp = params.copy_array_to_vec(x, n);
                c_bound = params.aslice(x, 0, number_of_diss_spec + number_of_surf_conf); // get the concentrations at the electrode
                current.push_back(0.); // the first current is always 0 
                E_Corr.push_back(x[n-2]); // the first potential is always the one given by the potential program
                surface_profile.push_back(vslice(xp, 0, number_of_surf_conf)); // store initial contidions
                concentration_profile.push_back(vslice(xp, number_of_surf_conf, n-2));
                i_kin_prev = 0.; // current is currently 0

                
                for (int tt = 1; tt < time.size(); tt++){ // starting from second point
                        x[n-2] = potential_program[tt]; // update potential guess
                        params.update_time_step(potential_program[tt], xp, delta_e[tt-1], i_kin_prev); // store previous solution
                        info = hybrd1(n, x, f, 1e-16, wa, lw, params); // solve the DAE set fot the time point
                        xp = params.copy_array_to_vec(x, n); // copy the solution into a vector
                        c_bound = params.aslice(x, 0, number_of_diss_spec + number_of_surf_conf); // get the concentrations at the boundary

                        i_kin_prev = F*params.A*params.calc_current(2, c_bound, x[n-2]); // calculate the faradaic current
                        current.push_back(i_kin_prev + x[n-1]); // store faradaic + capacitive current
                        E_Corr.push_back(x[n-2]); // store actual potential at time 
                        surface_profile.push_back(vslice(xp, 0, number_of_surf_conf)); // store concentrations
                        concentration_profile.push_back(vslice(xp, number_of_surf_conf, n-2));
                }
        }

private:
        Params params; // create the Params object
        int nx; // number of spatial points
        
//functions to call
        vector<vector<double>> get_kinetic_constants(vector<vector<double>> k_vector, vector<int> kinetic_types){
                /*
                function to update the kinetic constants for irreversibility

                Input:
                - vector<vector<double>> k_vector: contains rate constants for specified reactions
                - vector<int> kinetic_types: contains what type of irrversiblity or reversibility it is

                Output:
                - vector<vector<double>> k_vector: corrected constant vector
                */
                for(int i = 0; i < k_vector.size(); i++){
                        if (kinetic_types[i] == 1){
                                if (k_vector[i].size() == 1){
                                        k_vector[i].push_back(0.);
                                }
                                k_vector[i][1] = k_vector[i][0];
                                k_vector[i][0] = 0.;
                        }
                        else if (kinetic_types[i] == 2){
                                if (k_vector[i].size() == 1){
                                        k_vector[i].push_back(0.);
                                }
                                k_vector[i][0] = k_vector[i][0];
                                k_vector[i][1] = 0.;
                        }
                }
                return k_vector;
        }

        double find_gama(double dx, double xmax, int nx){
                /*
                using the bisection method to find the value of gama that is used to construct the exponentially expanding grid
                
                Input:
                - double dx: the distance to the first point in the diffusion layer
                - double xmax: the maximum distance of the diffusion layer
                - int nx: number of spatial points

                Output:
                - double gama: value of gama to construct the diffusion layer points
                */

                double nnx = static_cast<double>(nx);

                double a = 1.; // intial guesses
                double b = 2.;
                double gama;
                double f;

                for (int it = 0; it<=50; it++){ // bisection method to find gama
                        gama = (a+b)/2.;
                        f = dx*(pow(gama,nnx)-1.)/(gama-1.) - xmax;
                        if (f<=0){
                                a = gama;
                        }else{
                                b = gama;
                        }
                        if (abs(b-a)<=1e-8){
                                break;
                        }
                }
                gama = (a+b)/2.;
                if (gama > 2.){ // gama needs to be between 1 and 2
                        throw("bad gama value");
                }
                return gama;
        }

        vector<vector<double>> fornberg_weights(double z, vector<double> x, int n, int m){
                /*

                From Bengt Fornbergs (1998) SIAM Review paper.

                Input Parameters
                z location where approximations are to be accurate,
                x(0:nd) grid point locations, found in x(0:n)
                n one less than total number of grid points; n must
                not exceed the parameter nd below,
                nd dimension of x- and c-arrays in calling program
                x(0:nd) and c(0:nd,0:m), respectively,
                m highest derivative for which weights are sought,

                Output Parameter
                c(0:nd,0:m) weights at grid locations x(0:n) for derivatives
                of order 0:m, found in c(0:n,0:m)
                dimension x(0:nd),c(0:nd,0:m)

                */
                vector<vector<double>>c(n+1, vector<double>(m+1));
                int mn;
                double c1, c2, c3, c4, c5;
                c1 = 1.;
                c4 = x[0] - z;

                c[0][0] = 1.;

                for (int i = 1; i<n; i++){
                        mn = min(i, m);
                        c2 = 1.;
                        c5 = c4;
                        c4 = x[i] - z;
                        for(int j = 0; j<i; j++){
                        c3 = x[i] - x[j];
                        c2 = c3*c2;
                        if (j==(i-1)){
                                for (int k = mn; k>0; k--){
                                c[i][k] = c1*( k*c[i-1][k-1] - c5*c[i-1][k] )/c2;
                                }
                                c[i][0] = -c1*c5*c[i-1][0]/c2;
                        }

                        for (int k = mn; k>0; k--){
                                c[j][k] = ( c4*c[j][k] - k*c[j][k-1] )/c3;
                        }
                        c[j][0] = c4*c[j][0]/c3;
                        }
                        c1 = c2;
                }
                return c;
        }

        vector<double> space_ranges(double tmax, double D, double fraction, int nx){
                /*
                function to calculate the diffusion layer points, based on the exponentially expanding grid

                Input:
                - double tmax: maximum time
                - double D: diffusion constant, should be the highest possible
                - double fraction: fraction of the maximum, where the first point is located
                - int nx: number of points in the diffusion layer

                Output:
                - vector<double> x: spatial grid for the simulation
                */
                double xmax = 6.*sqrt(tmax*D); // calculate the maximum distance
                double dx = fraction*xmax; // distance to the first point
                double gama = find_gama(dx, xmax, nx); // find gama
                vector<double> x(nx+2); 

                for(int i = 0; i<nx+2; i++){
                        x[i] = dx*(pow(gama, i)-1.)/(gama-1.); // calculate the grid point locations
                }
                return x;
        }

        vector<double> vslice(vector<double> copy_from, int i1, int i2){
                /*
                function to copy a part of a vector into a new vector

                Input:
                - vector<double> copy_from: vector to copy from
                - int i1: index where the copy starts
                - int i2: index to copy to

                Output:
                - vector<double> copy: copy of the slice from the original vector
                */
                vector<double> copy; 
                for(int i = i1; i < i2; i++){
                        copy.push_back(copy_from[i]);
                }
                return copy;
        }

        vector<vector<vector<double>>> calc_main_coeficients(vector<double> x, double dt, vector<double> D, int nx, double B){
                /*
                function to calculate the coefficient for estimating the derivatives for solving the diffusion equation

                Input:
                - vector<double> x: the spatial grid points
                - double dt: time difference
                - vector<double> D: diffusion coefficients for each dissolved species
                - int nx: number of spatial points
                - double B: fluid velocity constant

                Output:
                -vector<vector<vector<double>>> main_array: coeficients for estimating the derivative
                */

                // the drivative is estimated using four points 
                vector<double> a1; // allocate to estimate for each point and species
                vector<double> a2; // second point and so on
                vector<double> a3;
                vector<double> a4;
                vector<vector<double>> weights; // calculated weights from fornberg algorithm
                vector<double> xinbtw;
                double alfa1d, alfa2d, alfa3d, alfa4d, alfa1v, alfa2v, alfa3v, alfa4v;

                vector<vector<vector<double>>> main_array(nx, vector<vector<double>>(D.size(), vector<double>(4))); // the array

                for (int i = 1; i<nx; i++){
                        xinbtw = vslice(x, i-1, i+3); // get the points for each grid point to estimate the derivative
                        weights = fornberg_weights(x[i], xinbtw, 4, 2); // calculate the weights
                        alfa1d = weights[0][2]; // copy them out for the second derivative
                        alfa2d = weights[1][2];
                        alfa3d = weights[2][2];
                        alfa4d = weights[3][2];

                        alfa1v = -(B*pow(x[i],2))*weights[0][1]; // use the wieghts to calculate the first derivative and multiply with fluid velocity
                        alfa2v = -(B*pow(x[i],2))*weights[1][1];
                        alfa3v = -(B*pow(x[i],2))*weights[2][1];
                        alfa4v = -(B*pow(x[i],2))*weights[3][1];

                        for (int j = 0; j < D.size(); j++){
                                main_array[i-1][j][0] = dt*(-alfa1d*D[j] - alfa1v); // store the coefficients into the array
                                main_array[i-1][j][1] = dt*(-alfa2d*D[j] - alfa2v)+1.;
                                main_array[i-1][j][2] = dt*(-alfa3d*D[j] - alfa3v);
                                main_array[i-1][j][3] = dt*(-alfa4d*D[j] - alfa4v);
                        }
                }
                return main_array;
        }

        vector<vector<double>> calc_boundary_condition(vector<double> x, int i, vector<double> D, int nx, double B){
                /*
                function to create the boundary condition coefficients. calculates the first derivative coefficients

                Input:
                - vector<double> x: spatial grid points
                - int i: index for where to calculate the boundary condition
                - vector<double> D: diffusion coefficient for the species
                - int nx: number of spatial points
                - double B: constant to calculate fluid velocity
                */

                vector<vector<double>> bound_array(D.size(), vector<double>(3));

                double alfa1, alfa2, alfa3;
                vector<double> xinbtw;
                vector<vector<double>> weights;

                if (i==0){ // at the electrode
                        xinbtw = vslice(x, 0, 3);
                        weights = fornberg_weights(x[i], xinbtw, 3, 1);
                }else if (i == -1){ // away from the electrode
                        xinbtw = vslice(x, nx-3, nx);
                        weights = fornberg_weights(x[nx-1], xinbtw, 3, 1);
                }else{
                        throw("Boundary Error: boundary indexed incorrectly");
                }

                alfa1 = weights[0][1] - B*pow(x[i],2.);
                alfa2 = weights[1][1] - B*pow(x[i],2.);
                alfa3 = weights[2][1] - B*pow(x[i],2.);

                for (int i = 0; i<D.size(); i++){ // store coefficients
                        bound_array[i][0] = -alfa1*D[i];
                        bound_array[i][1] = -alfa2*D[i];
                        bound_array[i][2] = -alfa3*D[i];
                }
                return bound_array;
        }

        vector<vector<double>> create_constant_list(vector<int> indexs, vector<vector<double>> consts){
                /*
                function uses reaction indexes to pull the constants from the kinetic constants 

                Input:
                - vector<int> indexs: indexes of reactions of a certain type to pull
                - vector<vector<double>> consts: all rate constants

                Output:
                - vector<vector<double>> c: rate constants of a certain type
                */
                vector<vector<double>> c;
                for(int i = 0; i < indexs.size(); i++){
                        c.push_back(consts[indexs[i]]); // get the constants from the list and store them
                }
                return c;
        }

        void fcn(int n, double x[], double f[], int &iflag, Params params){
                /*
                function for the DAE set of equation to be solved at each time step
                contains all the ODEs and PDEs, and two equations to calculate the potential drop and capacitive current

                Input:
                - int n: length of x
                - double x[]: array of solutions
                - double f[]: array of difference for RHS and LHS
                - int &iflag: output of hybrd method
                - Params params: the params object to calculate kinetics
                */
                vector<double> kinetics; // kinetics matrix, RHS of the ode, given as a sum
                vector<double> conc_slice; // concentrations at a point in the diffusion layer
                vector<double> bound_slice; // concentration slice at the electrode
                int spec_num;
                double ga;

                bound_slice = params.aslice(x, 0, params.n_ads + params.n_dis); // get the concentrations at the electrode
                kinetics = params.sum_two_vectors(params.calc_kinetics(0, bound_slice, params.isotherm),
                        params.calc_EC_kinetics(2, bound_slice, x[n-2])); // calculate the kinetics at the boundary/electrode

                for(int i = 0; i < params.n_ads; i++){ // ODEs for the adsorbed species
                        f[i] = (x[i] - params.cp[i])/params.dt*params.eqilib - kinetics[i]*params.spectator[i];
                }
                
                if (params.n_dis > 0){ // PDE setup 
                        
                        for(int i = params.n_ads; i<params.n_ads+params.n_dis; i++){ // calculate the boundary condition
                                conc_slice = params.get_conc_at_pos(x, i-params.n_ads, 0, 3, params.n_dis, params.n_ads);
                                f[i] = params.vector_dot_product(params.bound1[i-params.n_ads], conc_slice) - kinetics[i]*params.spectator[i];
                        }

                        for (int xx = 1; xx < params.nx; xx++){ // calculate over the whole diffusion layer
                                conc_slice = params.aslice(x, params.n_ads + params.n_dis*xx, params.n_ads+params.n_dis + xx*params.n_dis);
                                kinetics = params.calc_kinetics(1, conc_slice, params.null);
                                for (int i = params.n_ads + params.n_dis*xx; i < params.n_ads + params.n_dis + xx*params.n_dis; i++){
                                        spec_num = i - params.n_ads - params.n_dis*xx;
                                        conc_slice = params.get_conc_at_pos(x, spec_num, xx-1, xx+3, params.n_dis, params.n_ads);
                                        f[i] = params.vector_dot_product(params.a[xx-1][spec_num], conc_slice)*params.eqilib - 
                                                params.dt*kinetics[spec_num]*params.spectator[params.n_ads + spec_num] - params.cp[i];
                                }
                        }
                        for (int i = n-2*params.n_dis-2; i<n-params.n_dis-2; i++){ // calculate the other boundary condition
                                f[i] = x[i] - params.bound2[i-n+2*params.n_dis+2];
                        }
                        for (int i = n-params.n_dis-2; i<n-2; i++){
                                f[i] = x[i] - x[i-params.n_dis];
                        }
                }
                ga = params.A*F*params.calc_current(2, bound_slice, x[n-2]); // calculate the current

                f[n-2] = params.pnom - x[n-2] - params.Ru*ga - params.Ru*x[n-1]; // IR drop

                f[n-1] = (1+params.A*params.Ru*params.Cdl/params.dt)*x[n-1] - params.A*params.Cdl*params.delta
                         - params.Ru*params.A*params.Cdl*params.cp[n-1]/params.dt 
                         + params.Ru*params.A*params.Cdl*(ga - params.ga_prev)/params.dt; // the capacitive current, from the differential equation
                // f[n-1] = x[n-1] - params.A*params.Cdl*(x[n-2] - params.cp[n-2])/params.dt;
                
        }
};

PYBIND11_MODULE(cpp_ekitty_simulator, m){ // pybind11 binding so that Python can access C++ functions
    py::class_<Electrokitty_simulator>(m, "cpp_ekitty_simulator")
    .def(py::init()) // not all functions are accessible to Python
    .def("set_parameters", &Electrokitty_simulator::set_params)
    .def("set_simulation_programm", &Electrokitty_simulator::set_sim_prog)
    .def("simulator_main_loop", &Electrokitty_simulator::simulator_main_loop)
    .def("give_current", [](Electrokitty_simulator &self){
        py::array current = py::cast(self.give_current());
        return current;
    })
    .def("give_E_corr", [](Electrokitty_simulator &self){
        py::array E_corr = py::cast(self.give_E_Corr());
        return E_corr;
    })
    .def("give_surf_profile", [](Electrokitty_simulator &self){
        py::array surf_p = py::cast(self.give_surf_profile());
        return surf_p;
    })
    .def("give_concentration_profile", [](Electrokitty_simulator &self){
        py::array c_p = py::cast(self.give_concentration_profile());
        return c_p;
    })
    .def("simulate", [](Electrokitty_simulator &self){
        py::array i_sim = py::cast(self.simulate());
        return i_sim; 
    })
    .def_readwrite("current", &Electrokitty_simulator::current)
    .def_readwrite("t", &Electrokitty_simulator::t)
    .def_readwrite("E_generated", &Electrokitty_simulator::E_generated)
    .def_readwrite("concentration_profile", &Electrokitty_simulator::concentration_profile)
    .def_readwrite("surface_profile", &Electrokitty_simulator::surface_profile)
    .def("__getstate__", [](Electrokitty_simulator &self){})
    .def("__setstate__", [](Electrokitty_simulator &self){})
    ;
}
