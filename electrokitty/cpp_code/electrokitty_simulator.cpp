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
pomembni zapiski za c++ simulator:

-nujno mora biti spectator pravilono nastavljen
-podatki iz parserja delajo vredu
-natančnost je primerljiva z Py verzijo
-tol je 1e-8 kar je glede na testiranja dost
-popravljen, da za enosmerne reakcije vzame 1 konst.
*/

class Electrokitty_simulator: MINPAC{
public:
        double F;
        double R;
        double PI;

        vector<double> t;
        vector<double> E_generated;
        vector<double> current;

        vector<vector<double>> concentration_profile;
        vector<vector<double>> surface_profile;

        vector<double> cell_const;

        vector<double> diffusion_const;
        vector<double> isotherm;

        vector<double> spectators;
        vector<double> spatial_information;

        vector<vector<double>> species_information;
        vector<vector<double>> kin;

        vector<double> x_dir;
        int number_of_diss_spec;
        int number_of_surf_conf;
        vector<double> E_Corr;
        //mechanism list
        vector<vector<string>> spec;
        vector<vector<vector<vector<int>>>> index;
        vector<vector<int>> types;
        vector<vector<int>> r_ind;
        vector<double> num_el;

        vector<int> tells;
        int gammaposition;
        Fit_Params fit_params;

// ######## funkcije
        Electrokitty_simulator(){
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
        vector<double> snum_el){
                cell_const = scell_const;
                diffusion_const = sdiffusion_const;
                isotherm = sisotherm;
                for(int i = 0; i<2; i++){
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
        }

        void set_sim_prog(vector<double> Time, vector<double> E){
                t = Time;
                E_generated = E;
        }

        vector<double> simulate(){
                current.clear();
                E_Corr.clear();
                surface_profile.clear();
                concentration_profile.clear();

                vector<vector<double>> cons;
                cons.push_back(cell_const);
                cons.push_back(diffusion_const);
                cons.push_back(isotherm);
                cons.push_back(spectators);

                simulator_main_loop(
                        spec, index, types, r_ind, num_el, kin, cons,spatial_information, t, species_information,
                        E_generated, 0);
                return current;
        }

        void create_optimization_problem(vector<int> ttells, int gampos){
                tells = ttells;
                gammaposition = gampos;
        }

        vector<double> calc_from_guess(vector<double> guess){
                current.clear();
                E_Corr.clear();
                surface_profile.clear();
                concentration_profile.clear();
                Fit_Params fparams = unpack_fit_params(guess, tells, gammaposition);
                vector<vector<double>> cons;
                cons.push_back(fparams.cell_params);
                cons.push_back(diffusion_const);
                cons.push_back(fparams.isotherm);
                cons.push_back(spectators);

                simulator_main_loop(
                        spec, index, types, r_ind, num_el, fparams.kinetics, cons,spatial_information, t, fparams.spec_info,
                        E_generated, 0);
                return current;
        }

        vector<vector<double>> give_surf_profile(){
                return surface_profile;
        }

        vector<vector<double>> give_concentration_profile(){
                return concentration_profile;
        }

        vector<double> give_current(){
                return current;
        }

        vector<double> give_E_Corr(){
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
                vector<double> isotherm_cons_burner;
                vector<vector<double>> ads_cons;
                vector<vector<double>> bulk_cons;
                vector<vector<double>> EC_cons;
                vector<vector<double>> bound1;
                vector<double> bound2;
                vector<double> xp;
                vector<double> delta_e;
                vector<vector<vector<double>>> a;
                double dt, velocity_c;
                int n;
                int iflag = 1;
                int info;
                t = time;
                spec = sspec;
                index = iindex;
                types = ttypes;
                r_ind = rr_ind;
                num_el = nnum_el;
                number_of_diss_spec = int(spec[1].size());
                number_of_surf_conf = int(spec[0].size());

                //spectators more bit nujno 1d vector!!!!!!!!!!!
                
                cell_const = constants[0];
                diffusion_const = constants[1];
                isotherm = constants[2];
                isotherm_cons_burner = constants[2];
                spectators = constants[3];

                species_information = species_info;

                vector<double> null(number_of_diss_spec);
                for (int i = 0; i<null.size(); i++){
                        null[i] = 0.;
                }
                
                if (number_of_surf_conf>0){
                        double max_surf_conc;
                        max_surf_conc = *max_element(species_info[0].begin(), species_info[0].end());
                        for (int i = 0; i<isotherm_cons_burner.size(); i++){
                                isotherm_cons_burner[i]/=max_surf_conc;
                        }
                }
                for(int i = 0; i<number_of_diss_spec; i++){
                        isotherm_cons_burner.push_back(0.);
                }
                
                ads_cons = create_constant_list(r_ind[0], kinetic_cons);
                bulk_cons = create_constant_list(r_ind[1], kinetic_cons);
                EC_cons = create_constant_list(r_ind[2], kinetic_cons);
                
                ads_cons = get_kinetic_constants(ads_cons, types[0]);
                bulk_cons = get_kinetic_constants(bulk_cons, types[1]);
                
                vector<vector<vector<double>>> various_constants = {
                        ads_cons,
                        bulk_cons,
                        EC_cons};
                
                dt = time[1] - time[0];

                if(number_of_diss_spec>0){
                        x_dir = space_ranges(time[time.size()-1], *max_element(diffusion_const.begin(), diffusion_const.end()),
                        spatial_info[0], int(spatial_info[1]));
                }else{
                        x_dir = space_ranges(time[time.size()-1], 1.,
                        spatial_info[0], int(spatial_info[1]));
                }

                velocity_c = -0.51/sqrt(spatial_info[2])*pow(2*PI*spatial_info[3],1.5);

                a = calc_main_coeficients(x_dir, dt, diffusion_const, int(x_dir.size())-2, velocity_c);

                bound1 = calc_boundary_condition(x_dir, 0, diffusion_const, 3, velocity_c);

                n = number_of_surf_conf+number_of_diss_spec*(int(x_dir.size()))+2;
                /* double wa[( n * ( 3 * n + 13 ) ) / 2 + 100];
                int lw = ( n * ( 3 * n + 13 ) ) / 2 + 100;
                double x[number_of_surf_conf+number_of_diss_spec*(x_dir.size())+2];
                double f[number_of_surf_conf+number_of_diss_spec*(x_dir.size())+2]; */

                double wa[100000];
                int lw = 100000;
                double x[1000];
                double f[1000];
                
                for (int i = 0; i<number_of_surf_conf; i++){
                        x[i] = species_info[0][i];
                }

                for (int i = 0; i < x_dir.size(); i++){
                        for (int j = 0; j<number_of_diss_spec; j++){
                                x[number_of_surf_conf+number_of_diss_spec*i+j] = species_info[1][j];
                        }
                        
                }

                x[n-2] = potential_program[0];
                x[n-1] = 0.;
                for (int i = 0; i<number_of_diss_spec; i++){
                        bound2.push_back(species_info[1][i]);
                }

                for (int i = 1; i<time.size(); i++){
                        delta_e.push_back((potential_program[i]-potential_program[i-1])/dt);
                }

                if (eqilibration == 0){
                        params.set_params(int(spatial_info[1]), dt, number_of_surf_conf, number_of_diss_spec, 
                        bound1, bound2, a, null, various_constants, index, isotherm_cons_burner, spectators, 1., cell_const);
                        params.set_ec_params(cell_const[0], num_el, types[2]);
                }else{
                        params.set_params(int(spatial_info[1]), dt, number_of_surf_conf, number_of_diss_spec, 
                        bound1, bound2, a, null, various_constants, index, isotherm_cons_burner, spectators, 0., cell_const);
                        params.set_ec_params(cell_const[0], num_el, types[2]);
                }
                
                xp = params.copy_array_to_vec(x, n);

                vector<double> c_bound;
                for (int tt = 0; tt<time.size(); tt++){
                        x[n-2] = potential_program[tt];
                        params.update_time_step(potential_program[tt], xp, delta_e[tt-1]);
                        info = hybrd1(n, x, f, 1e-15, wa, lw, params);
                        xp = params.copy_array_to_vec(x, n);
                        c_bound = params.aslice(x, 0, number_of_diss_spec+number_of_surf_conf);
                        current.push_back(F*params.A*params.calc_current(2,c_bound,x[n-2]) + x[n-1]);
                        E_Corr.push_back(x[n-2]);
                        surface_profile.push_back(vslice(xp, 0, number_of_surf_conf));
                        concentration_profile.push_back(vslice(xp, number_of_surf_conf, n-2));
                }
        }

private:
        Params params;
        int nx;

        Fit_Params unpack_fit_params(vector<double> gues, vector<int> tell, int gama_pos){
                Fit_Params fparams;
                vector<vector<double>> kinetics;
                vector<double> cell_params;
                vector<vector<double>> spec_info;
                vector<double> isot;
                int index1 = 0;
                int index2;

                cell_params.push_back(cell_const[0]);
                spec_info = species_information;

                for (int i = 0; i<tell[0]; i++){
                        index2 = tell[i+1];
                        kinetics.push_back(vslice(gues, index1, index2));
                        index1 = index2;
                }

                if (tell[tell[0]+1] != 0){
                        cell_params.push_back(gues[tell[tell[0]+1]]); //Ru
                }else{
                        cell_params.push_back(cell_const[1]);
                }

                if (tell[tell[0]+2] != 0){
                        cell_params.push_back(gues[tell[tell[0]+2]]); //Cdl
                }else{
                        cell_params.push_back(cell_const[2]);
                }

                if (tell[tell[0]+3] != 0){
                        cell_params.push_back(gues[tell[tell[0]+3]]); //A
                }else{
                        cell_params.push_back(cell_const[3]);
                }

                if (tell[tell[0]+4] != 0){
                        spec_info[0][gama_pos] = gues[tell[tell[0]+4]];
                }

                if (tell[tell[0]+5] != 0){
                        isot = vslice(gues, tell[tell[0]+5], int(gues.size()));
                }else{
                        isot = isotherm;
                }

                fparams.insert_params(kinetics, cell_params, spec_info, isot);
                return fparams;
        }
        
//functions to call
        vector<vector<double>> get_kinetic_constants(vector<vector<double>> k_vector, vector<int> kinetic_types){
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
                double nnx = static_cast<double>(nx);

                double a = 1.;
                double b = 2.;
                double gama;
                double f;

                for (int it = 0; it<=50; it++){
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
                if (gama>2.){
                        throw("bad gama value");
                }
                return gama;
        }

        vector<vector<double>> fornberg_weights(double z, vector<double> x, int n, int m){
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
                double xmax = 6.*sqrt(tmax*D);
                double dx = fraction*xmax;
                double gama = find_gama(dx, xmax, nx);
                vector<double> x(nx+2);

                for(int i = 0; i<nx+2; i++){
                        x[i] = dx*(pow(gama, i)-1.)/(gama-1.);
                }
                return x;
        }

        vector<double> vslice(vector<double> copy_from, int i1, int i2){
                vector<double> copy; 
                for(int i = i1; i<i2; i++){
                        copy.push_back(copy_from[i]);
                }
                return copy;
        }

        vector<vector<vector<double>>> calc_main_coeficients(vector<double> x, double dt, vector<double> D, int nx, double B){
                vector<double> a1;
                vector<double> a2;
                vector<double> a3;
                vector<double> a4;
                vector<vector<double>> weights;
                vector<double> xinbtw;
                double alfa1d, alfa2d, alfa3d, alfa4d, alfa1v, alfa2v, alfa3v, alfa4v;

                vector<vector<vector<double>>> main_array(nx, vector<vector<double>>(D.size(), vector<double>(4)));

                for (int i = 1; i<nx; i++){
                        xinbtw = vslice(x, i-1, i+3);
                        weights = fornberg_weights(x[i], xinbtw, 4, 2);
                        alfa1d = weights[0][2];
                        alfa2d = weights[1][2];
                        alfa3d = weights[2][2];
                        alfa4d = weights[3][2];

                        alfa1v = -(B*pow(x[i],2))*weights[0][1];
                        alfa2v = -(B*pow(x[i],2))*weights[1][1];
                        alfa3v = -(B*pow(x[i],2))*weights[2][1];
                        alfa4v = -(B*pow(x[i],2))*weights[3][1];

                        for (int j = 0; j<D.size(); j++){
                                main_array[i-1][j][0] = dt*(-alfa1d*D[j] - alfa1v);
                                main_array[i-1][j][1] = dt*(-alfa2d*D[j] - alfa2v)+1.;
                                main_array[i-1][j][2] = dt*(-alfa3d*D[j] - alfa3v);
                                main_array[i-1][j][3] = dt*(-alfa4d*D[j] - alfa4v);
                        }
                }
                return main_array;
        }

        vector<vector<double>> calc_boundary_condition(vector<double> x, int i, vector<double> D, int nx, double B){
                vector<vector<double>> bound_array(D.size(), vector<double>(3));

                double alfa1, alfa2, alfa3;
                vector<double> xinbtw;
                vector<vector<double>> weights;

                if (i==0){
                        xinbtw = vslice(x, 0, 3);
                        weights = fornberg_weights(x[i], xinbtw, 3, 1);
                }else if (i == -1){
                        xinbtw = vslice(x, nx-3, nx);
                        weights = fornberg_weights(x[nx-1], xinbtw, 3, 1);
                }else{
                        throw("Boundary Error: boundary indexed incorrectly");
                }

                alfa1 = weights[0][1] - B*pow(x[i],2.);
                alfa2 = weights[1][1] - B*pow(x[i],2.);
                alfa3 = weights[2][1] - B*pow(x[i],2.);

                for (int i = 0; i<D.size(); i++){
                        bound_array[i][0] = -alfa1*D[i];
                        bound_array[i][1] = -alfa2*D[i];
                        bound_array[i][2] = -alfa3*D[i];
                }
                return bound_array;
        }

        vector<vector<double>> create_constant_list(vector<int> indexs, vector<vector<double>> consts){
                vector<vector<double>> c;
                for(int i = 0; i<indexs.size(); i++){
                        c.push_back(consts[indexs[i]]);
                }
                return c;
        }

        void fcn(int n, double x[], double f[], int &iflag, Params params){
                vector<double> kinetics;
                vector<double> conc_slice;
                vector<double> bound_slice;
                int spec_num;
                double ga;
                bound_slice = params.aslice(x, 0, params.n_ads+params.n_dis);
                kinetics = params.sum_two_vectors(params.calc_kinetics(0, bound_slice, params.isotherm),
                        params.calc_EC_kinetics(2,bound_slice, x[n-2]));
                for(int i = 0; i < params.n_ads; i++){
                        f[i] = (x[i] - params.cp[i])/params.dt*params.eqilib - kinetics[i]*params.spectator[i];
                }
                
                if (params.n_dis > 0){
                        
                        for(int i = params.n_ads; i<params.n_ads+params.n_dis; i++){
                                conc_slice = params.get_conc_at_pos(x, i-params.n_ads, 0, 3, params.n_dis, params.n_ads);
                                f[i] = params.vector_dot_product(params.bound1[i-params.n_ads], conc_slice) - kinetics[i]*params.spectator[i];
                        }

                        for (int xx = 1; xx < params.nx; xx++){
                                conc_slice = params.aslice(x,params.n_ads+params.n_dis*xx, params.n_ads+params.n_dis+xx*params.n_dis);
                                kinetics = params.calc_kinetics(1, conc_slice, params.null);
                                for (int i = params.n_ads+params.n_dis*xx; i < params.n_ads+params.n_dis+xx*params.n_dis; i++){
                                        spec_num = i - params.n_ads - params.n_dis*xx;
                                        conc_slice = params.get_conc_at_pos(x, spec_num,xx-1, xx+3, params.n_dis, params.n_ads);
                                        f[i] = params.vector_dot_product(params.a[xx-1][spec_num], conc_slice)*params.eqilib - 
                                                params.dt*kinetics[spec_num]*params.spectator[params.n_ads+spec_num] - params.cp[i];
                                }
                        }
                        for (int i = n-2*params.n_dis-2; i<n-params.n_dis-2; i++){
                                f[i] = x[i] - params.bound2[i-n+2*params.n_dis+2];
                        }
                        for (int i = n-params.n_dis-2; i<n-2; i++){
                                f[i] = x[i] - x[i-params.n_dis];
                        }
                }
                ga = params.A*F*params.calc_current(2, bound_slice, x[n-2]);

                f[n-2] = params.pnom - x[n-2] - params.Ru*ga - params.Ru*x[n-1];
                f[n-1] = (1+params.Ru*params.Cdl/params.dt)*x[n-1] - params.Cdl*params.delta
                         - params.Ru*params.Cdl*params.cp[n-1]/params.dt; 
                
        }
};

PYBIND11_MODULE(cpp_ekitty_simulator, m){
    py::class_<Electrokitty_simulator>(m, "cpp_ekitty_simulator")
    .def(py::init())
    .def("set_parameters", &Electrokitty_simulator::set_params)
    .def("set_simulation_programm", &Electrokitty_simulator::set_sim_prog)
    .def("create_optimization_problem", &Electrokitty_simulator::create_optimization_problem)
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
    .def("calc_from_guess", [](Electrokitty_simulator &self, vector<double> guess){
        py::array i_sim = py::cast(self.calc_from_guess(guess));
        return i_sim;
    })
    .def_readwrite("current", &Electrokitty_simulator::current)
    .def_readwrite("t", &Electrokitty_simulator::t)
    .def_readwrite("E_generated", &Electrokitty_simulator::E_generated)
    .def_readwrite("concentration_profile", &Electrokitty_simulator::concentration_profile)
    .def_readwrite("surface_profile", &Electrokitty_simulator::surface_profile)
    ;
}
