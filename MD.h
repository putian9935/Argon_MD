#ifndef MD_H_INCLUDED
#define MD_H_INCLUDED

#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <iomanip>

struct Particle
{
    double x, y, z;
    double px, py, pz;
};

struct Force
{
    double fx, fy, fz;
};

class MD_system
{
public:
    //parameters of system
    double temperature;
    double a; // length of box
    int N;    //number of particles
    int N_hoover;
    static constexpr double k = 1;

    std::vector<Particle> particles;
    std::vector<Force> force;
    std::vector<Particle> daemons;

    std::vector<double> Q; //daemons mass
    double m = 1;          //particle mass

    //parameters for solution
    double dt;

    double r_cutoff_big = 5.0;   //cut off when calculating forces among particles
    double r_cutoff_small = 0.9; //cut off when calculating forces among particles

    void initialize_system();
    void update();
    void calculate_force();

    MD_system(int N, double temperature, double a, int N_hoover = 3, double dt = 0.01, int = 200, bool = false);
    ~MD_system();

    // temperature
    double get_temperature();
    bool shift_momentum;

    // radial distribution
    bool record_distance();
    void record_done() { output_stream.close(); };

    // pressure
    bool calculate_pressure;
    Particle accumulate_momentum_crossed; //using px,py,pz only
    void clear_pressure();
    int test_counter = 0;
    double pressure_virial();

    void append_current_state();

    // Auto-correlation, diffusion, etc.
    std::vector<std::vector<Particle>> trajectory;
    std::vector<double> velocity_auto_correlation;
    int every_save;
    bool has_velocity_auto_correlation_calced;
    void calculate_velocity_auto_correlation(int = 500, const char *const = "correlation.dat");
    double calculate_self_diffusion_constant(bool = true, int = -1);

    // Viscosity
    std::vector<std::vector<Force>> full_pair_force; // fpf[i][j] from j to i
    std::vector<double> full_potential;
    void accumulate_full_pair_force_and_potential();
    std::vector<Particle> stress_tensor_traj;
    std::vector<double> stress_tensor_auto_correlation;
    bool has_stress_tensor_auto_correlation_calced;
    void calculate_stress_tensor_auto_correlation(int = 500, const char *const = "viscosity_correlation.dat");
    double calculate_shear_viscosity_coefficient(bool = true, int = -1);

    // Thermal conductivity
    std::vector<Particle> heat_flux_traj;
    std::vector<double> heat_flux_auto_correlation;
    bool has_heat_flux_auto_correlation_calced;
    void calculate_heat_flux_auto_correlation(int = 500, const char *const = "heat_flux_correlation.dat");
    double calculate_thermal_conductivity(bool = true, int = -1);

    // Some conversion constants
    static const double pressure_conversion_constant, time_conversion_constant, velocity_conversion_constant,
        temperature_conversion_constant, length_conversion_constant;

private:
    std::fstream output_stream;
    bool stream_opened;

    // maintain periodic boundary condition
    void move_back_to_cube(Particle &particle);

    // nearest distance (along a specific axis) under periodic boundary condition
    inline double nearest_dist(double orig_dist)
    {
        return orig_dist > (a / 2.0) ? orig_dist - a : (orig_dist < (-a / 2.0) ? orig_dist + a : orig_dist);
    }

    // functions to compute forces (particles & daemon)
    void clear_force();
    void accumulate_pair_force(int i, int j);
    void accumulate_daemon_force(int i);
    void compute_dd_force();
    double interaction_force(double dist);
    double interaction_potential(double dist);

    Particle accumulate_heat_flux();
    Particle accumulate_stress_tensor();
};

Particle get_pressure_collision(MD_system &sys, int init_steps, int simulation_steps);
Particle get_pressure_virial(MD_system &sys, int init_steps, int simulation_steps); //px=P,py=var(P)

Particle calculate_transport_properties(MD_system &, int, int, double=0.01);

#endif // MD_H_INCLUDED
