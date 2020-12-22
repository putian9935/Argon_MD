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

struct Particle
{
    double x,y,z;
    double px,py,pz;
};

struct Force
{
    double fx,fy,fz;
};


class MD_system
{
public:
    //parameters of system
    double temperature;
    double a; // length of box
    int N; //number of particles
    int N_hoover;
    static constexpr double k=1;

    std::vector<Particle> particles;
    std::vector<Force> force;
    std::vector<Particle> daemons;

    std::vector<double> Q;//daemons mass
    double m=1;//particle mass

    //parameters for solution
    double dt;

    double r_cutoff=5.0;//cut off when calculating forces among particles

    void initialize_system();
    void update();
    void calculate_force();

    MD_system(int N, double temperature, double a, int N_hoover=3, double dt=0.01,int=200);
    ~MD_system();


    // temperature
    double get_temperature();

    // radial distribution
    bool record_distance();
    void record_done() {output_stream.close();};

    // diffusion
    void calculate_auto_correlation(int=500,char * const="correlation.dat");

    // pressure
    bool calculate_pressure;
    Particle accumulate_momentum_crossed; //using px,py,pz only
    void clear_pressure();
    int test_counter=0;

    
    // Auto-correlation, etc.
    std::vector<std::vector<Particle>> trajectory;
    std::vector<double> auto_correlation;
    int every_save;
    void append_current_state();

private:

    std::fstream output_stream;
    bool stream_opened;

    // maintain periodic boundary condition
    void move_back_to_cube(Particle& particle);

    // nearest distance (along a specific axis) under periodic boundary condtion
    inline double nearest_dist(double orig_dist) {
        return orig_dist > (a/2.0) ? orig_dist-a : (orig_dist < (-a/2.0) ? orig_dist + a : orig_dist);
    }

    // functions to compute forces (particles & daemon)
    void clear_force();
    void accumulate_pair_force(int i, int j);
    void accumulate_daemon_force(int i);
    void compute_dd_force();
    double interaction_force(double dist);


};

Particle get_pressure(MD_system& sys, int init_steps, int simulation_steps);

#endif // MD_H_INCLUDED
