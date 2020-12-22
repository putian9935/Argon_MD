#include "MD.h"

// hello world!
// is my username correct
//cyc join in

MD_system::MD_system(int N, double temperature, double a, int N_hoover, double dt, int every_save) : N(N), temperature(temperature), a(a), N_hoover(N_hoover), dt(dt), every_save(every_save), stream_opened(false), calculate_pressure(false)
{
    particles = std::vector<Particle>(N);
    daemons = std::vector<Particle>(N_hoover);
    force = std::vector<Force>(N + N_hoover);
    trajectory = std::vector<std::vector<Particle>>(0);
    Q = std::vector<double>(N_hoover);

    // initialize coordinates & momentum
    initialize_system();
}

MD_system::~MD_system()
{
    record_done(); // close fstream
}

void MD_system::initialize_system()
{

    // distribution generators
    std::default_random_engine generator((unsigned)time(NULL));
    std::uniform_real_distribution<double> unif_dist(0, a);
    std::normal_distribution<double> stdnorm_dist(0, 1);
    double sigma_p = k * temperature * m; // TODO: CHANGE THIS

    // particles
    for (int i = 0; i < N; i++)
    {
        particles[i].x = unif_dist(generator);
        particles[i].y = unif_dist(generator);
        particles[i].z = unif_dist(generator);
        particles[i].px = sigma_p * stdnorm_dist(generator);
        particles[i].py = sigma_p * stdnorm_dist(generator);
        particles[i].pz = sigma_p * stdnorm_dist(generator);
    }

    // daemons
    for (int i = 0; i < N_hoover; i++)
    {
        daemons[i].x = unif_dist(generator);
        daemons[i].px = 0;

        Q[i] = k * temperature * 15; // This has to be made clear!!
    }
    Q[0] *= 3 * N;

    // force (should be computed before `update`)
    calculate_force();
}

void MD_system::update()
{
    // REQUIRE: Calculate F(0) first
    // after updating position, move the particle back to the original cube
    for (int i = 0; i < N; i++)
    {
        particles[i].px += dt / 2 * force[i].fx;
        particles[i].x += dt / m * particles[i].px;

        particles[i].py += dt / 2 * force[i].fy;
        particles[i].y += dt / m * particles[i].py;

        particles[i].pz += dt / 2 * force[i].fz;
        particles[i].z += dt / m * particles[i].pz;

        move_back_to_cube(particles[i]);
    }
    for (int i = 0; i < N_hoover; i++)
    {
        daemons[i].px += dt / 2 * force[i + N].fx;
        daemons[i].x += dt / Q[i] * daemons[i].px;
    }

    calculate_force();

    //update p
    for (int i = 0; i < N; i++)
    {
        particles[i].px += dt / 2 * force[i].fx;
        particles[i].py += dt / 2 * force[i].fy;
        particles[i].pz += dt / 2 * force[i].fz;
    }
    for (int i = 0; i < N_hoover; i++)
    {
        daemons[i].px += dt / 2 * force[i + N].fx;
    }
}

void MD_system::move_back_to_cube(Particle &particle)
{
    if ((particle.x < 0) || (particle.x > a))
    {
        int cross_cubes = floor(particle.x / a);
        particle.x -= a * cross_cubes;

        if (calculate_pressure)
        {
            test_counter += 1;
            accumulate_momentum_crossed.px += fabs(particle.px * cross_cubes);
        }
    }

    if ((particle.y < 0) || (particle.y > a))
    {
        int cross_cubes = floor(particle.y / a);
        particle.y -= a * cross_cubes;

        if (calculate_pressure)
        {
            test_counter += 1;
            accumulate_momentum_crossed.py += fabs(particle.py * cross_cubes);
        }
    }

    if ((particle.z < 0) || (particle.z > a))
    {
        int cross_cubes = floor(particle.z / a);
        particle.z -= a * cross_cubes;

        if (calculate_pressure)
        {
            test_counter += 1;
            accumulate_momentum_crossed.pz += fabs(particle.pz * cross_cubes);
        }
    }
}

void MD_system::calculate_force()
{
    clear_force(); // reset to 0

    // interaction between molecules
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            accumulate_pair_force(i, j);
        }
    }
    // add contribution from daemon
    for (int i = 0; i < N; i++)
    {
        accumulate_daemon_force(i);
    }

    compute_dd_force();
}

void MD_system::clear_force()
{
    for (int i = 0; i < N + N_hoover; i++)
    {
        force[i].fx = force[i].fy = force[i].fz = 0;
    }
}

void MD_system::accumulate_daemon_force(int i)
{
    double scale_factor = daemons[0].px / Q[0];
    force[i].fx -= scale_factor * particles[i].px;
    force[i].fy -= scale_factor * particles[i].py;
    force[i].fz -= scale_factor * particles[i].pz;
}

double MD_system::interaction_force(double dist)
{
    // V(r)=4*epsilon*((sigma/r)^12-(sigma/r)^6)
    // F(r)=4*epsilon*(12(sigma/r)^12/r-6(sigma/r)^6/r)
    double epsilon = 1; // TODO: change the parameters
    double sigma = 1;
    double sr6 = std::pow(sigma / dist, 6);

    return epsilon * (48 * sr6 * sr6 / dist - 24 * sr6 / dist);
}

void MD_system::accumulate_pair_force(int i, int j)
{
    double rx = nearest_dist(particles[j].x - particles[i].x);
    double ry = nearest_dist(particles[j].y - particles[i].y);
    double rz = nearest_dist(particles[j].z - particles[i].z);

    double r = std::sqrt(rx * rx + ry * ry + rz * rz);

    if (r < r_cutoff)
    {
        double F = interaction_force(r);
        force[j].fx += F * rx / r;
        force[j].fy += F * ry / r;
        force[j].fz += F * rz / r;

        force[i].fx -= F * rx / r;
        force[i].fy -= F * ry / r;
        force[i].fz -= F * rz / r;
    }
}

void MD_system::compute_dd_force()
{
    // Lecture 6, page 36, Nose-Hoover Chain thermostat
    // Note that for the 1st daemon, the index of it is N instead of 0

    // Step 1. Calc p_{\eta_1}
    force[N].fx = 0.;
    for (int i = 0; i < N; ++i)
    {
        force[N].fx += particles[i].px * particles[i].px + particles[i].py * particles[i].py + particles[i].pz * particles[i].pz; // accumulate squares of momentum first, this might be of use elsewhere, so maybe it's a good idea to save this quantity
    }
    force[N].fx /= m;                                                                // divide by mass at last stage
    force[N].fx -= (3 * N * k * temperature + daemons[0].px * daemons[1].px / Q[1]); // this is definitely a 3-D system

    // Step 2. Calc p_{\eta_j}, with j = 2, ..., N_hoover-1
    for (int j = 1; j < (N_hoover - 1); ++j)
    {
        force[N + j].fx = daemons[j - 1].px * daemons[j - 1].px / Q[j - 1] - daemons[j + 1].px * daemons[j].px / Q[j + 1] - k * temperature;
    }

    // Step 3. Calc p_{\eta_{H_hoover}} (oops, I'm curious how TeX would be typesetting this)
    force[N + N_hoover - 1].fx = daemons[N_hoover - 2].px * daemons[N_hoover - 2].px / Q[N_hoover - 2] - k * temperature;
}

bool MD_system::record_distance()
{
    if (!stream_opened)
    {
        output_stream.open("../distance_stat.txt", std::fstream::app | std::fstream::in);
        stream_opened = true;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            double rx = nearest_dist(particles[j].x - particles[i].x);
            double ry = nearest_dist(particles[j].y - particles[i].y);
            double rz = nearest_dist(particles[j].z - particles[i].z);

            output_stream << std::sqrt(rx * rx + ry * ry + rz * rz) << std::endl;
        }
    }

    return !output_stream.fail();
}

void MD_system::clear_pressure()
{
    accumulate_momentum_crossed.px = 0;
    accumulate_momentum_crossed.py = 0;
    accumulate_momentum_crossed.pz = 0;
}

void MD_system::calculate_auto_correlation(int max_time, char * const file_name)
{
    std::ofstream save_correl_stream(file_name);
    auto_correlation.reserve(max_time);
    int percent = 0;
    printf("Calculating auto-correlation: 0%%");
    for (int tau = 0; tau < max_time; ++tau)
    { // Iterate over different time difference
        if ((tau + 1) * 100 / max_time > percent)
        {
            percent = (tau + 1) * 100 / max_time;
            printf("\rCalculating auto-correlation %d%%    ", percent);
            fflush(stdout);
        }
        auto_correlation.push_back(0.); // clear
        for (int t = 0; t < trajectory.size() - tau; ++t)
        { // Iterate over time
            for (int j = 0; j < N; ++j)
            { // Iterate over particles
                auto_correlation[tau] += (trajectory[t][j].px * trajectory[t + tau][j].px + trajectory[t][j].py * trajectory[t + tau][j].py + trajectory[t][j].pz * trajectory[t + tau][j].pz);
            }
        }
        auto_correlation[tau] /= (3. * N * (trajectory.size() - tau));
        save_correl_stream << auto_correlation[tau] << '\n';
    }

    printf("\n");
    save_correl_stream.close();
}

void MD_system::append_current_state()
{
    trajectory.push_back(particles);
}
// non class members

Particle get_pressure(MD_system &sys, int init_steps, int simulation_steps)
{

    sys.calculate_pressure = false;

    int percent = 0;
    printf("Calculating pressure:----------------------------------------\n");
    printf("preparing: 0%%");
    for (int i = 0; i < init_steps; i++)
    {
        sys.update();
        if ((i + 1) * 100 / init_steps > percent)
        {
            percent = (i + 1) * 100 / init_steps;
            printf("\rpreparing %d%%    ", percent);
            fflush(stdout);
        }
    }
    printf("\n");

    sys.clear_pressure();
    sys.calculate_pressure = true;

    printf("simulating: 0%%");
    percent = 0;
    for (int i = 0; i < simulation_steps; i++)
    {
        sys.update();
        if (!(i % sys.every_save))  // copy less frequently
            sys.append_current_state();
        if ((i + 1) * 100 / simulation_steps > percent)
        {
            percent = (i + 1) * 100 / simulation_steps;
            printf("\rsimulating %d%%   ", percent);
            fflush(stdout);
        }
    }
    printf("\n");

    Particle pressure;
    pressure.px = sys.accumulate_momentum_crossed.px / 2 / sys.a / sys.a / sys.dt / simulation_steps;
    pressure.py = sys.accumulate_momentum_crossed.py / 2 / sys.a / sys.a / sys.dt / simulation_steps;
    pressure.pz = sys.accumulate_momentum_crossed.pz / 2 / sys.a / sys.a / sys.dt / simulation_steps;

    sys.calculate_pressure = false;

    std::cout << "The pressure of the system is: " << pressure.px << ";" << pressure.py << ";" << pressure.pz << std::endl;
    std::cout << "(crossed " << sys.test_counter << " times in total)" << std::endl;
    printf("-------------------------------------------------------------\n");

    return pressure;
}
