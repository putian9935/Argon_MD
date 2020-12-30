#include "MD.h"
#include "helper_funcs.h"
// hello world!
// is my username correct
//cyc join in


MD_system::MD_system(int N, double temperature, double a, int N_hoover, double dt, int every_save, bool shift_momentum)
    : N(N), temperature(temperature), a(a), N_hoover(N_hoover), dt(dt), every_save(every_save),
      shift_momentum(shift_momentum), stream_opened(false), calculate_pressure(false),
      has_velocity_auto_correlation_calced(false), has_stress_tensor_auto_correlation_calced(false), has_heat_flux_auto_correlation_calced(false)
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
    double sigma_p = sqrt(k * temperature * m);

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
    if (shift_momentum)
    {
        double average_px = 0., average_py = 0., average_pz = 0.;
        for (const auto &particle : particles)
        {
            average_px += particle.px;
            average_py += particle.py;
            average_pz += particle.pz;
        }
        average_px /= N;
        average_py /= N;
        average_pz /= N;
        for (auto &particle : particles)
        {
            particle.px -= average_px;
            particle.py -= average_py;
            particle.pz -= average_pz;
        }
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
    double sr6 = std::pow(1.0 / dist, 6);

    return (48 * sr6 * sr6 / dist - 24 * sr6 / dist);
}

double MD_system::interaction_potential(double dist)
{
    // V(r)=4*epsilon*((sigma/r)^12-(sigma/r)^6)
    double sr6 = std::pow(1.0 / dist, 6);
    return 4 * sr6 * (sr6 - 1.0);
}

void MD_system::accumulate_pair_force(int i, int j)
{
    double rx = nearest_dist(particles[j].x - particles[i].x);
    double ry = nearest_dist(particles[j].y - particles[i].y);
    double rz = nearest_dist(particles[j].z - particles[i].z);

    double r = std::sqrt(rx * rx + ry * ry + rz * rz);

    if (r < r_cutoff_big)
    {
        double F = interaction_force(std::max(r_cutoff_small, r));
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

void MD_system::append_current_state()
{
    accumulate_full_pair_force_and_potential();
    stress_tensor_traj.push_back(accumulate_stress_tensor());
    heat_flux_traj.push_back(accumulate_heat_flux());
    trajectory.push_back(particles);
}

void MD_system::accumulate_full_pair_force_and_potential()
{
    // accumulate the force and potential of each particle without any cutoff
    // expensive function, yet called less frequently
    full_pair_force.reserve(N);
    full_potential.reserve(N);
    for (int i = 0; i < N; ++i)
    {
        full_pair_force[i].reserve(N);
        full_potential[i] = 0.;
        for (int j = 0; j < N; ++j)
            full_pair_force[i][j].fx = full_pair_force[i][j].fy = full_pair_force[i][j].fz = 0.;
    }

    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            double rx = nearest_dist(particles[j].x - particles[i].x);
            double ry = nearest_dist(particles[j].y - particles[i].y);
            double rz = nearest_dist(particles[j].z - particles[i].z);

            double r = std::sqrt(rx * rx + ry * ry + rz * rz);
            double F = interaction_force(r);
            full_pair_force[j][i].fx = F * rx / r;
            full_pair_force[j][i].fy = F * ry / r;
            full_pair_force[j][i].fz = F * rz / r;

            full_pair_force[i][j].fx = -full_pair_force[j][i].fx;
            full_pair_force[i][j].fy = -full_pair_force[j][i].fy;
            full_pair_force[i][j].fz = -full_pair_force[j][i].fz;

            full_potential[i] += interaction_potential(r);
            full_potential[j] += interaction_potential(r);
        }
    }
}


Particle MD_system::accumulate_stress_tensor()
{
    // The three components represents: sigma_{xy}, sigma_{yz}, and sigma_{zx}, respectively
    Particle ret;
    ret.x = ret.y = ret.z = 0.;
    for (int i = 0; i < N; ++i)
    {
        ret.x += (particles[i].px * particles[i].py);
        ret.y += (particles[i].py * particles[i].pz);
        ret.z += (particles[i].pz * particles[i].px);
    }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N ;++j)
        {
            // if (j == i) continue;
            ret.x += .5 * (particles[i].x - particles[j].x) * full_pair_force[i][j].fy;
            ret.y += .5 * (particles[i].y - particles[j].y) * full_pair_force[i][j].fz;
            ret.z += .5 * (particles[i].z - particles[j].z) * full_pair_force[i][j].fx;
        }
    return ret;
}


Particle MD_system::accumulate_heat_flux()
{
    // Heat flux consists of two parts: 1. particle energy 2. work done by force
    // The three components represents: sigma_{xy}, sigma_{yz}, and sigma_{zx}, respectively
    Particle ret;
    ret.x = ret.y = ret.z = 0.;

    // 1. calculate particle energy
    for (int i = 0; i < N; ++i)
    {
        double e = .5 * (particles[i].px * particles[i].px + particles[i].py * particles[i].py + particles[i].pz * particles[i].pz) + .5 * full_potential[i];
        ret.x += e * particles[i].px;
        ret.y += e * particles[i].py;
        ret.z += e * particles[i].pz;
    }

    // 2. calculate work done by force
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            double buf = (full_pair_force[i][j].fx * particles[i].px + full_pair_force[i][j].fy * particles[i].py + full_pair_force[i][j].fz * particles[i].pz) / 2.;
            ret.x += buf * (particles[i].x - particles[j].x);
            ret.y += buf * (particles[i].y - particles[j].y);
            ret.z += buf * (particles[i].z - particles[j].z);
        }
    }
    return ret;
}

void MD_system::calculate_velocity_auto_correlation(int max_time, const char *const file_name)
{
    fprintf(stderr, "----------------------------------------\nCalculating velocity auto-correlation:\n");
    // Calculate auto-correlation at several time difference tau
    std::ofstream save_correl_stream(file_name);
    velocity_auto_correlation.reserve(max_time);

    int percent = 0;
    fprintf(stderr, "Processing: 0%%");
    for (int tau = 0; tau < max_time; ++tau)
    { // Iterate over different time difference
        if ((tau + 1) * 100 / max_time > percent)
        {
            percent = (tau + 1) * 100 / max_time;
            fprintf(stderr, "\rProcessing: %d%%    ", percent);
            fflush(stderr);
        }
        velocity_auto_correlation.push_back(0.); // clear
        for (int t = 0; t < trajectory.size() - tau; ++t)
        { // Iterate over time
            for (int j = 0; j < N; ++j)
            { // Iterate over particles
                velocity_auto_correlation[tau] += (trajectory[t][j].px * trajectory[t + tau][j].px + trajectory[t][j].py * trajectory[t + tau][j].py + trajectory[t][j].pz * trajectory[t + tau][j].pz);
            }
        }
        velocity_auto_correlation[tau] /= (3. * N * (trajectory.size() - tau));
        save_correl_stream << std::setprecision(17) << velocity_auto_correlation[tau] << '\n';
    }
    fprintf(stderr, "\n");
    save_correl_stream.close();
    has_velocity_auto_correlation_calced = true;
    fprintf(stderr, "----------------------------------------\n");
}

double MD_system::calculate_self_diffusion_constant(bool use_coarse_estimate, int cut_off)
{
    // Integrate auto-correlation over tau
    if (!has_velocity_auto_correlation_calced)
        calculate_velocity_auto_correlation();

    double rescaled_dt = dt * every_save * time_conversion_constant; // Difference between tau's
    if (use_coarse_estimate)
    { // use exponential decay to estimate
        printf("Using coarse estimate...\n");
        return velocity_auto_correlation[0] / (log(velocity_auto_correlation[0]) - log(velocity_auto_correlation[1])) * rescaled_dt * velocity_conversion_constant * velocity_conversion_constant;
    }

    if (cut_off > velocity_auto_correlation.size())
    { // too large, give a comparison
        printf("Cut off set too large, really should use coarse estimate, which yields %.6e\n", velocity_auto_correlation[0] / (log(velocity_auto_correlation[0]) - log(velocity_auto_correlation[1])) * rescaled_dt * velocity_conversion_constant * velocity_conversion_constant);
        cut_off = velocity_auto_correlation.size();
    }

    printf("Using Simpson rule...\n");
    // In order to use Simpson rule, total tau's must be an odd number
    cut_off -= (cut_off & 1);
    double self_diffusion_constant = velocity_auto_correlation[0] + velocity_auto_correlation[cut_off - 1];
    for (int i = 1; i < cut_off - 1; ++i)
        self_diffusion_constant += (2 << (i & 1)) * velocity_auto_correlation[i];

    return self_diffusion_constant * rescaled_dt / 3. * velocity_conversion_constant * velocity_conversion_constant; // convert velocity into SI units!
}

void MD_system::calculate_stress_tensor_auto_correlation(int max_time, const char *const file_name)
{
    fprintf(stderr, "----------------------------------------\nCalculating stress-tensor auto-correlation:\n");
    // Calculate auto-correlation at several time difference tau
    std::ofstream save_correl_stream(file_name);

    stress_tensor_auto_correlation.reserve(max_time);

    int percent = 0;

    auto bufx = 0., bufy = 0., bufz = 0.;
    for (auto &x : stress_tensor_traj)
    {
        bufx += x.x;
        bufy += x.y;
        bufz += x.z;
    }
    bufx /= stress_tensor_traj.size();
    bufy /= stress_tensor_traj.size();
    bufz /= stress_tensor_traj.size();
    for (auto &x : stress_tensor_traj)
    {
        x.x -= bufx;
        x.y -= bufy;
        x.z -= bufz;
    }

    fprintf(stderr, "Processing: 0%%");

    for (int tau = 0; tau < max_time; ++tau)
    { // Iterate over different time difference
        if ((tau + 1) * 100 / max_time > percent)
        {
            percent = (tau + 1) * 100 / max_time;
            fprintf(stderr, "\rProcessing: %d%%    ", percent);
            fflush(stderr);
        }
        stress_tensor_auto_correlation.push_back(0.); // clear
        for (int t = 0; t < stress_tensor_traj.size() - tau; ++t)
        { // Iterate over time
            for (int j = 0; j < N; ++j)
            { // Iterate over particles
                stress_tensor_auto_correlation[tau] += (stress_tensor_traj[t].x * stress_tensor_traj[t + tau].x + stress_tensor_traj[t].y * stress_tensor_traj[t + tau].y + stress_tensor_traj[t].z * stress_tensor_traj[t + tau].z);
            }
        }
        stress_tensor_auto_correlation[tau] /= (3. * (stress_tensor_traj.size() - tau));
        save_correl_stream << std::setprecision(17) << stress_tensor_auto_correlation[tau] << '\n';
    }
    fprintf(stderr, "\n");
    save_correl_stream.close();
    has_stress_tensor_auto_correlation_calced = true;
    fprintf(stderr, "----------------------------------------\n");
}

double MD_system::calculate_shear_viscosity_coefficient(bool use_coarse_estimate, int cut_off)
{
    // Integrate auto-correlation over tau
    if (!has_stress_tensor_auto_correlation_calced)
        calculate_stress_tensor_auto_correlation();

    double rescaled_dt = dt * every_save * time_conversion_constant; // Difference between tau's
    if (use_coarse_estimate)
    { // use exponential decay to estimate
        printf("Using coarse estimate...\n");
        return stress_tensor_auto_correlation[0] / (log(stress_tensor_auto_correlation[0]) - log(stress_tensor_auto_correlation[1])) * rescaled_dt / pow(a, 3) / temperature * pressure_conversion_constant;
    }

    if (cut_off > stress_tensor_auto_correlation.size())
    { // too large, give a comparison
        printf("Cut off set too large, really should use coarse estimate, which yields %.6e\n", stress_tensor_auto_correlation[0] / (log(stress_tensor_auto_correlation[0]) - log(stress_tensor_auto_correlation[1])) * rescaled_dt / pow(a, 3) / temperature * pressure_conversion_constant);
        cut_off = stress_tensor_auto_correlation.size();
    }

    printf("Using Simpson rule...\n");
    // In order to use Simpson rule, total tau's must be an odd number
    cut_off -= (cut_off & 1);
    double shear_viscosity_coefficient = stress_tensor_auto_correlation[0] + stress_tensor_auto_correlation[cut_off - 1];
    for (int i = 1; i < cut_off - 1; ++i)
    {
        shear_viscosity_coefficient += (2 << (i & 1)) * stress_tensor_auto_correlation[i];
    }

    return shear_viscosity_coefficient * rescaled_dt / 3. / pow(a, 3) / temperature * pressure_conversion_constant; // convert viscosity into SI units!
}

void MD_system::calculate_heat_flux_auto_correlation(int max_time, const char *const file_name)
{
    fprintf(stderr, "----------------------------------------\nCalculating heat-flux auto-correlation:\n");
    // Calculate auto-correlation at several time difference tau
    std::ofstream save_correl_stream(file_name);
    heat_flux_auto_correlation.reserve(max_time);

    int percent = 0;

    auto bufx = 0., bufy = 0., bufz = 0.;
    for (auto &x : heat_flux_traj)
    {
        bufx += x.x;
        bufy += x.y;
        bufz += x.z;
    }
    bufx /= heat_flux_traj.size();
    bufy /= heat_flux_traj.size();
    bufz /= heat_flux_traj.size();
    for (auto &x : heat_flux_traj)
    {
        x.x -= bufx;
        x.y -= bufy;
        x.z -= bufz;
    }

    fprintf(stderr, "Processing: 0%%");

    for (int tau = 0; tau < max_time; ++tau)
    { // Iterate over different time difference
        if ((tau + 1) * 100 / max_time > percent)
        {
            percent = (tau + 1) * 100 / max_time;
            fprintf(stderr, "\rProcessing: %d%%    ", percent);
            fflush(stderr);
        }
        heat_flux_auto_correlation.push_back(0.); // clear
        for (int t = 0; t < heat_flux_traj.size() - tau; ++t)
        { // Iterate over time
            for (int j = 0; j < N; ++j)
            { // Iterate over particles
                heat_flux_auto_correlation[tau] += (heat_flux_traj[t].x * heat_flux_traj[t + tau].x + heat_flux_traj[t].y * heat_flux_traj[t + tau].y + heat_flux_traj[t].z * heat_flux_traj[t + tau].z);
            }
        }
        heat_flux_auto_correlation[tau] /= (3. * (heat_flux_traj.size() - tau));
        save_correl_stream << std::setprecision(17) << heat_flux_auto_correlation[tau] << '\n';
    }
    fprintf(stderr, "\n");
    save_correl_stream.close();
    has_heat_flux_auto_correlation_calced = true;
    fprintf(stderr, "----------------------------------------\n");
}

double MD_system::calculate_thermal_conductivity(bool use_coarse_estimate, int cut_off)
{
    // Integrate auto-correlation over tau
    if (!has_heat_flux_auto_correlation_calced)
        calculate_heat_flux_auto_correlation();

    double rescaled_dt = dt * every_save; // Difference between tau's
    if (use_coarse_estimate)
    { // use exponential decay to estimate
        printf("Using coarse estimate...\n");
        return heat_flux_auto_correlation[0] / (log(heat_flux_auto_correlation[0]) - log(heat_flux_auto_correlation[1])) * rescaled_dt / pow(a, 3) / temperature / temperature * (pressure_conversion_constant * length_conversion_constant * length_conversion_constant / (temperature_conversion_constant * time_conversion_constant ));
    }

    if (cut_off > heat_flux_auto_correlation.size())
    { // too large, give a comparison
        printf("Cut off set too large, really should use coarse estimate, which yields %.6e\n", heat_flux_auto_correlation[0] / (log(heat_flux_auto_correlation[0]) - log(heat_flux_auto_correlation[1])) * rescaled_dt / pow(a, 3) / temperature / temperature * (pressure_conversion_constant * length_conversion_constant * length_conversion_constant / (temperature_conversion_constant * time_conversion_constant )));
    }

    printf("Using Simpson rule...\n");
    // In order to use Simpson rule, total tau's must be an odd number
    cut_off -= (cut_off & 1);
    double thermal_conductivity = heat_flux_auto_correlation[0] + heat_flux_auto_correlation[cut_off - 1];
    for (int i = 1; i < cut_off - 1; ++i)
        thermal_conductivity += (2 << (i & 1)) * heat_flux_auto_correlation[i];

    return thermal_conductivity * rescaled_dt / 3. / pow(a, 3) / temperature / temperature * (pressure_conversion_constant * length_conversion_constant * length_conversion_constant / (temperature_conversion_constant * time_conversion_constant )); // convert viscosity into SI units!
}

double MD_system::get_temperature()
{
    double p_squared = 0.0;
    for (Particle &par : particles)
    {
        p_squared += par.px * par.px + par.py * par.py + par.pz * par.pz;
    }
    return p_squared / (3 * m * N * k);
}

// Conversion constants
// please refer for details to http://cms.sjtu.edu.cn/doc/courseware/2019/LJ.pdf
const double MD_system::pressure_conversion_constant = 4.2e7;
const double MD_system::time_conversion_constant = 2.17e-12;
const double MD_system::velocity_conversion_constant = 1.57e2;
const double MD_system::temperature_conversion_constant = 1.2e2;
const double MD_system::length_conversion_constant = 3.4e-10;

// non class members

Particle get_pressure_collision(MD_system &sys, int init_steps, int simulation_steps)
{
    printf("The simulation runs for %.1f ps in total, with %.1f ps burn-in\n", ((init_steps + simulation_steps) * sys.dt * MD_system::time_conversion_constant) / 1e-12, (init_steps * sys.dt * MD_system::time_conversion_constant) / 1e-12);
    printf("system parameter: T=%.1f K, V=(%.1f A)^3\n", sys.temperature * sys.temperature_conversion_constant, 1e10 * sys.a * sys.length_conversion_constant);
    sys.calculate_pressure = false;

    int percent = 0;
    printf("Calculating pressure(collision):----------------------------------------\n");
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
        if ((i + 1) * 100 / simulation_steps > percent)
        {
            percent = (i + 1) * 100 / simulation_steps;
            printf("\rsimulating %d%%   ", percent);
            fflush(stdout);
        }
    }
    printf("\n");

    Particle pressure;
    pressure.px = sys.accumulate_momentum_crossed.px / sys.a / sys.a / sys.dt / simulation_steps;
    pressure.py = sys.accumulate_momentum_crossed.py / sys.a / sys.a / sys.dt / simulation_steps;
    pressure.pz = sys.accumulate_momentum_crossed.pz / sys.a / sys.a / sys.dt / simulation_steps;

    sys.calculate_pressure = false;

    std::cout << "The pressure of the system is: " << pressure.px << ";" << pressure.py << ";" << pressure.pz << std::endl;
    std::cout << "Averaged pressure: " << (pressure.px + pressure.py + pressure.pz) / 3 * sys.pressure_conversion_constant << "Pa" << std::endl;
    std::cout << "(crossed " << sys.test_counter << " times in total)" << std::endl;
    printf("-------------------------------------------------------------\n\n");

    return pressure;
}

double MD_system::pressure_virial()
{
    double pressure_i = 0;
    for (int j = 0; j < N; j++)
    {
        pressure_i += particles[j].px * particles[j].px + particles[j].py * particles[j].py + particles[j].pz * particles[j].pz;

        for (int k = 0; k < j; k++)
        {
            double rx = nearest_dist(particles[j].x - particles[k].x);
            double ry = nearest_dist(particles[j].y - particles[k].y);
            double rz = nearest_dist(particles[j].z - particles[k].z);

            double r = std::sqrt(rx * rx + ry * ry + rz * rz);
            if (r < r_cutoff_big)
            {
                pressure_i += r * interaction_force(std::max(r_cutoff_small, r));
            }
        }
    }
    pressure_i /= 3 * a * a * a;

    return pressure_i;
}

Particle get_pressure_virial(MD_system &sys, int init_steps, int simulation_steps)
{
    printf("The simulation runs for %.1f ps in total, with %.1f ps burn-in\n", ((init_steps + simulation_steps) * sys.dt * MD_system::time_conversion_constant) / 1e-12, (init_steps * sys.dt * MD_system::time_conversion_constant) / 1e-12);
    printf("system parameter: T=%.1f K, V=(%.1f A)^3\n", sys.temperature * sys.temperature_conversion_constant, 1e10 * sys.a * sys.length_conversion_constant);
    sys.calculate_pressure = false;

    int percent = 0;
    printf("Calculating pressure(viral):----------------------------------------\n");
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

    double sum_pressure = 0;
    double sum_pressure2 = 0;

    printf("simulating: 0%%");
    percent = 0;
    for (int i = 0; i < simulation_steps; i++)
    {
        sys.update();

        double pressure_i = sys.pressure_virial();
        sum_pressure += pressure_i;
        sum_pressure2 += pressure_i * pressure_i;

        if ((i + 1) * 100 / simulation_steps > percent)
        {
            percent = (i + 1) * 100 / simulation_steps;
            printf("\rsimulating %d%%   ", percent);
            fflush(stdout);
        }
    }
    printf("\n");

    Particle result;
    result.px = sum_pressure / simulation_steps;
    result.py = sum_pressure2 / simulation_steps - result.px * result.px;

    std::cout << "Averaged pressure: " << result.px * sys.pressure_conversion_constant << "Pa; "
              << "std(P): " << std::sqrt(result.py) * sys.pressure_conversion_constant << "Pa; " << std::endl;
    std::cout << "(simulation steps: " << simulation_steps << ")" << std::endl;
    printf("-------------------------------------------------------------\n\n");

    return result;
}

/*
Same parameter signature and meaning as above, print several transport properties.
Still, pressure is returned.
*/
Particle calculate_transport_properties(MD_system &sys, int init_steps, int simulation_steps, double burn_in_dt)
{
    printf("The simulation runs for %.1f ps in total, with %.1f ps burn-in\n", ((init_steps * burn_in_dt + simulation_steps * sys.dt)  * MD_system::time_conversion_constant) / 1e-12, (init_steps * burn_in_dt * MD_system::time_conversion_constant) / 1e-12);
    printf("System parameters: T=%.1f K, V=(%.1f A)^3\n", sys.temperature * sys.temperature_conversion_constant, 1e10 * sys.a * sys.length_conversion_constant);
    sys.calculate_pressure = false;

    int percent = 0;
    fprintf(stderr, "Calculating pressure:----------------------------------------\n");
    fprintf(stderr, "preparing: 0%%");
    double simulation_dt = sys.dt; 
    sys.dt = burn_in_dt;
    for (int i = 0; i < init_steps; i++)
    {
        sys.update();
        if ((i + 1) * 100 / init_steps > percent)
        {
            percent = (i + 1) * 100 / init_steps;
            fprintf(stderr, "\rpreparing %d%%    ", percent);
            fflush(stderr);
        }
    }
    sys.dt = simulation_dt;
    fprintf(stderr, "\n");

    sys.clear_pressure();
    sys.calculate_pressure = true;

    fprintf(stderr, "simulating: 0%%");
    percent = 0;
    for (int i = 0; i < simulation_steps; i++)
    {
        sys.update();
        if (!(i % sys.every_save)) // copy less frequently
            sys.append_current_state();
        if ((i + 1) * 100 / simulation_steps > percent)
        {
            percent = (i + 1) * 100 / simulation_steps;
            fprintf(stderr, "\rsimulating %d%%   ", percent);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\n");

    Particle pressure;
    pressure.px = sys.accumulate_momentum_crossed.px / sys.a / sys.a / sys.dt / simulation_steps;
    pressure.py = sys.accumulate_momentum_crossed.py / sys.a / sys.a / sys.dt / simulation_steps;
    pressure.pz = sys.accumulate_momentum_crossed.pz / sys.a / sys.a / sys.dt / simulation_steps;

    sys.calculate_pressure = false;

    std::cout << "The pressure of the system is: " << pressure.px << ";" << pressure.py << ";" << pressure.pz << std::endl;
    std::cout << "Averaged pressure: " << (pressure.px + pressure.py + pressure.pz) / 3 * sys.pressure_conversion_constant << "Pa" << std::endl;
    std::cout << "(crossed " << sys.test_counter << " times in total)" << std::endl;
    printf("-------------------------------------------------------------\n\n");

    return pressure;
}
