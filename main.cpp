#include "util.hpp"
#include "InputData.hpp"
#include "force.hpp"
#include "ensemble.hpp"

#include <fstream>
#include <ctime>
#include <string> 

using namespace cl::sycl;


void writeData(buffer<Atom> &atoms_buf, int step)
{
	auto atoms = atoms_buf.get_access<access::mode::read>();
	std::string fname("./log/out_"+std::to_string(step)+".dat");
    std::ofstream outputFile(fname, std::ios::out);

    if (!outputFile) {
        std::cerr << "Failed to open the file for writing!" << std::endl;
        return ;
    }

    outputFile << "IDX  x  y  z  vx  vy  vz  fx  fy  fz\n";
    
    for (int i = 0; i < atoms.size(); i++) {
		outputFile << i << ", " 
					<< atoms[i].position[0] << ", " 
					<< atoms[i].position[1]  << ", " 
					<< atoms[i].position[2] << ", " 
					<< atoms[i].velocity[0] << ", " 
					<< atoms[i].velocity[1]  << ", " 
					<< atoms[i].velocity[2] << ", " 
					<< atoms[i].force[0] << ", " 
					<< atoms[i].force[1]  << ", " 
					<< atoms[i].force[2] << std::endl;
    }


    outputFile.close();


    fname = "./log/vel_"+std::to_string(step)+".dat";
    std::ofstream outputFile2(fname, std::ios::out);

    if (!outputFile2) {
        std::cerr << "Failed to open the file for writing!" << std::endl;
        return ;
    }

    outputFile2 << "IDX  vx  vy  vz\n";
    
    for (int i = 0; i < atoms.size(); i++) {
		outputFile2 << i << " " 
					<< atoms[i].velocity[0] << " " 
					<< atoms[i].velocity[1]  << " " 
					<< atoms[i].velocity[2] << std::endl;
    }


    outputFile2.close();


    fname = "./log/force_"+std::to_string(step)+".dat";
    std::ofstream outputFile3(fname, std::ios::out);

    if (!outputFile3) {
        std::cerr << "Failed to open the file for writing!" << std::endl;
        return ;
    }

    outputFile3 << "IDX  fx  fy  fz\n";
    
    for (int i = 0; i < atoms.size(); i++) {
		outputFile3 << i << " " 
					<< atoms[i].force[0] << " " 
					<< atoms[i].force[1]  << " " 
					<< atoms[i].force[2] << std::endl;
    }


    outputFile3.close();
}

void write_xyz(buffer<Atom> &atoms_buf, int step)
{
	auto atoms = atoms_buf.get_access<access::mode::read>();
	std::string fname("./log/trj_"+std::to_string(step)+".xyz");
    std::ofstream outputFile(fname, std::ios::out);

    if (!outputFile) {
        std::cerr << "Failed to open the file for writing!" << std::endl;
        return ;
    }

    outputFile << atoms.size() << "\n";
    outputFile << "\n";
    
    for (int i = 0; i < atoms.size(); i++) {
		outputFile << "C" << "\t" 
					<< atoms[i].position[0] << "\t" 
					<< atoms[i].position[1]  << "\t" 
					<< atoms[i].position[2] << "\t" << std::endl;
    }


    outputFile.close();

}

void DisplayData(queue &q, buffer<Atom> &atoms, SimInfo sim_info, int step) {
	int numAtoms = sim_info.numAtoms;
	double pe = 0.0, ke = 0.0;

    // Apply PBC to atom positions

	auto acc = atoms.get_access<access::mode::read>();

    for(int idx=0; idx < numAtoms; idx++){
		double v = sqrt(acc[idx].velocity[0]*acc[idx].velocity[0]
					+ acc[idx].velocity[1]*acc[idx].velocity[1]
					+ acc[idx].velocity[2]*acc[idx].velocity[2]);
		
		ke += 0.5*acc[idx].mass*v*v;
		pe += acc[idx].pe;
	}

    std::cout<<"== " << step << "\t" << ke << "\t" << pe << std::endl;
}


void ComputeForces(queue &q, buffer<Atom> &atoms, buffer<Bond> &bonds, buffer<Angle> &angles, buffer<Dihedral> &dihedrals, EnvInfo env_info, SimInfo sim_info) {
	int numAtoms = sim_info.numAtoms;

    // First clear the force on each atom
    q.submit([&](handler &cgh) {
        auto acc = atoms.get_access<access::mode::read_write>(cgh);
        
        cgh.parallel_for(range<1>(numAtoms), [=](id<1> idx) {
            acc[idx].force = {0, 0, 0};
            acc[idx].pe = 0.0;
        });
    }).wait_and_throw();
		
	//AngleKernel(q, atoms, angles, sim_info);
	//DihedralKernel(q, atoms, dihedrals, sim_info);
	LennardJonesKernel(q, atoms, env_info, sim_info);
	//EwaldKernel(q, atoms, numAtoms, info);

}


void UpdatePositions(queue &q, buffer<Atom> &atoms, EnvInfo env_info, SimInfo sim_info) {
	double dt = env_info.dt;
	int numAtoms = sim_info.numAtoms;
	//float volume = sim_info.boxX * sim_info.boxY * sim_info.boxZ;
	float box[] = {sim_info.boxX, sim_info.boxY, sim_info.boxZ};

	q.submit([&](handler &cgh) {
		//sycl::stream syclout(65535,255,cgh);

		auto acc = atoms.get_access<access::mode::read_write>(cgh);

		//cgh.parallel_for(range<1>(numAtoms), [=](id<1> idx) {
		cgh.parallel_for(numAtoms, [=](auto idx)
		{
			for(int i = 0; i < 3; i++)
			{
				//Velocity Verlet

				double new_accel = acc[idx[0]].force[i] / acc[idx[0]].mass;
				acc[idx[0]].position[i] += acc[idx[0]].velocity[i] * dt+0.5*acc[idx[0]].accel[i]*dt*dt;			
				acc[idx[0]].velocity[i] += 0.5*(acc[idx[0]].accel[i]+new_accel)*dt;
				acc[idx[0]].accel[i] = new_accel;

        		if(acc[idx[0]].position[i] < 0)
        			acc[idx[0]].position[i] += box[i];
        		else if(acc[idx[0]].position[i] > box[i])
        			acc[idx[0]].position[i] -= box[i];
				
			}
		});
	}).wait_and_throw();
}


int main() {
    clock_t start, finish;
    double duration, stepPerSec;
    for (auto platform : sycl::platform::get_platforms())
    {
        std::cout << "Platform: "
                  << platform.get_info<sycl::info::platform::name>()
                  << std::endl;

        for (auto device : platform.get_devices())
        {
            std::cout << "\tDevice: "
                      << device.get_info<sycl::info::device::name>()
                      << std::endl;
        }
    }

     sycl::queue q{sycl::gpu_selector_v};
	//sycl::queue q{sycl::cpu_selector_v};

	std::cout << "Running on "
				<< q.get_device().get_info<sycl::info::device::name>()
				<< "\n";

	EnvInfo env_info{.logInfo.temp = false, 
					 .logInfo.pe = false,
					 .logInfo.ke = false,
					 .logInfo.etotal = false,
					 .logInfo.evdwl = false,
					 .logInfo.ecoul = false,
					 .logInfo.epair = false,
					 .logInfo.ebond = false,
					 .logInfo.eangle = false,
					 .logInfo.edihed = false,
					 .logInfo.emol = false};
	SimInfo sim_info;

	
	std::vector<Atom> atoms;
	std::vector<Bond> bonds;
	std::vector<Angle> angles;
	std::vector<Dihedral> dihedrals;

    
	
	std::string name_info = "in.mine";
	std::string name_data = "data.lj";
	InputData inputData(name_info, name_data);

	if(!inputData.read_data(atoms, bonds, angles, dihedrals, sim_info))
	{
		std::cout << "Error occured while reading data file: " << name_data << std::endl;
		return 0;
	}

	if(!inputData.read_info(env_info, sim_info))
	{
		std::cout << "Error occured while reading information file: " << name_info << std::endl;
		return 0;
	}

	inputData.init_velocity(atoms, env_info.temp, sim_info);
	
	buffer<Atom> atoms_buf(atoms.data(), atoms.size());
	buffer<Bond> bonds_buf(bonds.data(), bonds.size());
	buffer<Angle> angles_buf(angles.data(), angles.size());
	buffer<Dihedral> dihedrals_buf(dihedrals.data(), dihedrals.size());

	Thermostat thermostat{0.0f};


	//for (int step = 0; step < env_info.numSteps; ++step) {
	std::cout << "==  Step  KE \t PE" << "\n";

	start = clock();
	for (int step = 0; step < env_info.numSteps; ++step) {
		ComputeForces(q, atoms_buf, bonds_buf, angles_buf, dihedrals_buf, env_info, sim_info);

		if(step % env_info.freq_log == 0)
		{
			writeData(atoms_buf, step);
			DisplayData(q, atoms_buf, sim_info, step);
		}

		UpdatePositions(q, atoms_buf, env_info, sim_info);

	}
	finish = clock();

	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	stepPerSec = env_info.numSteps/duration;
	std::cout << stepPerSec << " timesteps/s" << std::endl;

	//ofs.close();
	
	return 0;
}
