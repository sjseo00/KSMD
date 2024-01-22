#include "InputData.hpp"


bool InputData::read_data(std::vector<Atom>& atoms, std::vector<Bond>& bonds, std::vector<Angle>& angles, std::vector<Dihedral>& dihedrals, SimInfo &sim_info)
{
    int num_atoms = 0;
    
    std::ifstream inputFile(m_strDataFile);

    if (!inputFile.is_open()) {
        std::cout << "Error opening file: " << m_strDataFile << std::endl;
        return false;
    }

    std::string line;
    std::string boxBoundsLine;
    
    sim_info.numAtoms = 0;
    sim_info.numBonds = 0;
    sim_info.numAngles = 0;
    sim_info.numDihedrals = 0;
    sim_info.numLJs = 0;
    
    // Read header information
    ReadType readtype = NONE;
    
    while (std::getline(inputFile, line))
    {

        std::istringstream iss(line);
        if(iss.str().length() < 1) continue;
        
        if (line.find("Atoms") != std::string::npos)
        {
            readtype = ATOM;
            continue;
        }
        else if (line.find("Velocities") != std::string::npos)
        {
            readtype = VEL;
            continue;
        }        
        else if (line.find("Bonds") != std::string::npos)
        {
            readtype = BOND;
            continue;
        }
        else if (line.find("Angles") != std::string::npos)
        {
            readtype = ANGLE;
            continue;
        }
        else if (line.find("xlo xhi") != std::string::npos)
        {
            float len1, len2;
            iss >> len1 >> len2 ;
            sim_info.boxX = len2-len1;
            continue;
        } 
        else if (line.find("ylo yhi") != std::string::npos)
        {
            float len1, len2;
            iss >> len1 >> len2 ;
            sim_info.boxY = len2-len1;
            continue;
        } 
        else if (line.find("zlo zhi") != std::string::npos)
        {
            float len1, len2;
            iss >> len1 >> len2 ;
            sim_info.boxZ = len2-len1;
            continue;
        } 


       
        if ( readtype == ATOM )
        {
            
            Atom atom;
//            "atom_style	bond = atom-ID molecule-ID atom-type x y z"
            iss >> atom.id >> atom.molId >> atom.position[0] >> atom.position[1] >> atom.position[2];
            atom.mass = 1.0; /// Should be changed...
            atom.velocity[0] = 0.0;
            atom.velocity[1] = 0.0;
            atom.velocity[2] = 0.0;
            atom.accel[0] = 0.0;
            atom.accel[1] = 0.0;
            atom.accel[2] = 0.0;
            atom.force[0] = 0.0;
            atom.force[1] = 0.0;
            atom.force[2] = 0.0;
            atoms.push_back(atom);
            sim_info.numAtoms++;
        }
        else if ( readtype == BOND )
        {
            Bond bond;
            iss >> bond.id >> bond.type >> bond.atom1 >> bond.atom2;
            for (auto& parm : sim_info.bondParm) {
                if(parm.type == bond.type)
                {
                    bond.k = parm.k;
                    bond.length = parm.length;
                }
            }
            
            bonds.push_back(bond);
            sim_info.numBonds++;
        } 
        else if ( readtype == ANGLE )
        {
            // Not implemented yet!
        } 

    }

    // Print the read header information
    std::cout << "Number of atoms: " << sim_info.numAtoms << std::endl;
    std::cout << "Number of bonds: " << sim_info.numBonds << std::endl;
    std::cout << "BOX: " << sim_info.boxX << " " << sim_info.boxY << " " << sim_info.boxZ << std::endl;

    inputFile.close();
    return true;
}


bool InputData::read_info(EnvInfo &env_info, SimInfo &sim_info)
{
   
    std::ifstream inputFile(m_strInfoFile);

    if (!inputFile.is_open()) {
        std::cout << "Error opening file: " << m_strInfoFile << std::endl;
        return false;
    }

    std::string line;
    
    
    while (std::getline(inputFile, line))
    {
        std::string tag;
        std::istringstream iss(line);
        
        if(iss.str().length() < 1) continue;
       
        if (line.find("timestep") != std::string::npos)
        {
            iss >> tag >> env_info.dt;        
        }
        else if (line.find("run") != std::string::npos)
        {
            iss >> tag >> env_info.numSteps;        
        }
        else if (line.find("temperature") != std::string::npos)
        {
            iss >> tag >> env_info.temp;        
        }
        else if (line.find("freq_log") != std::string::npos)
        {
            iss >> tag >> env_info.freq_log;        
        }
        else if (line.find("freq_traj") != std::string::npos)
        {
            iss >> tag >> env_info.freq_traj;        
        }
        else if (line.find("bond_coeff") != std::string::npos)
        {
            // Should read bond_style first..
            BondParm bond;

            iss >> tag >> bond.type >> bond.k >> bond.length;
            sim_info.bondParm.push_back(bond);
        }
        else if (line.find("pair_coeff") != std::string::npos)
        {
            // Should read pair_style first..
            LJParm lj;
            iss >> tag >> lj.type1 >> lj.type2 >> lj.epsilon >> lj.sigma >> lj.cutoff;

            sim_info.ljParm.push_back(lj);
        }
        else if (line.find("log_info") != std::string::npos)
        {
            while (iss >> tag) {
                if (tag.find("temp") != std::string::npos) env_info.logInfo.temp = true;
                else if (tag.find("pe") != std::string::npos) env_info.logInfo.pe = true;
                else if (tag.find("ke") != std::string::npos) env_info.logInfo.ke = true;
                else if (tag.find("etotal") != std::string::npos) env_info.logInfo.etotal = true;
                else if (tag.find("evdwl") != std::string::npos) env_info.logInfo.evdwl = true;
                else if (tag.find("ecoul") != std::string::npos) env_info.logInfo.ecoul = true;
                else if (tag.find("epair") != std::string::npos) env_info.logInfo.epair = true;
                else if (tag.find("ebond") != std::string::npos) env_info.logInfo.ebond = true;
                else if (tag.find("eangle") != std::string::npos) env_info.logInfo.eangle = true;
                else if (tag.find("edihed") != std::string::npos) env_info.logInfo.edihed = true;
                else if (tag.find("emol") != std::string::npos) env_info.logInfo.emol = true;
            }
        }

        /* Following should be included for temperature and pressure
        info.Tdamp = 100.0f;  // Thermostat damping time scale
        info.T0 = 300.0f;  // Reference temperature
        info.tauT = 1.0f;  // Temperature relaxation time

        info.P0 = 1.0f;  // Reference pressure
        info.tauP = 1.0f;  // Pressure relaxation time
        info.alpha = 0.5f;
        */
    }
    return true;
}

struct Vector3D {
    double x, y, z;
};

double magnitude(float3 vec) {
    return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}


void setInitialTemperature_Normal(std::vector<Atom>& atoms, double desiredTemperature) {
    std::random_device rd;
    std::mt19937 gen(rd());
    double std_dev = std::sqrt(kB * desiredTemperature);
    std::normal_distribution<double> dist(0.0, std_dev);

    Vector3D totalMomentum = {0.0, 0.0, 0.0};
    
    for (Atom& atom : atoms) {
        atom.velocity[0] = dist(gen) / std::sqrt(atom.mass);
        atom.velocity[1] = dist(gen) / std::sqrt(atom.mass);
        atom.velocity[2] = dist(gen) / std::sqrt(atom.mass);
        
        totalMomentum.x += atom.mass * atom.velocity[0];
        totalMomentum.y += atom.mass * atom.velocity[1];
        totalMomentum.z += atom.mass * atom.velocity[2];
    }

    // Remove any drift (total momentum should be zero)
    size_t numAtoms = atoms.size();
    totalMomentum.x /= numAtoms;
    totalMomentum.y /= numAtoms;
    totalMomentum.z /= numAtoms;
    
    for (Atom& atom : atoms) {
        atom.velocity[0] -= totalMomentum.x / atom.mass;
        atom.velocity[1] -= totalMomentum.y / atom.mass;
        atom.velocity[2] -= totalMomentum.z / atom.mass;
    }

    // Correct for any slight deviations in temperature
    double totalKE = 0.0;
    for (const Atom& atom : atoms) {
        double v = magnitude(atom.velocity);
        totalKE += 0.5 * atom.mass * v * v;
    }
    
    double currentTemperature = (2.0 / (3.0 * kB)) * (totalKE / numAtoms);
    std::cout << "TotalKE: " << totalKE << std::endl;
    std::cout << "CurrentT: " << currentTemperature << std::endl;
    std::cout << "DesiredT: " << desiredTemperature << std::endl;
    double scalingFactor = std::sqrt(desiredTemperature / currentTemperature);
    
    for (Atom& atom : atoms) {
        atom.velocity[0] *= scalingFactor;
        atom.velocity[1] *= scalingFactor;
        atom.velocity[2] *= scalingFactor;
    }
}



void setInitialTemperature_Maxwell(std::vector<Atom>& atoms, double desiredTemperature) {
    std::random_device rd;
    std::mt19937 gen(rd());
    double m = atoms[0].mass; // Atomic mass is set to 1.0. Should be corrected.
    double variance = std::sqrt(kB * desiredTemperature / m);
    std::normal_distribution<double> distribution(0.0, variance);
    
    for (Atom& atom : atoms) {
        atom.velocity[0] = distribution(gen);
        atom.velocity[1] = distribution(gen);
        atom.velocity[2] = distribution(gen);
    }

    // Correct for any slight deviations in temperature
    double totalKE = 0.0;
    for (const Atom& atom : atoms) {
        double v = magnitude(atom.velocity);
        totalKE += 0.5 * atom.mass * v * v;
    }
    
    size_t numAtoms = atoms.size();
    double currentTemperature = (2.0 / (3.0 * kB)) * (totalKE / numAtoms);
    std::cout << "TotalKE: " << totalKE << std::endl;
    std::cout << "CurrentT: " << currentTemperature << std::endl;
}


void setInitialTemperature_Rescale(std::vector<Atom>& atoms, double desiredTemperature) {
    std::random_device rd;
    std::mt19937 gen(rd());
    size_t numAtoms = atoms.size();
    std::uniform_real_distribution<double> dist(-1*std::sqrt(desiredTemperature), std::sqrt(desiredTemperature));
    double currentKineticEnergy = 0.0;

    for (Atom& atom : atoms) {
        atom.velocity[0] = dist(gen);
        atom.velocity[1] = dist(gen);
        atom.velocity[2] = dist(gen);

        double v = magnitude(atom.velocity);
        
        currentKineticEnergy += 0.5 * atom.mass * v * v;
    }

    double desiredKineticEnergy = 1.5 * numAtoms * kB * desiredTemperature;
    double scalingFactor = std::sqrt(desiredKineticEnergy / currentKineticEnergy);

    std::cout << "Current " << currentKineticEnergy << '\n';
    std::cout << "Desired " << desiredKineticEnergy << '\n';
    std::cout << "Scaling " << scalingFactor << '\n';
    //std::cin.get();

    for (Atom& atom : atoms) {
        atom.velocity[0] *= scalingFactor;
        atom.velocity[1] *= scalingFactor;
        atom.velocity[2] *= scalingFactor;
    } 


    double totalKE = 0.0;
    for (const Atom& atom : atoms) {
        double v = magnitude(atom.velocity);
        totalKE += 0.5 * atom.mass * v * v;
    }
     
    double currentTemperature = (2.0 / (3.0 * kB)) * (totalKE / numAtoms);
    std::cout << "TotalKE: " << totalKE << std::endl;
    std::cout << "CurrentT: " << currentTemperature << std::endl;
}


void InputData::init_velocity(std::vector<Atom>& atoms, float temperature, SimInfo sim_info) {
	//setInitialTemperature(atoms, temperature);
    //setInitialTemperature_Maxwell(atoms, temperature);
    std::cout << "======= Initialize velocity =======" << std::endl;
    setInitialTemperature_Rescale(atoms, temperature);
    std::cout << "===================================" << std::endl;
}
