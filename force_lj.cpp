#include "force.hpp"



// This is a simple LJ interaction kernel
void LennardJonesKernel(queue &q, buffer<Atom> &atoms, EnvInfo env_info, SimInfo sim_info) {

  int numAtoms = sim_info.numAtoms;
  // This should be corrected!!!!
  double sigma = sim_info.ljParm[0].sigma;
  double epsilon = sim_info.ljParm[0].epsilon;
  float cutoff = sim_info.ljParm[0].cutoff;
  float box_x = sim_info.boxX;  
  float box_y = sim_info.boxY;
  float box_z = sim_info.boxZ;
  
  q.submit([&](handler &cgh) {
    auto atoms_acc = atoms.get_access<access::mode::read_write>(cgh);

    cgh.parallel_for(numAtoms, [=](auto idx)
    {
      for (int j = idx[0]+1; j < numAtoms; j++)
      {
        double rx = atoms_acc[idx[0]].position[0] - atoms_acc[j].position[0];
        double ry = atoms_acc[idx[0]].position[1] - atoms_acc[j].position[1];
        double rz = atoms_acc[idx[0]].position[2] - atoms_acc[j].position[2];
            
        // To apply PBC
//    	syclout << atoms_acc[j].position[0] << " " << atoms_acc[idx[0]].position[0] << " rx = " << rx << cl::sycl::endl;

        double dx = (sycl::abs(rx) > box_x/2.0)? box_x-rx : rx;
        double dy = (sycl::abs(ry) > box_y/2.0)? box_y-ry : ry;
        double dz = (sycl::abs(rz) > box_z/2.0)? box_z-rz : rz;

		double rsq = dx*dx + dy*dy + dz*dz;
        double r = sycl::sqrt(rsq);

        double energy = 0.0f;
        double force = 0.0f;
        double cutoffsq = cutoff*cutoff;
        
        if (rsq < cutoffsq)
        {
          double sr2 = 1.0/rsq;
          double sr6 = sr2 * sr2 * sr2 * sigma;
          double force = 48.0 * sr6 * (sr6 - 0.5) * sr2 * epsilon;
          energy = (4.0 * sr6 * (sr6-1.0)) * epsilon;

          atoms_acc[idx[0]].force[0] += force*dx;
          atoms_acc[idx[0]].force[1] += force*dy;
          atoms_acc[idx[0]].force[2] += force*dz;
          atoms_acc[j].force[0] -= force*dx;
          atoms_acc[j].force[1] -= force*dy;
          atoms_acc[j].force[2] -= force*dz;
        }
        atoms_acc[idx[0]].pe += energy;
      }
      
    }); 
    
  }).wait();
  
}


