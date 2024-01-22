#include "force.hpp"

// Bond potential kernel
void BondKernel(queue &q, buffer<Atom> &atoms, buffer<Bond> &bonds, EnvInfo env_info, SimInfo sim_info) {
    // Compute forces from bonds
    int numBonds = sim_info.numBonds;
  
    
    q.submit([&](handler &cgh) {
        auto atoms_acc = atoms.get_access<access::mode::read_write>(cgh);
        auto bonds_acc = bonds.get_access<access::mode::read>(cgh);
        
        cgh.parallel_for(range<1>(numBonds), [=](id<1> idx) {
            auto bond = bonds_acc[idx];
            auto atom1 = atoms_acc[bond.atom1];
            auto atom2 = atoms_acc[bond.atom2];
            float3 delta = atom1.position - atom2.position;
            float distance = length(delta);
            float3 force = bond.k * (bond.length - distance) * normalize(delta);
            atoms_acc[bond.atom1].force += -force;
            atoms_acc[bond.atom2].force += force;
        });
    });
    q.wait();
    
}
