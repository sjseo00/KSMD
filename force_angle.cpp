#include "force.hpp"

typedef struct _ANGLE_PARM
{
	float k;
	float theta;
}AngleParm; // sj: This should be corrected as force_bond


// Angle potential kernel
void AngleKernel(queue &q, buffer<Atom> &atoms, buffer<Angle> &angles, int numAtoms, int numAngles) {
	// Compute forces from angles
	AngleParm angle_parm; // sj: This should be corrected
	
	q.submit([&](handler &cgh) {
		auto atoms_acc = atoms.get_access<access::mode::read_write>(cgh);
		auto angles_acc = angles.get_access<access::mode::read>(cgh);
		cgh.parallel_for(range<1>(numAngles), [=](id<1> idx) {
			auto angle = angles_acc[idx];
			auto atom1 = atoms_acc[angle.atom1];
			auto atom2 = atoms_acc[angle.atom2];
			auto atom3 = atoms_acc[angle.atom3];

			// Compute the angle
			float3 r21 = atom1.position - atom2.position;  // Vector from atom2 to atom1
			float3 r23 = atom3.position - atom2.position;  // Vector from atom2 to atom3
			float length_r21 = length(r21);
			float length_r23 = length(r23);

			float theta = acos(dot(r21, r23) / (length_r21 * length_r23));

			// Compute the forces
			float dV_dtheta = -2.0f * angle_parm.k * (theta - angle_parm.theta);  // Derivative of potential energy
			float3 force1 = dV_dtheta * cross(cross(r21, r23), r21) / (length_r21 * length_r21);
			float3 force3 = dV_dtheta * cross(cross(r21, r23), r23) / (length_r23 * length_r23);
			float3 force2 = -force1 - force3;

			// Add the forces to the atoms
			atoms_acc[angle.atom1].force += force1;
			atoms_acc[angle.atom2].force += force2;
			atoms_acc[angle.atom3].force += force3;
		});
	});

}
