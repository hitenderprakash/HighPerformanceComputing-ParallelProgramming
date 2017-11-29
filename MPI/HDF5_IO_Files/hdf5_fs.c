/*
 * Program to demonstarte HDF5 Paralle I/O
 * Hitender Prakash
 * 
 */

#include <stdio.h>
#include <hdf5.h>
#include <math.h>

//particle data structure
typedef struct particle3D{
	double x;
	double y;
	double z;
	int type;
}particle_t;

#define PARTICLE_COUNT 20

int main(int argc, char **argv){
	double pi=4.0*atan(1.0);
	
	particle_t particles[PARTICLE_COUNT];
	int i;
	for(i=0;i<PARTICLE_COUNT;i++){
		particles[i].x=sin(2*pi*i/(PARTICLE_COUNT-1));
		particles[i].y=cos(2*pi*i/(PARTICLE_COUNT-1));
		particles[i].z=0.1*i;
		particles[i].type = (i>=PARTICLE_COUNT/2);
	}
	
	//create HDF5 type for data layout in memory
	int mtype = H5Tcreate(H5T_COMPOUND, sizeof(particle_t));
	
	H5Tinsert(mtype,"x coordinate",HOFFSET(particle_t,x), H5T_NATIVE_DOUBLE);
	H5Tinsert(mtype,"y coordinate",HOFFSET(particle_t,y), H5T_NATIVE_DOUBLE);
	H5Tinsert(mtype,"z coordinate",HOFFSET(particle_t,z), H5T_NATIVE_DOUBLE);
	H5Tinsert(mtype,"particle type",HOFFSET(particle_t,type), H5T_NATIVE_INT);
	
	//create HDF5 data layout in file 
	int ftype = H5Tcreate(H5T_COMPOUND,3*sizeof(double)+sizeof(int));

	H5Tinsert(ftype, "x coordinate",HOFFSET(particle_t,x),H5T_NATIVE_DOUBLE);
	H5Tinsert(ftype, "y coordinate",HOFFSET(particle_t,y),H5T_NATIVE_DOUBLE);
	H5Tinsert(ftype, "z coordinate",HOFFSET(particle_t,z),H5T_NATIVE_DOUBLE);
	H5Tinsert(ftype, "particle type",HOFFSET(particle_t,type),H5T_NATIVE_INT);
	//create data space
	hsize_t dim = PARTICLE_COUNT;
	int space =H5Screate_simple(1,&dim,NULL);
	
	//create new file with default proerties
	int fd = H5Fcreate("particles.h5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
	//create data set
	int dset = H5Dcreate(fd, "particle data", mtype,space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	//write the entire data set and close the file 
	H5Dwrite(dset,mtype,H5S_ALL, H5S_ALL, H5P_DEFAULT,particles);
	H5Fclose(fd);
	return 0;
}
