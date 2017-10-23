
//mark is an array that stores the array sizes at different levels
int mark[3] ={nx*ny*nz, ny*nz, nz};

//The following function takes in the velocity field as a 1D double array and computes the vorticity at various grid points
int arr(int dir, int i, int j, int k, int dir1, int diff);
void vort_calculator(double *velr, double ****vort)
{
	int i,j,k,dir; //Integers for running iterations
	int n[3] ={nx, ny, nz}; //array declaring inverse of step sizes in the 3d grid
	
	//Different loops for different ways to compute the different order and type of derivatives numerically
	//Forward difference
	if(typed=="f"){
		//First order forward difference
		if(order==1){
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					for(k=0;k<nz;k++){
						//Compute all three elements of vorticity in the following loop
						//A function arr is called, which returns the required position 
						//of the point required, in the 1d array, since velr is originally
						// a 4d double array
						for(dir=0;dir<nv;dir++){
							vort[dir][i][j][k] = 		
							//Since step sizes are 1/n[i], direct multiplication by n[i]
							//gives us 1/step size, which we need in finite difference formulae
							 					
							n[(dir+1)%nv]*(velr[arr(dir+2,i,j,k,dir+1,1)]
							-velr[arr(dir+2,i,j,k,dir+1,0)])
							-n[(dir+2)%nv]*(velr[arr(dir+1,i,j,k,dir+2,1)]
							-velr[arr(dir+1,i,j,k,dir+2,0)]);
						}
					}				
				}
			}
		}	
		//Second order forward difference
		if(order==2){
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					for(k=0;k<nz;k++){
						for(dir=0;dir<nv;dir++){
							vort[dir][i][j][k] = 
							(0.5)*(n[(dir+1)%nv]*(-velr[arr(dir+2,i,j,k,dir+1,2)]
							+4*velr[arr(dir+2,i,j,k,dir+1,1)]-3*velr[arr(dir+2,i,j,k,dir+1,0)])
							-n[(dir+2)%nv]*(-velr[arr(dir+1,i,j,k,dir+2,2)]
							+4*velr[arr(dir+1,i,j,k,dir+2,1)]-3*velr[arr(dir+1,i,j,k,dir+2,0)]));
						}
					}				
				}
			}
		}		
	}
	//Backward difference
	if(typed=="b"){
		if(order==1){
		//First order backward difference
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					for(k=0;k<nz;k++){
						for(dir=0;dir<nv;dir++){
							vort[dir][i][j][k] = 
							n[(dir+1)%nv]*(-velr[arr(dir+2,i,j,k,dir+1,-1)]
							+velr[arr(dir+2,i,j,k,dir+1,0)])
							-n[(dir+2)%nv]*(-velr[arr(dir+1,i,j,k,dir+2,-1)]
							+velr[arr(dir+1,i,j,k,dir+2,0)]);
						}
					}				
				}
			}
		}	
		if(order==2){
		//Second order backward difference
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					for(k=0;k<nz;k++){
						for(dir=0;dir<nv;dir++){
							vort[dir][i][j][k] = 
							(0.5)*(n[(dir+1)%nv]*(velr[arr(dir+2,i,j,k,dir+1,-2)]
							-4*velr[arr(dir+2,i,j,k,dir+1,-1)]+3*velr[arr(dir+2,i,j,k,dir+1,0)])
							-n[(dir+2)%nv]*(velr[arr(dir+1,i,j,k,dir+2,-2)]
							-4*velr[arr(dir+1,i,j,k,dir+2,-1)]+3*velr[arr(dir+1,i,j,k,dir+2,0)]));																
						}
					}				
				}
			}
		}	
	}
	//Central difference
	if(typed=="c"){
	//First order central difference
		if(order==1){
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					for(k=0;k<nz;k++){
						for(dir=0;dir<nv;dir++){
							vort[dir][i][j][k] =0.5*(n[(dir+1)%nv]*(velr[arr(dir+2,i,j,k,dir+1,1)]
							-velr[arr(dir+2,i,j,k,dir+1,-1)])-n[(dir+2)%nv]*(velr[arr(dir+1,i,j,k,dir+2,1)]
							-velr[arr(dir+1,i,j,k,dir+2,-1)]));
							//printf("%lf\n",vort[dir][i][j][k]);
							//printf("%lf\n",n[(dir+2)%nv]*(velr[arr(dir+1,i,j,k,dir+2,1)]-velr[arr(dir+1,i,j,k,dir+2,-1)]));
						}
					}				
				}
			}
		}		
		if(order==2){
		//Second order central difference
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					for(k=0;k<nz;k++){
						for(dir=0;dir<nv;dir++){
							vort[dir][i][j][k] = 
							(n[(dir+1)%nv]*(-velr[arr(dir+2,i,j,k,dir+1,2)]
							+8*velr[arr(dir+2,i,j,k,dir+1,1)]
							-8*velr[arr(dir+2,i,j,k,dir+1,-1)]
							+velr[arr(dir+2,i,j,k,dir+1,-2)])
							-n[(dir+2)%nv]*(-velr[arr(dir+1,i,j,k,dir+2,2)]
							+8*velr[arr(dir+1,i,j,k,dir+2,1)]
							-8*velr[arr(dir+1,i,j,k,dir+2,-1)]
							+velr[arr(dir+1,i,j,k,dir+2,-2)]))/12;
						}
					}				
				}
			}
		}		
	}
}


//Function that returns the position of the forward/backward/original element in the array (for any order formula)
int arr(int dir, int i, int j, int k, int dir1, int diff)
// The arguments are defined like this:
// dir refers to the the component of velocity we are differentiating
//i, j, k refer to the co-ordinates of the point of differentiation in the 3d grid
// dir1 reers to the direction of partial differentiation
// If we want to compute f(x+n*h) in finite difference, then "n" is noted as diff in this fuction.
// Diff can be used for both forward and backward difference formulae
// Diff=0 refers to the original grid point of computation.
// Periodic boundary conditions have been used. 
{
	int mark[4] ={nx*ny*nz, ny*nz, nz,1};
	int index= mark[0]*(dir%nv)+(mark[1]*i+mark[2]*j+k+diff*mark[(dir1)%nv+1])%(mark[0]);
	return index;
	//Whenever the return value or any of the index values exceed that of the highest possible value of the index, 
	//I have taken index mod(maximum possible index value), this method holds true, since we use periodic boundary conditions
}



