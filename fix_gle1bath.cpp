/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Herve Ness, Christian D. Lorenz (KCL - UK) 
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_gle1bath.h"
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "pair.h"

#define MAXLINE 10000
#define MAXINBATH 1000
#define MAXDOF 500

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   Parses parameters passed to the method, allocates some memory
------------------------------------------------------------------------- */

FixGLE1B::FixGLE1B(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  int narg_min = 12;
  int nnlocal ;


  // Check to make sure we have the minimal number of inputs
  if (narg < narg_min) error->all(FLERR,"Illegal fix gle command");

  time_integrate = 1;
  restart_peratom = 1;
  MPI_Comm_rank(world,&me);

  // Parse the first set of required input arguments
  // 0 = Fix ID           (e.g., 1)
  // 1 = Group ID         (e.g., all) will be used to define bath atoms
  // 2 = gle              (name of this fix)
  // 3 = Group ID of 'system' atoms

  bathgroupL = group->find(arg[5]);
  if (bathgroupL == -1) error->all(FLERR,"Fix gle could not find fix group ID");
  bathgroupbitL = group->bitmask[bathgroupL];
  numinbathL = group->count(bathgroupL);

  sysgroup = group->find(arg[1]);
  if (sysgroup == -1) error->all(FLERR,"Fix gle could not find fix group ID");
  sysgroupbit = group->bitmask[sysgroup];

  nnlocal = group->count(sysgroup);


  //for now have temperature of bath as an input but maybe will want to use the result of a compute temperature in LAMMPS
  gle_temp_bathL = atof(arg[6]);

  // gle_delta: inputted increment with which to vibrate atoms
  gle_delta = atof(arg[3]);

  //gle_mubar: mass of virtual degrees of freedom
  gle_mubarL = atof(arg[7]);

  // gle_rseed: random number generator seed
  gle_rseed = atoi(arg[4]);
  random = new RanMars(lmp,gle_rseed + comm->me);


	printf(" \n") ;
	printf("HN read input param\n")	;
	printf("HN read input param : delta for force %f\n",gle_delta) ;
	printf("HN read input param : mubar %f\n",gle_mubarL) ;
	printf("HN read input param : rseed %i\n",gle_rseed) ;
	printf("HN read input param : nnlocal %i\n",nnlocal) ;
	printf(" \n") ;


  if (gle_rseed <= 0) error->all(FLERR,"Illegal fix gle command: must have seed > 0");
  infile = NULL;
  if (strcmp(arg[8],"doffile") == 0) {
     int n = strlen(arg[9])+1;
     infile = new char[n];
     strcpy(infile,arg[9]);
     restart_file = 1;
  }
  infile2 = NULL;
  if (strcmp(arg[10],"bathfile") == 0) {
     int n = strlen(arg[11])+1;
     infile2 = new char[n];
     strcpy(infile2,arg[11]);
     restart_file = 1;
  }
  int numlines,ncolumns,idof,idir,bathid,counter,id;
  char line[MAXLINE],*chunk;
  memory->create(gle_tauL,MAXDOF,"gle:gle_tauL");
  memory->create(gle_omegaL,MAXDOF,"gle:gle_omegaL");
  memory->create(gle_cL,MAXDOF,MAXLINE,3,"gle:gle_cL");
  memory->create(gle_mapL,atom->nmax,"gle:gle_mapL");
  memory->create(gle_mbathL,MAXINBATH,"gle:gle_mbathL");
  counter = 0;
  nvdofL = 0;

  printf("HN read doffile: %i\n",nvdofL);


  if (comm->me == 0) {
     fp1 = fopen(infile,"r");
     if (fp1 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix gle doffile %s",infile);
        error->one(FLERR,str);
     }
     while(fgets(line,MAXLINE,fp1) != NULL)
     {
       chunk = strtok(line," \t\r\n");
       while (chunk != NULL) {
          if (counter == 0) {
             nvdofL++;
             id = atoi(chunk);
             chunk = strtok(NULL," \r\n");
             counter++;
          }
          else if (counter == 1) {
             gle_tauL[id] = atof(chunk);
             chunk = strtok(NULL," \r\n");
             counter++;
          } 
          else if (counter == 2) {
             gle_omegaL[id] = atof(chunk);
             chunk = strtok(NULL," \r\n");
             counter = 0;
          }
       }
     }
     fclose(fp1);
  }

  printf("HN read doffile - DONE: %i\n",nvdofL);

  printf("HN read bathfile: %i\n",counter);


  counter = 0;
  if (comm->me == 0) {
     fp2 = fopen(infile2,"r");
     if (fp2 == NULL) {
        char str[128];
        sprintf(str,"Cannot open fix gle bathfile %s",infile2);
        error->one(FLERR,str);
     }
     id = 0;
     idir = 0;
     idof = 0;
     numlines = 0;
     ncolumns = 0;
     while(fgets(line,MAXLINE,fp2) != NULL)
     {

//  printf("HN read bathfile: line = %s%c%i\n",line," ",MAXLINE);

       chunk = strtok(line, " \t\r\n");
       while (chunk != NULL) {
//          printf("CDL counter id: %d %d %d %d \n",bathid,ncolumns,counter,id);
          if (counter == ncolumns+1) {
             counter = 1;
             numlines++;
          }
          if (counter == 0) {
             ncolumns = atoi(chunk);
             chunk = strtok(NULL, " \r\n");
             counter++;
          }
          else if (counter == 1) { 
           
             id = atoi(chunk);  
//             printf("CDL id test: %d\n",id);
             chunk = strtok(NULL, " \r\n");
             counter++;
          }
          else if (counter == 2) {
//             printf("CDL bathid test: %d\n",atoi(chunk));
             bathid = atoi(chunk);
             gle_mapL[bathid] = id;
             chunk = strtok(NULL, " \r\n");
             counter++;
          }
          else if (counter == 3) {
//             printf("CDL mbath test: %12.6f\n",atof(chunk));
             gle_mbathL[id] = atof(chunk);
             chunk = strtok(NULL," \r\n");
             counter++; 
             idof = 0;
          }
          else if (counter > 3  && counter <= ncolumns) {
             if (idir == 3) {
                idir = 0;
                idof++;
             }
//             printf("CDL c test: %12.6f\n",atof(chunk));
             gle_cL[idof][id][idir] = atof(chunk);

//HN// [30.07.14] Rescaling for ev == 23 kcal/mol , i.e 1/sqrt(23.xxx)

             gle_cL[idof][id][idir] *= 0.208239 ;

//HN//       printf("HN read bathfile: put read coef in gle_c %d %d %d %d %d %d %12.6f\n",idof,counter,ncolumns,numlines,id,idir,atof(chunk) );

             idir++;
             chunk = strtok(NULL," \r\n");
             counter++;
          } 
          else break;
       }
     }
     fclose(fp2);
  }

  printf("HN read bathfile - DONE: %i\n",counter);


  memory->create(gle_skL,MAXDOF,2,"gle:gle_skL");
  memory->create(gle_g1L,atom->nmax,3,MAXINBATH,3,"gle:gle_g1L");
  memory->create(gle_akL,MAXDOF,"gle:gle_akL");
  memory->create(gle_bkL,MAXDOF,"gle:gle_bkL");
  memory->create(gle_fp1,MAXINBATH,3,"gle:gle_fp1");
  memory->create(gle_fm1,MAXINBATH,3,"gle:gle_fm1");
  memory->create(gle_farray,MAXINBATH,3,"gle:gle_farray");
  memory->create(gle_flagrangeL,MAXINBATH,3,"gle:gle_flagrangeL");
}

/* ----------------------------------------------------------------------
   Destroys memory allocated by the method
------------------------------------------------------------------------- */
FixGLE1B::~FixGLE1B()
{
  memory->destroy(gle_tauL);
  memory->destroy(gle_omegaL);
  memory->destroy(gle_cL);
  memory->destroy(gle_skL);
  memory->destroy(gle_mapL);
  memory->destroy(gle_g1L);
  memory->destroy(gle_akL);
  memory->destroy(gle_bkL);
  memory->destroy(gle_mbathL);
  memory->destroy(gle_fp1);
  memory->destroy(gle_fm1);
  memory->destroy(gle_farray);
  memory->destroy(gle_flagrangeL);
  delete [] infile;
  delete [] infile2;
  // remove callbacks to fix, so atom class stops calling it
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);
}

/* ----------------------------------------------------------------------
   Specifies when the fix is called during the timestep
------------------------------------------------------------------------- */

int FixGLE1B::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ----------------------------------------------------------------------
   Initialize the method parameters before a run
------------------------------------------------------------------------- */

void FixGLE1B::init()
{
  double kTL;
 
  	dtv = update->dt;
  	dtf = 0.5 * update->dt * (force->ftm2v) ;
  	dtv2 = dtv * 0.5;
  	kTL = (force->boltz) * gle_temp_bathL * (force->ftm2v) ; 

//HN//	works only for units == metal
//HN//
//HN//	kT = (force->boltz) * gle_temp_bath * ftm2v : correct, so that
//HN//  units of sqrt(kT mubar) are in [g/mol][\AA][ps]^-1 as required
//HN//  for the gle_sk{k][0] vDOF, same units as momentum p_i\alpha

	printf(" \n");
	printf("HN in SUBROUTINE init() \n");
	printf("HN : bath temp = %f\n",gle_temp_bathL );
	printf("HN : local kT = %f\n",kTL );
	printf("HN : mubar = %f\n",gle_mubarL );
	printf(" \n");
	printf("HN : (force->boltz) = %f\n", (force->boltz) );
	printf("HN : (force->mvv2e) = %f\n", (force->mvv2e) );
	printf("HN : (force->ftm2v) = %f\n", (force->ftm2v) );
	printf("HN : (ftm2v * mvv2e) = %f\n", (force->ftm2v*force->mvv2e) );
	printf(" \n");
	printf("HN : dtv = %f\n", dtv );
	printf("HN : dtf = %f\n", dtf );
	printf(" \n");

	printf(" \n");
//	printf("Rescaled RANDOM variables \n");
//	printf("NO rescaling of RANDOM variables \n");
//	printf(" \n");
//	printf("UNITS Polarisation forces \n");
//	printf(" \n");

  	for (int k = 0; k < nvdofL; k++) {
      		gle_akL[k] = exp((-1.0*dtv2)/gle_tauL[k]);
      		gle_bkL[k] = sqrt( kTL * gle_mubarL * ( 1.0-(gle_akL[k]*gle_akL[k]) ) );
      		gle_skL[k][0] = 0.0;
      		gle_skL[k][1] = 0.0;
    	}
//HN//	[24.06.2014]
//HN//	following some of integration made in the fix_gld.cpp
//HN//	all coef. associated with vDOF are divided by ftm2v !
//HN//	especially the on in sqrt(kT ck)
//HN//
//HN//      gle_bk[k] = gle_bk[k] / (force->ftm2v) ;
//HN//	
//HN//	this seems to give the same resutls for Ekin as
//HN//	when working with kT = (force->boltz) * gle_temp_bath in [eV]
//HN//  and gle_bk[k] as before
//HN//  so clearly some relations between 1 / (force->ftm2v)
//HN//	and 1 / (force->mvv2e)
//HN//
//HN//	25.06.2014
//HN//  indeed ftm2v * mvv2e = 1 !
//HN//


  	if (force->pair->single_enable == 0)
    		error->all(FLERR,"Fix gle incompatible with given pair_style");
  
  	if (strstr(update->integrate_style,"respa"))
    		step_respa = ((Respa *) update->integrate)->step;

}

/* ----------------------------------------------------------------------
   First half of a timestep (V^{n} -> V^{n+1/2}; X^{n} -> X^{n+1})
------------------------------------------------------------------------- */

void FixGLE1B::initial_integrate(int vflag)
{
  double dtfm,sqmbmu, tmp1, tmp2, tmp3;
  double ftm2v = force->ftm2v;
  double fpolar[3];
  int j2, bathid;

  bigint count;

  // update v and x of atoms in group
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;

  int nlocal = atom->nlocal;

  	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

//	printf(" \n");
//	printf("HN in SUBROUTINE initial_integrate\n");
//	printf("HN : nvdof, nlocal, numinbath = %i %i %i\n",nvdof,nlocal, numinbath);
//	printf("HN : mubar = %f\n",gle_mubar);
//      printf(" \n");
//HN//  	for (int i=0; i < nlocal; i++) {
//HN//		printf("HN : atom i : xyz \n");
//HN//		printf("HN : %i  %f %f %f\n",i,x[i][0],x[i][1],x[i][2] );
//HN//		printf("HN : atom i, gle_map[i]= %i %i\n",i,gle_map[i] );
//HN//        }


  // Advance vDOF with random noise by dt/2
  	for (int k=0; k < nvdofL; k++) {
    		tmp1 = random->gaussian()  ;
    		tmp2 = random->gaussian()  ;
    		gle_skL[k][0] *= gle_akL[k];
    		gle_skL[k][0] += gle_bkL[k] * tmp1 ;
    		gle_skL[k][1] *= gle_akL[k];
    		gle_skL[k][1] += gle_bkL[k] * tmp2 ;

//HN//	printf("HN : k  = %i\n",k);
//HN//	printf("HN : gle_ak[k]  = %f\n",gle_ak[k]);
//HN//	printf("HN : gle_bk[k]  = %f\n",gle_bk[k]);
//HN//	printf("HN : random 1  = %f\n",tmp1);
//HN//	printf("HN : random 2  = %f\n",tmp2);
//HN//	printf("HN : gle_sk[k][0]  = %f\n",gle_sk[k][0]);
//HN//	printf("HN : gle_sk[k][1]  = %f\n",gle_sk[k][1]);
//HN//	printf(" \n");

//HN//	printf("HN : gle_sk[k][0]  = %f\n",gle_sk[k][0]);
  	}
//HN//	printf(" \n");


  // call the routine to calculate g12
  	g12_system_bath(gle_delta,gle_flagrangeL);

//	printf(" \n");
//	for (int j = 0; j < numinbath; j++) {
//		printf("HN : bath atom j = %i\n",j); 
//		printf("HN : bath atom j : force         = %f  %f  %f\n",gle_flagrange[j][0],gle_flagrange[j][1],gle_flagrange[j][2]); 
//	}
//	printf(" \n");

  	if (rmass) {
		printf("Inside rmass if-test \n");
		printf(" \n");

    		for (int i = 0; i < nlocal; i++) {
      			if (mask[i] & sysgroupbit) {
        			dtfm = dtf / rmass[i];

//	printf("HN : i, rmass = %i %f\n",i,rmass[i]);

	// Advance V by dt/2
        			v[i][0] += dtfm * f[i][0];
        			v[i][1] += dtfm * f[i][1];
        			v[i][2] += dtfm * f[i][2];
        // Add effects of polarisation from bath DOF
        			calc_polarisation_force(i,fpolar,gle_flagrangeL);
        			v[i][0] += dtfm * fpolar[0];
        			v[i][1] += dtfm * fpolar[1];
        			v[i][2] += dtfm * fpolar[2];
        
        // Add effects of extended DOF
				tmp1 = 0.0 ;
				tmp2 = 0.0 ;
				tmp3 = 0.0 ;

        			for (int k = 0; k <nvdofL; k++) {
          				for (int j = 0; j < numinbathL; j++) {
            					sqmbmu = sqrt(gle_mbathL[j]/gle_mubarL);
            					for (int jdir = 0; jdir <= 2; jdir++) { 
	      						tmp1 += sqmbmu * gle_g1L[i][0][j][jdir] * gle_cL[k][j][jdir] * gle_skL[k][0];
	      						tmp2 += sqmbmu * gle_g1L[i][1][j][jdir] * gle_cL[k][j][jdir] * gle_skL[k][0];
	      						tmp3 += sqmbmu * gle_g1L[i][2][j][jdir] * gle_cL[k][j][jdir] * gle_skL[k][0];
	    					}
          				}
        			}
				v[i][0] += dtfm * tmp1 ;
				v[i][1] += dtfm * tmp2 ;
				v[i][2] += dtfm * tmp3 ;

	// Advance X by dt
        			x[i][0] += dtv * v[i][0];
        			x[i][1] += dtv * v[i][1];
        			x[i][2] += dtv * v[i][2];

      			}	// END if mask
		}	// END loop on atom i

      // Advance gle_sk_2 (DOF) by dt/2
      		for (int k = 0; k < nvdofL; k++) {
        		gle_skL[k][1] += ( (-1.0*dtv2) * gle_omegaL[k] * gle_skL[k][0] );
      		}
      // Call routine to calculate g12
      		g12_system_bath(gle_delta,gle_flagrangeL);

      // Advance gle_sk_1 (DOF) by dt
      		for (int k = 0; k < nvdofL; k++) {

        		gle_skL[k][0] += 2.0*dtv2 * gle_omegaL[k] * gle_skL[k][1];

        		for (int i2 = 0; i2 < nlocal; i2++) {
          			if (mask[i2] & sysgroupbit) {
           				for (int idir = 0; idir <= 2; idir++) {
             					for (int j = 0; j < numinbathL; j++) {

               						sqmbmu = sqrt( gle_mbathL[j] * gle_mubarL );

               						for (int jdir = 0; jdir <= 2; jdir++) {
                 			gle_skL[k][0] += -2.0*dtv2 * sqmbmu * gle_g1L[i2][idir][j][jdir] * gle_cL[k][j][jdir] * v[i2][idir];
               						}
             					}
           				}
          			}
        		}
      		}      // END loop : Advance gle_sk_1 (DOF) by dt
    	} 	// END if rmass
	else
	{
//		printf("Inside other rmass option if-test CORRECT \n");
//		printf(" \n");

		for (int i = 0; i < nlocal; i++) {
			if (mask[i] & sysgroupbit) {
//
//	FOR atom i in the system region ONLY
//
        			dtfm = dtf / mass[type[i]];

//				printf("HN : i_sys, tag[i], mask[i], mass = %i  %i  %i  %f\n",i,tag[i],mask[i],mass[type[i]]);

	// Advance V by dt/2
        			v[i][0] += dtfm * f[i][0];
        			v[i][1] += dtfm * f[i][1];
        			v[i][2] += dtfm * f[i][2];

//HN//				printf("HN : dtf   = %f\n",dtf);
//HN//				printf("HN : dtfm  = %f\n",dtfm);
//				printf("HN : fx[i_sys] = %f\n",f[i][0]);
//				printf("HN : fy[i_sys] = %f\n",f[i][1]);
//				printf("HN : fz[i_sys] = %f\n",f[i][2]);
//				printf(" \n");

        // Add effects of polarisation from bath DOF
        			calc_polarisation_force(i,fpolar,gle_flagrangeL);
//        			calc_polarisation_force_diag(i,fpolar,gle_flagrangeL);
        			v[i][0] += dtfm * fpolar[0];
        			v[i][1] += dtfm * fpolar[1];
        			v[i][2] += dtfm * fpolar[2];
        
//HN//	printf("HN : dtfm  = %f\n",dtfm);
//				printf("HN : fpol_x[i_sys] = %f\n",fpolar[0]);
//				printf("HN : fpol_y[i_sys] = %f\n",fpolar[1]);
//				printf("HN : fpol_z[i_sys] = %f\n",fpolar[2]);
//				printf(" \n");

        // Add effects of extended DOF
				tmp1 = 0.0 ;
				tmp2 = 0.0 ;
				tmp3 = 0.0 ;

				for (int k = 0; k <nvdofL; k++) {
          				for (int j = 0; j < numinbathL; j++) {
            					sqmbmu = sqrt( gle_mbathL[j] / gle_mubarL );

//HN//	printf("HN : l_bath, mass = %i %f\n",i,gle_mbath[j]);

            					for (int jdir = 0; jdir <= 2; jdir++) { 
	      						tmp1 += dtfm * sqmbmu * gle_g1L[i][0][j][jdir] * gle_cL[k][j][jdir] * gle_skL[k][0];
	      						tmp2 += dtfm * sqmbmu * gle_g1L[i][1][j][jdir] * gle_cL[k][j][jdir] * gle_skL[k][0];
	      						tmp3 += dtfm * sqmbmu * gle_g1L[i][2][j][jdir] * gle_cL[k][j][jdir] * gle_skL[k][0];

//HN//	printf("HN : dtfm * sqmbmu  = %f\n",dtfm * sqmbmu);
//HN//	printf("HN : i, j, jdir  = %i %i %i\n",i,j,jdir);
//HN//	printf("HN : gle_g1L[i][0][j][jdir]  = %f\n",gle_g1L[i][0][j][jdir]);
//HN//	printf("HN : gle_g1L[i][1][j][jdir]  = %f\n",gle_g1L[i][1][j][jdir]);
//HN//	printf("HN : gle_g1L[i][2][j][jdir]  = %f\n",gle_g1L[i][2][j][jdir]);
//HN//	printf("HN : k, j, jdir  = %i %i %i\n",i,j,jdir);
//HN//	printf("HN : gle_cL[j][jdir]  = %f\n",gle_cL[k][j][jdir]);
//HN//	printf("HN : k  = %i\n",k);
//HN//	printf("HN : gle_skL[k][0]  = %f\n",gle_skL[k][0]);

	    					}
          				}
        			}
//				printf("HN : effects of vDOF xdir  = %f\n",tmp1);
//				printf("HN : effects of vDOF ydir  = %f\n",tmp2);
//				printf("HN : effects of vDOF zdir  = %f\n",tmp3);
//				printf(" \n");

				v[i][0] += tmp1 ;
				v[i][1] += tmp2 ;
				v[i][2] += tmp3 ;

	// Advance X by dt
        			x[i][0] += dtv * v[i][0];
        			x[i][1] += dtv * v[i][1];
        			x[i][2] += dtv * v[i][2];

//HN//	printf("HN : dtv  = %f\n",dtv);
//HN//	printf("HN : x[i_sys],  = %i %f\n",i,x[i][0]);
//HN//	printf("HN : y[i_sys],  = %i %f\n",i,x[i][1]);
//HN//	printf("HN : y[i_sys],  = %i %f\n",i,x[i][2]);
//HN//	printf(" \n");
			}	// END if mask
		} 	// END loop on atom i
	}	// END if-else 
	
      // Advance gle_skL_2 (DOF) by dt/2
	for (int k = 0; k < nvdofL; k++) {
//HN// [30.07.14] Rescaling by (2\pi) coming from \omega = 2PI \nu
        	gle_skL[k][1] += ( (-1.0*dtv2) * gle_omegaL[k] * gle_skL[k][0] );

      	}

      // Call routine to calculate g12
      	g12_system_bath(gle_delta,gle_flagrangeL);

      // Advance gle_skL_1 (DOF) by dt
      	for (int k = 0; k < nvdofL; k++) {

//HN// [30.07.14] Rescaling by (2\pi) coming from \omega = 2PI \nu
	       	gle_skL[k][0] += dtv * gle_omegaL[k] * gle_skL[k][1];


      // Advance the rest of gle_sk_1 (DOF) by dt
        	tmp1 = 0.0 ;
		for (int i2 = 0; i2 < nlocal; i2++) {
          		if (mask[i2] & sysgroupbit) {
           			for (int idir = 0; idir <= 2; idir++) {
             				for (int j = 0; j < numinbathL; j++) {
                				sqmbmu = sqrt( gle_mbathL[j] * gle_mubarL );
               					for (int jdir = 0; jdir <= 2; jdir++) {

                 					tmp1 += sqmbmu * gle_g1L[i2][idir][j][jdir] * gle_cL[k][j][jdir] * v[i2][idir];
               					}
             				}	// END loop bath atom j
           			}	// END idir atom i2
          		}	// END test atom i2 in sys grp
        	}	// END loop atom i2
//HN//	[27.06.14]
//HN//  See page 174, 175 logbook
//HN//	Looks like I need the same conversion factor in step time
//HN//	as for the propagation of the velocity, presence of mass factor and co
//HN//  need to use 2* dtf = 2* (0.5 * dt * ftm2v), need full time step.
//HN//
		gle_skL[k][0] += -1.0 * dtv * ftm2v * tmp1 ;

//HN//		printf("HN : vDOF k  = %i \n",k);
//HN//		printf("HN : gle_skL[k][0]  = %f\n",gle_skL[k][0]);	

      	}	// END lool vDOF k

//HN//    	printf(" \n");

}	// END subroutine


/* ----------------------------------------------------------------------
   Second half of a timestep (V^{n+1/2} -> V^{n+1})
------------------------------------------------------------------------- */

void FixGLE1B::final_integrate()
{
  double dtfm,sqmbmu,tmp1, tmp2, tmp3 ;
  double fpolar[3];
  // update v of atoms in group

  double ftm2v = force->ftm2v;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int *tag = atom->tag;

  int nlocal = atom->nlocal;

  	if (igroup == atom->firstgroup) nlocal = atom->nfirst;

//	printf(" \n");
//	printf("HN in SUBROUTINE final_integrate \n");
//HN//	printf("HN : nvdof, numinbath = %i %i\n",nvdof,numinbath);
//HN//	printf("HN : mubar = %f\n",gle_mubar);
	printf(" \n");


	if (rmass) {
    		for (int i = 0; i < nlocal; i++) {
      			if (mask[i] & sysgroupbit) {
        			dtfm = dtf / rmass[i];
        			v[i][0] += dtfm * f[i][0];
        			v[i][1] += dtfm * f[i][1];
        			v[i][2] += dtfm * f[i][2];
        // Add effects of polarisation from bath DOF
        			calc_polarisation_force(i,fpolar,gle_flagrangeL);
        			v[i][0] += dtfm * fpolar[0];
        			v[i][1] += dtfm * fpolar[1];
        			v[i][2] += dtfm * fpolar[2];
        // Add effects of extended DOF 
        			for (int k = 0; k < nvdofL; k++) {
          				for (int j = 0; j < numinbathL; j++) {
            					sqmbmu = sqrt(gle_mbathL[j]/gle_mubarL);
            					for (int jdir = 0; jdir <= 2; jdir++) {
              						v[i][0] += dtfm*sqmbmu*gle_g1L[i][0][j][jdir]*gle_cL[k][j][jdir]*gle_skL[k][0];
              						v[i][1] += dtfm*sqmbmu*gle_g1L[i][1][j][jdir]*gle_cL[k][j][jdir]*gle_skL[k][0];
              						v[i][2] += dtfm*sqmbmu*gle_g1L[i][2][j][jdir]*gle_cL[k][j][jdir]*gle_skL[k][0];
            					}
          				}
        			}
      			}
		}
	} 
	else 
	{
    		for (int i = 0; i < nlocal; i++) {
      			if (mask[i] & sysgroupbit) {

        			dtfm = dtf / mass[type[i]];

//				printf("HN : i_sys, tag[i], mask[i], mass = %i  %i  %i  %f\n",i,tag[i],mask[i],mass[type[i]]);

	// Advance V by dt/2
        			v[i][0] += dtfm * f[i][0];
        			v[i][1] += dtfm * f[i][1];
        			v[i][2] += dtfm * f[i][2];

//				printf("HN : fx[i_sys] = %f\n",f[i][0]);
//				printf("HN : fy[i_sys] = %f\n",f[i][1]);
//				printf("HN : fz[i_sys] = %f\n",f[i][2]);
//				printf(" \n");

        // Add effects of polarisation from bath DOF
        			calc_polarisation_force(i,fpolar,gle_flagrangeL);
//        			calc_polarisation_force_diag(i,fpolar,gle_flagrangeL);
       	 			v[i][0] += dtfm * fpolar[0];
        			v[i][1] += dtfm * fpolar[1];
        			v[i][2] += dtfm * fpolar[2];

//				printf("HN : fpol_x[i_sys] = %f\n",fpolar[0]);
//				printf("HN : fpol_y[i_sys] = %f\n",fpolar[1]);
//				printf("HN : fpol_z[i_sys] = %f\n",fpolar[2]);
//				printf(" \n");

        // Add effects of extended DOF 
				tmp1 = 0.0 ;
				tmp2 = 0.0 ;
				tmp3 = 0.0 ;

        			for (int k = 0; k < nvdofL; k++) {
          				for (int j = 0; j < numinbathL; j++) {
            					sqmbmu = sqrt(gle_mbathL[j]/gle_mubarL);
            					for (int jdir = 0; jdir <= 2; jdir++) {
              						tmp1 += dtfm*sqmbmu*gle_g1L[i][0][j][jdir]*gle_cL[k][j][jdir]*gle_skL[k][0];
              						tmp2 += dtfm*sqmbmu*gle_g1L[i][1][j][jdir]*gle_cL[k][j][jdir]*gle_skL[k][0];
              						tmp3 += dtfm*sqmbmu*gle_g1L[i][2][j][jdir]*gle_cL[k][j][jdir]*gle_skL[k][0];
            					}
          				}
        			}
        			v[i][0] += tmp1 ;
        			v[i][1] += tmp2 ;
        			v[i][2] += tmp3 ;

//				printf("HN : effects of vDOF xdir  = %f\n",tmp1);
//				printf("HN : effects of vDOF ydir  = %f\n",tmp2);
//				printf("HN : effects of vDOF zdir  = %f\n",tmp3);
//				printf(" \n");

      			}	// END if mask
		}	// END loop atom i
	}	// END if-else-endif

  // Advance extended DOF-2 (gle_skL[][1]) by dt/2
  	for (int k = 0; k < nvdofL; k++) {
//HN// [30.07.14] Rescaling by (2\pi) coming from \omega = 2PI \nu
    		gle_skL[k][1] += (-1.0*dtv2) * gle_omegaL[k] * gle_skL[k][0];
  	}

  // Advance vDOF with random noise by dt/2
  	for (int k = 0; k < nvdofL; k++) {
    		gle_skL[k][0] *= gle_akL[k];
    		gle_skL[k][0] += gle_bkL[k] * (random->gaussian())  ;
    		gle_skL[k][1] *= gle_akL[k];
    		gle_skL[k][1] += gle_bkL[k] * (random->gaussian())  ;

//HN//		printf("HN : vDOF k  = %i \n",k);
//HN//		printf("HN : gle_skL[k][0]  = %f\n",gle_skL[k][0]);	
//HN//		printf("HN : gle_skL[k][1]  = %f\n",gle_skL[k][1]);	
  	}
//HN//	printf(" \n");

}	// END subroutine

/* ----------------------------------------------------------------------
   Calculation of the g_{i\alpha,l\gamma}(r) matrix elements between
   system DOF (i\alpha) and bath DOF (l\gamma), knowing that (r) is
   the set of all positions of the system DOFs
-------------------------------------------------------------------------
   Subroutine should return all matrix elements	g12[i][\alpha][j][\gamma]
	i : system atom 
	j : bath atom
	\alpha, \gamma : coord x,y,z [0,1,2] 
------------------------------------------------------------------------- */
void FixGLE1B::g12_system_bath(double gle_delta, double **gle_flagrange)
{
// local variables
   double delx,dely,delz,rsq, tmp;
   double fpair,gle_deltax2;
   int itype;
   double **x = atom->x;
   int *type = atom->type;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   if (igroup == atom->firstgroup) nlocal = atom->nfirst;
   
// fp1 = force on bath after moving atom in system in positive direction
// fm1 = force on bath after moving atom system in negative direction
// double fp1[0  to Nb_BATH_ATOM][0 to 2]
// double fm1[0  to Nb_BATH_ATOM][0 to 2]


   	gle_deltax2 = 2.0 * gle_delta ;

// Start loop on System atom i
   	for (int i = 0; i < nlocal ; i++) {
     		if (mask[i] & sysgroupbit) {

      			itype = type[i];        

			for (int idir = 0; idir <= 2; idir++) {

//    direction idir = x,y,z with (+gle_delta) for SYSTEM ATOM i 
      				x[i][idir] += gle_delta ;
//				calc_force_on_bath(gle_fp1) ;
				calc_force_on_bath_gen(gle_fp1) ;

//    direction idir = x,y,z with (-gle_delta) for SYSTEM ATOM i 
      				x[i][idir] -= (gle_deltax2) ;
//				calc_force_on_bath(gle_fm1) ;
				calc_force_on_bath_gen(gle_fm1) ;

//    Return to intial position
      				x[i][idir] += gle_delta ;

      				for (int j = 0; j < numinbathL ; j++) {

					tmp = gle_mbathL[j] * gle_deltax2 ;
					gle_g1L[i][idir][j][0] = ( gle_fp1[j][0] - gle_fm1[j][0] ) / tmp ;
					gle_g1L[i][idir][j][1] = ( gle_fp1[j][1] - gle_fm1[j][1] ) / tmp ;
					gle_g1L[i][idir][j][2] = ( gle_fp1[j][2] - gle_fm1[j][2] ) / tmp ;
      				}	
//    End loop on idir
			}	
//    End loop on System atom i
     		}
     	}

//    Need to know "force" (interac Langrangian) for the determination of the
//    polarisation force

//      calc_pair_force(i,itype,gle_flagrange);
//	calc_force_on_bath(gle_flagrange) ;
	calc_force_on_bath_gen(gle_flagrange) ;

      	for (int j = 0; j < numinbathL; j++) {
		gle_flagrange[j][0] = gle_flagrange[j][0] / gle_mbathL[j] ;
		gle_flagrange[j][1] = gle_flagrange[j][1] / gle_mbathL[j] ;
		gle_flagrange[j][2] = gle_flagrange[j][2] / gle_mbathL[j] ;
      	}	

/* ---------------------------------------------------------------------- */
  return;
}


/* ----------------------------------------------------------------------
   Calculate force f_b on bath atoms due to system region
   Should work for any type of interatomic potential
------------------------------------------------------------------------- */

void FixGLE1B::calc_force_on_bath_gen(double **gle_farray)
{
   int j,j2, bathid;
   double delx,dely,delz,rsq;
   double fpair = 0.0;
   double factor_coul = 1.0;
   double factor_lj = 1.0;
   double tmp = 0.0;
   double ffx, ffy, ffz ;
   pair = force->pair;
   cutsq = force->pair->cutsq;

   double **x = atom->x;
   double **f = atom->f;   
   int *type = atom->type;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   int *tag = atom->tag;

    if (igroup == atom->firstgroup) nlocal = atom->nfirst;
   
//HN// 	[14.10.14]
//HN//	Test on call of                 
//HN//  pair->compute(int eflag, int vflag)
//HN// 	should recalc all the forces, whatever the kind of interactomic potential?
//HN//

	for (int j = 0; j < nlocal ; j++) {
		f[j][0] = 0 ;
		f[j][1] = 0 ;
		f[j][2] = 0 ;
      	}
	pair->compute(0,1) ;

// Now array f[i][.] should be updated with proper forces, for corresponding atomic positions

    
// Transpose force components of f[i][.] into components of gle_farray[j2][.]
// corresponding to the entries of the bath reduced region DOF of cbk coeff

	for (int j = 0; j < nlocal; j++) {
     		if (mask[j]&bathgroupbitL) {

			bathid = tag[j] ; 
			j2 = gle_mapL[bathid];

      			gle_farray[j2][0] = f[j][0] ;
      			gle_farray[j2][1] = f[j][1] ;
      			gle_farray[j2][2] = f[j][2] ;
     		}
   	}
  	return;
}


/* ----------------------------------------------------------------------
   Calculate force f_b on bath atoms due to system region
   works only for pair potential of LJ type
------------------------------------------------------------------------- */

void FixGLE1B::calc_force_on_bath(double **gle_farray)
{
   int itype,jtype,j,j2, bathid;
   double delx,dely,delz,rsq;
   double fpair = 0.0;
   double factor_coul = 1.0;
   double factor_lj = 1.0;
   double tmp = 0.0;
   double ffx, ffy, ffz ;
   pair = force->pair;
   cutsq = force->pair->cutsq;

   double **x = atom->x;
   int *type = atom->type;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   int *tag = atom->tag;

    if (igroup == atom->firstgroup) nlocal = atom->nfirst;
   
//HN//	printf(" \n") ;
//HN//	printf("HN in subroutine calc_force_on_bath\n") ;

	for (int j = 0; j < nlocal; j++) {
     		if (mask[j]&bathgroupbitL) {

//HN// Here j is an atom of the bath-reduced region
      			jtype = type[j];
//HN//
//HN// [01.07.14]
//HN// What we actually need here, is gle_map[ ID of atom j] = j2
//HN// and that will fit with the entries of input file of cbk coefficients!
//HN//
			bathid = tag[j] ; 
			j2 = gle_mapL[bathid];

			ffx = 0.0 ;
			ffy = 0.0 ;
			ffz = 0.0 ;

//HN//			printf(" \n") ;
//HN//			printf("Atom j, tag[j]                    = %i  %i\n",j,bathid) ;
//HN//			printf("mask[j], type[j], gle_map[tag[j]] = %i  %i  %i\n",mask[j],type[j],j2) ;
//HN//			printf(" \n") ;

			for (int i = 0; i < nlocal; i++) {
     				if (mask[i] & sysgroupbit) {
//HN// Here i is an atom of the system region
					itype = type[i]; 

//HN//					printf("Atom i , tag[i] , mask[i] = %i  %i  %i\n",i,tag[i],mask[i]) ;

      					delx = x[i][0] - x[j][0];
      					dely = x[i][1] - x[j][1];
      					delz = x[i][2] - x[j][2];
      					rsq = delx*delx + dely*dely + delz*delz;
      					if (rsq < cutsq[itype][jtype])
       						tmp = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);	

					ffx += delx * fpair ;
					ffy += dely * fpair ;
					ffz += delz * fpair ;

//HN//					printf("fpair, delx, dely, delz = %f  %f  %f  %f\n",fpair,delx,dely,delz) ;
 				}
			}

//HN// 
//HN//      printf(" \n") ;
//HN//      printf("HN : in FixGLE::calc_pair_force\n") ;
//HN//      printf("HN : atom, gle_map  = %i %i\n",j,j2);

//HN//      CHECK SIGN : WE NEED THE FORCE ON ATOM J [Put negative sign, see my notes on pages 182, 189]
//HN//
      			gle_farray[j2][0] = -1.0 * ffx ;
      			gle_farray[j2][1] = -1.0 * ffy ;
      			gle_farray[j2][2] = -1.0 * ffz ;
     		}
   	}
  	return;
}


/* ----------------------------------------------------------------------
   Calculate polarisation force on atom i
------------------------------------------------------------------------- */

void FixGLE1B::calc_polarisation_force(int i,double *fpolar,double **gle_flagrange)
{
   double mass2, sumcoef;
   double ftm2v = force->ftm2v;

	for (int idir = 0; idir <= 2 ; idir++) {
     		sumcoef = 0.;
     		for (int b = 0; b < numinbathL ; b++) {
       			for (int bdir = 0; bdir <= 2 ; bdir++) {
	 			for (int bp = 0; bp < numinbathL ; bp++) {
	   				mass2 = sqrt(gle_mbathL[b]*gle_mbathL[bp]);				
           					for (int bpdir = 0; bpdir <= 2 ; bpdir++) {
	     						for (int k = 1; k <= nvdofL; k++) {
               	sumcoef += gle_g1L[i][idir][b][bdir]*gle_cL[k][b][bdir]*gle_cL[k][bp][bpdir]* gle_flagrange[bp][bpdir]*mass2;			
	     						}
	   					}
         			}
       			}
     		}
//     		fpolar[idir] = sumcoef * ftm2v ;
     		fpolar[idir] = sumcoef ;
   	}

  return;
}

/* ----------------------------------------------------------------------
   Calculate polarisation force on atom i, using only diag elements cbk
------------------------------------------------------------------------- */

void FixGLE1B::calc_polarisation_force_diag(int i,double *fpolar,double **gle_flagrange)
{
   double mass2, sumcoef;
   double ftm2v = force->ftm2v;

	for (int idir = 0; idir <= 2 ; idir++) {
     		sumcoef = 0.;
     		for (int b = 0; b < numinbathL ; b++) {
       			for (int bdir = 0; bdir <= 2 ; bdir++) {
	   			mass2 = gle_mbathL[b];				
	     			for (int k = 1; k <= nvdofL; k++) {
               	sumcoef += gle_g1L[i][idir][b][bdir]*gle_cL[k][b][bdir]*gle_cL[k][b][bdir]* gle_flagrange[b][bdir]*mass2;			
	     			}
	   		}
     		}
     		fpolar[idir] = sumcoef ;
   	}

  return;
}

/* ---------------------------------------------------------------------- */

void FixGLE1B::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * (force->ftm2v);

  // innermost level - GLE update of v and x
  // all other levels - GLE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixGLE1B::final_integrate_respa(int ilevel, int iloop)
{
  dtf = 0.5 * step_respa[ilevel] * (force->ftm2v);
  final_integrate();
}

/* ----------------------------------------------------------------------
   Called when a change to the timestep is requested mid-run
------------------------------------------------------------------------- */

void FixGLE1B::reset_dt()
{
  // set the time integration constants
  dtv = update->dt;
  dtv2 = dtv*0.5;
  dtf = 0.5 * update->dt * (force->ftm2v);
}

/* ----------------------------------------------------------------------
   Calculate force (on atom i) between pairs of atoms i,j
------------------------------------------------------------------------- */

void FixGLE1B::calc_pair_force(int i,int itype,double **gle_farray)
{
   int jtype,j,j2;
   double delx,dely,delz,rsq;
   double fpair = 0.0;
   double factor_coul = 1.0;
   double factor_lj = 1.0;
   double tmp = 0.0;
   pair = force->pair;
   cutsq = force->pair->cutsq;

   double **x = atom->x;
   int *type = atom->type;
   int *mask = atom->mask;
   int nlocal = atom->nlocal;
   if (igroup == atom->firstgroup) nlocal = atom->nfirst;
   
//HN// Here i is an atom of the system region
//HN// Here j is an atom of the bath-reduced region
	fpair = 0.0;
	for (int j = 0; j < nlocal; j++) {
     		if (mask[j]&bathgroupbitL) {
      			jtype = type[j];
      			delx = x[i][0] - x[j][0];
      			dely = x[i][1] - x[j][1];
      			delz = x[i][2] - x[j][2];
      			rsq = delx*delx + dely*dely + delz*delz;
      			if (rsq < cutsq[itype][jtype])
       				tmp = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);	

      			j2 = gle_mapL[j];

//HN// 
//HN//      printf(" \n") ;
//HN//      printf("HN : in FixGLE::calc_pair_force\n") ;
//HN//      printf("HN : atom, gle_map  = %i %i\n",j,j2);

//    CHECK SIGN : WE NEED THE FORCE ON ATOM J [Put negative 1 for the time being but revisit]
      			gle_farray[j2][0] += -1.0*delx*fpair;
      			gle_farray[j2][1] += -1.0*dely*fpair;
      			gle_farray[j2][2] += -1.0*delz*fpair;
     		}
   	}
  	return;
}


