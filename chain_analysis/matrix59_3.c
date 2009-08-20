/*
  Dump analysis for LAMMPS 2004
  Won-Ki Roh, Dec 2004/Jan 2005
  Modifications EL
  * revised Feb 1 2005 (add function to calculate bond length) *
  * revised Feb 4 2005 (add function to measure no. of adsorbed monomers)*v.54
  * revised Mar 11 2005 (add function to measure no. of adsorbed chains)*v.56
  * revised Mar 11-12 2005 (add function to measure no. of loop,train,tail's monomers)*v.57
  * revised Mar 22-24 2005 (add function to measure length & no. of loop,train,tail)*v.58
  * revised Mar 26    2005 (add function to calculate error bar of msd)*v.59
  * revised Apr 8     2005 (checking the dump file's disorder of configuration)*v.59_2
  * revised Apr 13    2005 (fix error bar of msd)*v.59_3
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "alloc2d.h"
#include "alloc3d.h"

#define TRUE  1
#define FALSE 0


#define KEYWORDS {"dumpname", "timesteps", "maxdeltat", "chainlength", "no_chains", "ad_condition_z", "compute_msd", "outputfile", "msdout", "adsorb"}
#define N_KEYWORDS 10


FILE *fp;
FILE *out1;
FILE *out2;
FILE *out3;

int main(int argc, char *argv[])
{
	FILE *cfg;
	char keyword[1024], value[1024], dumpname[1024], command[1024], outputfile[1024], msdout[1024], adsorb[1024];
	char *keywordlist[] = KEYWORDS;
	double **position_array, ***cm_array;
	int timestep, i, lineno, readstatus, compressed, length;
	int notimesteps, maxdelta, chainlength, nochains, msdon; 
	int total_nomon_ad, total_nochain_ad, no_train_mon, no_tail_mon, no_loop_mon, no_nonadsorb_mon, no_train, no_tail, no_loop, no_nonadsorb;
	double avetrain_mon, avetail_mon, aveloop_mon, avenonadsorb_mon;
	double ad_condition;
	double aveRg, aveaveRg, aveRgsum, aveRe, aveavel, coordaveRg[3], coordaveRe[3], coordaveavel[3];
	
	
	int readconfig(double **positions, int nchainis, int chainlength);
	void printconfig(double **positions, int atoms);
	double calcgyration(double ***cm_array, double **positions, int chainlength, int nchains, int timestep, double coordaveRg[3]);
	void calccentermass(double ***cm_array, double **positions, int chainlength, int nchains, int timestep);
	double calcmonomerdistance(double **positions, int chainlength, int nchains, double aveavel, double coordaveavel[3]);
	double calcendlength(double **positions, int chainlength, int nchains, double aveRe, double coordaveRe[3]);
	void calcmsd(double ***cm_array, int notimesteps, int nchains, int maxdelta);
	int adsorptioncondition(double **positions, int chainlength, int nchains, int total_nomon_ad, double ad_condition, int *total_nochain_ad);
	void looptraintail(double **positions, int chainlength, int nchains, double ad_condition, int *no_train_mon, int *no_tail_mon, int *no_loop_mon, int *no_nonadsorb_mon, double *avetrain_mon, double *avetail_mon, double *aveloop_mon, double *avenonadsorb_mon, int *no_train, int *no_tail, int *no_loop, int *no_nonadsorb);
	
	if (argc<2)
	{
		fprintf(stderr, "Usage: %s configurationfile\n", argv[0]);
		exit(0);
	}
	
	cfg=fopen(argv[1], "r");

	if (cfg == NULL)
	{
		fprintf(stderr, "Cannot open configuration file %s\n", argv[1]);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(stderr, "Reading configuration file %s\n", argv[1]);
	}

	readstatus = fscanf(cfg, "%s = %s\n", keyword, value);
	lineno = 1;

	while (readstatus != EOF)
	{
	    i=0;
	    while (i<N_KEYWORDS)
	    {
		if (strcmp(keyword, keywordlist[i])==0) break;
		i++;
	    }
	    if (i == N_KEYWORDS)
	    {
		fprintf(stderr, "Parse error in line %d\n", lineno);
		exit(EXIT_FAILURE);
	    }
	    switch(i)
	    {
		case(0):
		    strcpy(dumpname, value);
		    printf("\n");
		    printf("Dumpname: %s\n", dumpname);
		    break;
		case(1):
		    notimesteps=atoi(value);
		    printf("No. of timesteps: %d\n", notimesteps);
		    break;
		case(2):
		    maxdelta=atoi(value);
		    printf("Max. delta t: %d\n", maxdelta);
		    break;
		case(3):
		    chainlength=atoi(value);
		    printf("Chainlength: %d\n", chainlength);
		    break;
		case(4):
		    nochains=atoi(value);
		    printf("No. of chains: %d\n", nochains);
		    break;
		case(5):
		    ad_condition=atof(value);
		    printf("Adsorb condition z: %10.5f\n", ad_condition);
		    break;
		case(6):
		    msdon=atoi(value);
		    if(msdon==1)
			printf("Compute msd.\n");
		    else
			printf("Don't compute msd.\n");
		    break;
		case(7):
		    strcpy(outputfile, value);
		    printf("outputfile -> %s\n", outputfile);
		    break;
		case(8):
		    strcpy(msdout, value);
		    printf("msdout -> %s\n", msdout);
		    break;
		case(9):
		    strcpy(adsorb, value);
		    printf("adsorb -> %s\n", adsorb);
		    break;

		default:
		    fprintf(stderr, "Internal error\n");
		    exit(EXIT_FAILURE);
	    }			
	    readstatus = fscanf(cfg, "%s = %s\n", keyword, value);
	    lineno++;
	}
	printf("\n");
	
	fclose(cfg);
	
	position_array = allocate_double_matrix(nochains*chainlength, 3);
	cm_array = allocate_3d_double_matrix(notimesteps,  nochains, 3);
	
        length = strlen(dumpname);
        if (length > 3)
        {
                if (strcmp(&dumpname[length-3], ".gz") == 0)
                        compressed=TRUE;
                else
                        compressed=FALSE;
        }
        else
        {
                compressed=FALSE;
        }
	if ( compressed == TRUE )
        {
                strcpy(command, "gzip -dc ");
                strcat(command, dumpname);
                fp = popen(command, "r");
        }
        else
        {
                fp  = fopen(dumpname, "r");
        }

	if (fp==NULL)
	{
		fprintf(stderr, "Cannot open dumpfile %s\n", dumpname);
		exit(EXIT_FAILURE);
	}

	
	out1=fopen(outputfile, "w");
	
	fprintf(out1, "# timestep    xRg^2        yRg^2           zRg^2           Rg^2          xRe^2         yRe^2         zRe^2          Re^2            xl^2           yl^2          zl^2           l^2\n");
	    
	out3=fopen(adsorb, "w");
	fprintf(out3, "# tstep admon(/%d) ad_percent adchains trainmon tailmon loopmon nonadmon ave_trainmon ave_tailmon ave_loopmon ave_nonadmon no_train  no_tail   no_loop\n", chainlength*nochains);

	aveRgsum=0;
	
	for(timestep=0; timestep<notimesteps; timestep++)
	{
	    
	    readconfig(position_array, nochains, chainlength);
	    
	    /*printconfig(position_array, nochains*chainlength);*/
	    
	    calccentermass(cm_array, position_array, chainlength, nochains, timestep);

	    aveavel=calcmonomerdistance(position_array, chainlength, nochains, aveavel, coordaveavel);

	    aveRe=calcendlength(position_array, chainlength, nochains, aveRe, coordaveRe);	
	    
	    aveRg=calcgyration(cm_array, position_array, chainlength, nochains, timestep, coordaveRg);
	    
	    looptraintail(position_array, chainlength, nochains, ad_condition, &no_train_mon, &no_tail_mon, &no_loop_mon, &no_nonadsorb_mon, &avetrain_mon, &avetail_mon, &aveloop_mon, &avenonadsorb_mon, &no_train, &no_tail, &no_loop, &no_nonadsorb);

	    /*printf("%d %d %d %d\n", nonadsorb, train, tail, loop);*/

	    total_nomon_ad=adsorptioncondition(position_array, chainlength, nochains, total_nomon_ad, ad_condition, &total_nochain_ad);
	 
	    fprintf(out3, "    %d     %d %14.4f      %d        %d       %d        %d       %d    %10.4f  %10.4f  %10.4f   %10.4f %10.4f %10.4f %10.4f\n", timestep, total_nomon_ad, (double)total_nomon_ad/(chainlength*nochains), total_nochain_ad, no_train_mon, no_tail_mon, no_loop_mon, no_nonadsorb_mon, avetrain_mon, avetail_mon, aveloop_mon, avenonadsorb_mon, (double)no_train/nochains, (double)no_tail/nochains, (double)no_loop/nochains);

	    fprintf(out1, "  %d  %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f %14.10f\n", timestep, coordaveRg[0], coordaveRg[1], coordaveRg[2], aveRg, coordaveRe[0], coordaveRe[1], coordaveRe[2], aveRe, coordaveavel[0], coordaveavel[1], coordaveavel[2], aveavel);
	    
	    aveRgsum+=aveRg;
	}
	aveaveRg=aveRgsum/notimesteps;
	/*printf("Average radius of gyration over all timesteps: %14.10f\n", aveaveRg);*/   
	
       	if(msdon==1)
	{
	    out2=fopen(msdout, "w");
	    fprintf(out2, "# delta t    msd(x)        error            msd(y)      error            msd(z)         error             mad        error          msd(x+y)      error         no.timesteps(%d)\n", notimesteps);
	    
	    calcmsd(cm_array, notimesteps, nochains, maxdelta);
	    fclose(out2);
	}
		
	if (compressed == TRUE)
	{
	    pclose(fp);
	}
	else
	{
	    fclose(fp);
	}
	
	fclose(out1);
	fclose(out3);	

	free_3d_double_matrix(cm_array);
	free_double_matrix(position_array);
	
	exit(0);
}

void calcmsd(double ***cm_array, int notimesteps, int nchains, int maxdelta)
{
	int i, j, delta, chain, deltat, nodeltas, max_j, max_nodeltas; 
	double disp[3], **sumd, **msd, **xymsd, ***coordmsd; 
	double **sqrsumd, ***sqrcoordmsd, **sqrmsd, **sqrxymsd;
	double summsd, sumxymsd, coordsummsd[3]; 
	double sqrsummsd, sqrsumxymsd, sqrcoordsummsd[3];
	double varmsd, varxymsd, varcoordmsd[3]; 
	double errmsd, errxymsd, errcoordmsd[3];
	max_nodeltas=notimesteps-1;
	double number_of_samples[max_nodeltas];	
	double tempmsd[max_nodeltas], tempxymsd[max_nodeltas], tempsqrmsd[max_nodeltas], tempsqrxymsd[max_nodeltas], nodelta[max_nodeltas];
       

	sumd = allocate_double_matrix(notimesteps, 3);
	sqrsumd = allocate_double_matrix(notimesteps, 3);
	msd = allocate_double_matrix(nchains, notimesteps);
	sqrmsd = allocate_double_matrix(nchains, notimesteps);
	xymsd = allocate_double_matrix(nchains, notimesteps);
	sqrxymsd = allocate_double_matrix(nchains, notimesteps);
	coordmsd = allocate_3d_double_matrix(nchains, notimesteps, 3);
	sqrcoordmsd = allocate_3d_double_matrix(nchains, notimesteps, 3);

	for(chain=0; chain<nchains; chain++)
	{
	    for(delta=0; delta<maxdelta; delta++)
	    {
		sumd[delta][0]=0;
		sumd[delta][1]=0;
		sumd[delta][2]=0;
		sqrsumd[delta][0]=0;
		sqrsumd[delta][1]=0;
		sqrsumd[delta][2]=0;
		tempmsd[delta]=0;
		tempxymsd[delta]=0;
		tempsqrmsd[delta]=0;
		tempsqrxymsd[delta]=0;
		nodelta[delta]=0;

		number_of_samples[delta]=0;
	    }
	    /*printf("%g\n", number_of_samples[0]);*/

	    for(i=0; i<(notimesteps-1); i++)
	    {
		max_j=i+maxdelta;
		
		if(max_j > (notimesteps-1))
		    max_j=notimesteps-1;
	
		for(j=i; j<max_j; j++)
		{
		    deltat=j-i;
		    disp[0]=cm_array[i][chain][0]-cm_array[j+1][chain][0];
		    disp[1]=cm_array[i][chain][1]-cm_array[j+1][chain][1];
		    disp[2]=cm_array[i][chain][2]-cm_array[j+1][chain][2];
		    
		    sumd[deltat][0]+=disp[0]*disp[0];
		    sumd[deltat][1]+=disp[1]*disp[1];
		    sumd[deltat][2]+=disp[2]*disp[2];

		    sqrsumd[deltat][0]+=disp[0]*disp[0]*disp[0]*disp[0];
		    sqrsumd[deltat][1]+=disp[1]*disp[1]*disp[1]*disp[1];
		    sqrsumd[deltat][2]+=disp[2]*disp[2]*disp[2]*disp[2];
		 
		    tempmsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
		    tempxymsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1]);

		    tempsqrmsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2])*(disp[0]*disp[0]+disp[1]*disp[1]+disp[2]*disp[2]);
		    tempsqrxymsd[deltat]+=(disp[0]*disp[0]+disp[1]*disp[1])*(disp[0]*disp[0]+disp[1]*disp[1]);

		    number_of_samples[deltat]++;
		}    
	    }

	    for(delta=0; delta<maxdelta; delta++)
	    { 
		nodeltas=notimesteps-(delta+1);
		
		if(number_of_samples[delta]!=nodeltas)
		 printf("Ooops, It's wrong number of deltat!\n");		
		coordmsd[chain][delta][0]=sumd[delta][0]/nodeltas;
		coordmsd[chain][delta][1]=sumd[delta][1]/nodeltas;
		coordmsd[chain][delta][2]=sumd[delta][2]/nodeltas;

		sqrcoordmsd[chain][delta][0]=sqrsumd[delta][0]/nodeltas;
		sqrcoordmsd[chain][delta][1]=sqrsumd[delta][1]/nodeltas;
		sqrcoordmsd[chain][delta][2]=sqrsumd[delta][2]/nodeltas;
		
		msd[chain][delta]=tempmsd[delta]/nodeltas;
		sqrmsd[chain][delta]=tempsqrmsd[delta]/nodeltas;

		xymsd[chain][delta]=tempxymsd[delta]/nodeltas;
		sqrxymsd[chain][delta]=tempsqrxymsd[delta]/nodeltas;

		nodelta[delta]=nodeltas;
	    } 
	}
	
	for(delta=0; delta<maxdelta; delta++)
	{
	    coordsummsd[0]=0;
	    coordsummsd[1]=0;
	    coordsummsd[2]=0;
	    sqrcoordsummsd[0]=0;
	    sqrcoordsummsd[1]=0;
	    sqrcoordsummsd[2]=0;

	    summsd=0;
	    sumxymsd=0;
	    sqrsummsd=0;
	    sqrsumxymsd=0;
	    
	    for(chain=0; chain<nchains; chain++)
	    {  
		coordsummsd[0]+=coordmsd[chain][delta][0];
		coordsummsd[1]+=coordmsd[chain][delta][1];
		coordsummsd[2]+=coordmsd[chain][delta][2];

		sqrcoordsummsd[0]+=sqrcoordmsd[chain][delta][0];
		sqrcoordsummsd[1]+=sqrcoordmsd[chain][delta][1];
		sqrcoordsummsd[2]+=sqrcoordmsd[chain][delta][2];

		summsd+=msd[chain][delta];
		sumxymsd+=xymsd[chain][delta];

		sqrsummsd+=sqrmsd[chain][delta];
		sqrsumxymsd+=sqrxymsd[chain][delta];

	    }
	    
	    coordsummsd[0]/=nchains;
	    coordsummsd[1]/=nchains;
	    coordsummsd[2]/=nchains;

	    sqrcoordsummsd[0]/=nchains;
	    sqrcoordsummsd[1]/=nchains;
	    sqrcoordsummsd[2]/=nchains;

	    summsd/=nchains;
	    sumxymsd/=nchains;

	    sqrsummsd/=nchains;
	    sqrsumxymsd/=nchains;

	    varcoordmsd[0]=sqrcoordsummsd[0]-coordsummsd[0]*coordsummsd[0];
	    varcoordmsd[1]=sqrcoordsummsd[1]-coordsummsd[1]*coordsummsd[1];
	    varcoordmsd[2]=sqrcoordsummsd[2]-coordsummsd[2]*coordsummsd[2];

	    varmsd=sqrsummsd-summsd*summsd;
	    varxymsd=sqrsumxymsd-sumxymsd*sumxymsd;

	    if (nchains==1)
	    {
	    errcoordmsd[0]=sqrt(varcoordmsd[0])/sqrt(nchains);
	    errcoordmsd[1]=sqrt(varcoordmsd[1])/sqrt(nchains);
	    errcoordmsd[2]=sqrt(varcoordmsd[2])/sqrt(nchains);
	    errmsd=sqrt(varmsd)/sqrt(nchains);
	    errxymsd=sqrt(varxymsd)/sqrt(nchains);
	    }
	    else
	    {
	    errcoordmsd[0]=sqrt(varcoordmsd[0])/sqrt(nchains-1);
	    errcoordmsd[1]=sqrt(varcoordmsd[1])/sqrt(nchains-1);
	    errcoordmsd[2]=sqrt(varcoordmsd[2])/sqrt(nchains-1);
	    errmsd=sqrt(varmsd)/sqrt(nchains-1);
	    errxymsd=sqrt(varxymsd)/sqrt(nchains-1);
	    }
	    
	    /* printf("%g %g %g %g %g\n", sqrcoordsummsd[0], sqrcoordsummsd[1], sqrcoordsummsd[2], sqrsummsd, sqrsumxymsd);
	    printf("%g %g %g %g %g\n", varcoordmsd[0], varcoordmsd[1], varcoordmsd[2], varmsd, varxymsd);
	    printf("  %d  %14.10f %14.10f %14.10f %14.10f %14.10f\n", delta+1, coordsummsd[0], coordsummsd[1], coordsummsd[2], summsd, sumxymsd);
	    printf("%g %g %g %g %g\n", errcoordmsd[0], errcoordmsd[1], errcoordmsd[2], errmsd, errxymsd);*/

	    fprintf(out2, "    %d %14.5f %14.10f %14.5f %14.10f %14.5f %14.10f %14.5f %14.10f %14.5f %14.10f\n", delta+1, coordsummsd[0], errcoordmsd[0], coordsummsd[1], errcoordmsd[1], coordsummsd[2], errcoordmsd[2], summsd, errmsd, sumxymsd, errxymsd);
	    
	}	
	
	free_double_matrix(sumd);
	free_double_matrix(msd);
	free_double_matrix(xymsd);
	free_3d_double_matrix(coordmsd);
	free_double_matrix(sqrsumd);
	free_double_matrix(sqrmsd);
	free_double_matrix(sqrxymsd);
	free_3d_double_matrix(sqrcoordmsd);

	return;
}

int readconfig(double **positions, int nchains, int chainlength)
{
	int chain, monomer, current_atom;
	int i, headerstatus, timestep, natoms, atomtype;
	double box[3][2], xbox, ybox, zbox;
	
	int readatom(double *atom, double xbox, double ybox, double zbox, int *atomtype);
	int readheader(int *timestep, int *natoms, double box[3][2], int nchains, int chainlength);
	
	headerstatus = readheader(&timestep, &natoms, box, nchains, chainlength);
	if ( headerstatus != 0)
	{
		printf("Could not read header! Status=%d\n", headerstatus);
		exit(EXIT_FAILURE);
	}
	
	/*printf("\n");    
	printf("Processing timestep: %d\n", timestep);
	printf("Total no of atoms: %d\n", natoms);*/
	xbox=box[0][1]-box[0][0];
	ybox=box[1][1]-box[1][0];
	zbox=box[2][1]-box[2][0];

	
	for(chain=0; chain<nchains; chain++)
	{
		for(monomer=0; monomer<chainlength; monomer++)
		{
			current_atom=chain*chainlength+monomer;
			i = readatom(positions[current_atom], xbox, ybox, zbox, &atomtype);
			if (current_atom+1 != i)
			{
			    printf("The configuration of dump is disordered!!!\n Please sort it!!!\n");
			    exit(EXIT_FAILURE);
			}
			/*printf("I want to read atom %d. I actually read atom no. %d\n", current_atom, i);*/
		}
	}
/*remove the function to read surface from dump file*/

	/*for(surface=0; surface<(natoms-nchains*chainlength); surface++)
	{
	    fscanf(fp, "%d %d %lg %lg %lg %d %d %d\n", 
		   &surfaceatomno, &surftype, &surfx, &surfy, &surfz, &surfdummy[0], &surfdummy[1], &surfdummy[2]);
	    
		printf("%d\n", atomtype);
	    if (surftype==atomtype)
	      {
		  printf("Wrong input no. of chains or chainlength!\n");
		  exit(EXIT_FAILURE);
	      }

	}*/
	return(0);
}

int readheader(int *timestep, int *natoms, double box[3][2], int nchains, int chainlength)
{
	int i;
	char line[1024];
	
	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: TIMESTEP\n") != 0)
	{   
                /*printf("%s\n", line);*/
		return(-1);
	}
	fscanf(fp, "%d\n", timestep);
	
	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: NUMBER OF ATOMS\n") != 0)
	{
		return(-2);
	}
	fscanf(fp, "%d\n", natoms);

	if (*natoms != (nchains*chainlength))
	{
	    printf("Wrong chainlength or no_chains!\n");
	    exit(0);
	}

	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: BOX BOUNDS\n") != 0)
	{
		return(-3);
	}
	for(i=0; i<3; i++)
	{
		fscanf(fp, "%lg %lg\n", &box[i][0], &box[i][1]);
	}
	
	fgets(line, 1024, fp);
	if ( strcmp(line, "ITEM: ATOMS\n") != 0)
	{
		return(-4);
	}
	
	return(0);
}

void printconfig(double **positions, int atoms)
{
	int atomno;
	
	for(atomno=0; atomno<atoms; atomno++)
	{
		printf("Atom %d: x=%g y=%g z=%g\n",
		       atomno+1, positions[atomno][0], positions[atomno][1], positions[atomno][2]);
	}
	
	return;
}

int readatom(double *atom, double xbox, double ybox, double zbox, int *atomtype)
{
	int atomno, period[3];
	double x, y, z;
	
	fscanf(fp, "%d %d %lg %lg %lg %d %d %d\n",
	       &atomno, atomtype, &x, &y, &z, &period[0], &period[1], &period[2]);
	atom[0]=x+xbox*period[0];
	atom[1]=y+ybox*period[1];
	atom[2]=z+zbox*period[2];
	
	return(atomno);
}

double calcgyration(double ***cm_array, double **positions, int chainlength, int nchains, int timestep, double coordaveRg[3])
/*Referece: polymer physics; Rubinstein p.60*/
{
	int chain, i, atomno; 
	double aveRg, x, y, z, xsum, ysum, zsum, Rg, Rgsum, coordRgsum[3];

	Rgsum=0;
	coordRgsum[0]=0;
	coordRgsum[1]=0;
	coordRgsum[2]=0;

	for(chain=0; chain<nchains; chain++)
	{
		xsum=0;
		ysum=0;
		zsum=0;
		
		for(i=0; i<chainlength; i++)
		{
			atomno=i+chainlength*chain;
			
			x=positions[atomno][0]-cm_array[timestep][chain][0];
			y=positions[atomno][1]-cm_array[timestep][chain][1];
			z=positions[atomno][2]-cm_array[timestep][chain][2];
			
			xsum+=x*x;
			ysum+=y*y;
			zsum+=z*z;
		}

		Rg=(xsum+ysum+zsum)/chainlength;
		/*printf("Chain %d's Radius of gyration: %14.10f\n", chain+1, Rg);*/

		coordRgsum[0]+=xsum/chainlength;
		coordRgsum[1]+=ysum/chainlength;
		coordRgsum[2]+=zsum/chainlength;

		Rgsum+=Rg;
		/*printf("%g %g\n",Rgsum, Rg);*/
		
	}

	coordaveRg[0]=coordRgsum[0]/nchains;
	coordaveRg[1]=coordRgsum[1]/nchains;
	coordaveRg[2]=coordRgsum[2]/nchains;

	aveRg=Rgsum/nchains;
	/*printf("Average Radius of gyration: %14.10f\n", aveRg);*/
	
	return(aveRg);
}

void calccentermass(double ***cm_array, double **positions, int chainlength, int nchains, int timestep)
{
	int chain;
	void calccm(double *cm, double **positions, int chainlength);
	
	for(chain=0; chain<nchains; chain++)
	{
		calccm(cm_array[timestep][chain], positions+chain*chainlength, chainlength);

		/*printf("chain %d's center of mass: (%14.10f, %14.10f, %14.10f)\n", 
			chain+1, cm_array[timestep][chain][0], cm_array[timestep][chain][1], cm_array[timestep][chain][2]);*/
	}

	return;
}

void calccm(double cm[3], double **positions, int chainlength)
{
	int atomno;
	double xsum, ysum, zsum;
	xsum=0;
	ysum=0;
	zsum=0;
	
	for(atomno=0; atomno<chainlength; atomno++)
	{	
		xsum+=positions[atomno][0];
		ysum+=positions[atomno][1];
		zsum+=positions[atomno][2];
	}
	cm[0]=xsum/chainlength;
	cm[1]=ysum/chainlength;
	cm[2]=zsum/chainlength;

	return;
}

double calcendlength(double **positions, int chainlength, int nchains, double aveRe, double coordaveRe[3])
{
	double xend, yend, zend, Re, Resum, coordResum[3];
	int chain, n1, n2;
    

	coordResum[0]=0;
	coordResum[1]=0;
	coordResum[2]=0;
	
	Resum=0;
	
	for(chain=0; chain<nchains; chain++)
	{
		n1=chainlength*chain+chainlength-1;
		n2=chainlength*chain;	
		
		xend=positions[n1][0]-positions[n2][0];
		yend=positions[n1][1]-positions[n2][1];
		zend=positions[n1][2]-positions[n2][2];

		Re=xend*xend+yend*yend+zend*zend;
		
		/*printf("Chain %d's end-to-end length: %14.10f\n", chain+1, Re);*/

		coordResum[0]+=xend*xend;
		coordResum[1]+=yend*yend;
		coordResum[2]+=zend*zend;

		Resum+=Re;
	}

	aveRe=Resum/nchains;

	coordaveRe[0]=coordResum[0]/nchains;
	coordaveRe[1]=coordResum[1]/nchains;
	coordaveRe[2]=coordResum[2]/nchains;

	return(aveRe);
}	
double calcmonomerdistance(double **positions, int chainlength, int nchains, double aveavel, double coordaveavel[3])
{
    double xdist, ydist, zdist;
    double coordl[3], coordavel[3], distancel, suml, avel, sumcoordavel[3], sumavel;
    int chainno, atomno;

    sumcoordavel[0]=0;
    sumcoordavel[1]=0;
    sumcoordavel[2]=0;
    sumavel=0;
    
    for(chainno=0; chainno<nchains; chainno++)
    {
	coordl[0]=0;
	coordl[1]=0;
	coordl[2]=0;
	suml=0;

	for(atomno=0; atomno<(chainlength-1); atomno++)
	{
	    xdist=positions[chainno*chainlength+atomno+1][0]-positions[chainno*chainlength+atomno][0];
	    ydist=positions[chainno*chainlength+atomno+1][1]-positions[chainno*chainlength+atomno][1];
	    zdist=positions[chainno*chainlength+atomno+1][2]-positions[chainno*chainlength+atomno][2];

	    xdist=xdist*xdist;
	    ydist=ydist*ydist;
	    zdist=zdist*zdist;

	    distancel=xdist+ydist+zdist;

	    coordl[0]+=xdist;
	    coordl[1]+=ydist;
	    coordl[2]+=zdist;
	    suml+=distancel;

	}

	coordavel[0]=coordl[0]/(chainlength-1);
	coordavel[1]=coordl[1]/(chainlength-1);
	coordavel[2]=coordl[2]/(chainlength-1);
	avel=suml/(chainlength-1);

	sumcoordavel[0]+=coordavel[0];
	sumcoordavel[1]+=coordavel[1];
	sumcoordavel[2]+=coordavel[2];
	sumavel+=avel;

    }

    coordaveavel[0]=sumcoordavel[0]/nchains;
    coordaveavel[1]=sumcoordavel[1]/nchains;
    coordaveavel[2]=sumcoordavel[2]/nchains;
    aveavel=sumavel/nchains;

    return(aveavel);
}

int adsorptioncondition(double **positions, int chainlength, int nchains, int total_nomon_ad, double ad_condition, int *total_nochain_ad)
{
    int chainno, atomno, nomon_ad, nochain_ad;

    total_nomon_ad=0;
    *total_nochain_ad=0;

    for(chainno=0; chainno<nchains; chainno++)
    {
	nomon_ad=0;
	nochain_ad=0;

	for(atomno=0; atomno<chainlength; atomno++)
	{
	    if (positions[chainno*chainlength+atomno][2] <= ad_condition)
	    {
		nomon_ad++;
	    }
	}

	if (nomon_ad >= 1)
	{
	    nochain_ad++;
	}
	total_nomon_ad+=nomon_ad;
	*total_nochain_ad+=nochain_ad;
    }

    /*printf("%d\n", *total_nochain_ad);*/
    return(total_nomon_ad);
}

void looptraintail(double **positions, int chainlength, int nchains, double ad_condition, int *no_train_mon, int *no_tail_mon, int *no_loop_mon, int *no_nonadsorb_mon, double *avetrain_mon, double *avetail_mon, double *aveloop_mon, double *avenonadsorb_mon, int *no_train, int *no_tail, int *no_loop, int *no_nonadsorb)
{

    int chainno, atomno, index;
    int monomer[chainlength*nchains];
    int train[chainlength*nchains], tail[chainlength*nchains], loop[chainlength*nchains], nonadsorb[chainlength*nchains];
    int a, b, c, d;
        
    for(chainno=0; chainno<nchains; chainno++)
    {
	for(atomno=0; atomno<chainlength; atomno++)
	{
	    index=chainno*chainlength+atomno;

	    if (positions[index][2] > ad_condition)
	    {
		monomer[index]=0;

		if (monomer[index-1]==1)
		{
		    if (index == chainno*chainlength)
			goto end;
	
		    monomer[index]=2;
		}
		else if (monomer[index-1]==2)
		{
		    if (index == chainno*chainlength)
			goto end;	

		    monomer[index]=2;
		}
	    }

	    else if (positions[index][2] <= ad_condition)
	    {
		monomer[index]=1;
		
		if (monomer[index-1]==0)
		{
		    while (monomer[index-1]==0)
		    {
			if(index == chainno*chainlength)
			    goto end;
				
			monomer[index-1]=2;
			index--;
		    }
		}
		else if (monomer[index-1]==2)
		{
		    while (monomer[index-1]==2)
		    {
			if(index == chainno*chainlength)
			    goto end;
					
			monomer[index-1]=3;
			index--;
		    }
		}
	    }
	end: ;
	}
    }  

    *no_train=*no_tail=*no_loop=*no_nonadsorb=0;
    
    for(chainno=0; chainno<nchains; chainno++)
    {
    	for(atomno=0; atomno<chainlength; atomno++)
	{ 
	    index=chainno*chainlength+atomno;
	    
	    switch(monomer[index])
	    {
		case(0):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit0;
			
			nonadsorb[*no_nonadsorb]+=1;
		    }
		    else
		    {
		    exit0:
		
			*no_nonadsorb+=1;
			nonadsorb[*no_nonadsorb]=1;
		    }
		    break;

		case(1):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit1;
			
			train[*no_train]+=1;
		    }
	 
		    else 
		    {
		    exit1:
			
			*no_train+=1;
			train[*no_train]=1;
		    }
		    break;

		case(2):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit2;
			
			tail[*no_tail]+=1;
		    }
		    else
		    {
		    exit2:
			
			*no_tail+=1;
			tail[*no_tail]=1;
		    }
		    break;

		case(3):

		    if(monomer[index]== monomer[index-1])
		    {
			if (index == chainno*chainlength)
			    goto exit3;
			
			loop[*no_loop]+=1;
		    }
		    else
		    {
		    exit3:
			
			*no_loop+=1;
			loop[*no_loop]=1;
		    }
		    break;
	    }
	}
    }

    /*printf("%d %d %d %d\n", no_nonadsorb, no_train, no_tail, no_loop);*/
   
    *no_nonadsorb_mon=*no_train_mon=*no_tail_mon=*no_loop_mon=0;
    *avenonadsorb_mon=*avetrain_mon=*avetail_mon=*aveloop_mon=0;
   
    for(a=0; a<*no_nonadsorb; a++)
    {
	*no_nonadsorb_mon+=nonadsorb[a+1];
    }
    if(*no_nonadsorb != 0)
	*avenonadsorb_mon=(float) *no_nonadsorb_mon/ *no_nonadsorb;
 
    for(b=0; b<*no_train; b++)
    {
	*no_train_mon+=train[b+1];
    }
    if(*no_train != 0)
    *avetrain_mon=(float) *no_train_mon/ *no_train;
    
    for(c=0; c<*no_tail; c++)
    {
	*no_tail_mon+=tail[c+1];
    }
    if(*no_tail != 0)
    *avetail_mon=(float) *no_tail_mon/ *no_tail;

    for(d=0; d<*no_loop; d++)
    {
	*no_loop_mon+=loop[d+1];
    }
    if(*no_loop != 0)
    *aveloop_mon=(float) *no_loop_mon/ *no_loop;

    /*printf("%d %d %d %d\n", no_nonadsorb_mon, no_train_mon, no_tail_mon, no_loop_mon);*/

     /* printf("%10.4f %10.4f %10.4f %10.4f\n", avenonadsorb_mon, avetrain_mon, avetail_mon, aveloop_mon);*/

    return;
}
