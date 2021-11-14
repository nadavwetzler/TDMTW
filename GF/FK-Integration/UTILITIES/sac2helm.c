#include<sys/file.h>
#include<stdlib.h>
#include<fcntl.h>
#include<stdio.h>
#include <math.h>

FILE *fopen(), *fd4;

/*convert SAC binary to headerless binary file (fromHelm.c generated) */
/*or Helmberger ascii */
main(ac,av)
 int ac;
 char **av;
 {
 int i, npts, ihd[40],fd1, fd2, fd3;
 float dt, fhd[70], *tr;
 char chd[8][24], out[100], line[100], line1[100], line2[100];

setpar(ac,av);
mstpar("out","s",out);
endpar();

fd1=open("tmp1",O_RDONLY,0644);
fd2=open("tmp2",O_RDONLY,0644);
fd3=open("tmp3",O_RDONLY,0644);
fd4=fopen(out,"w");


read(fd1,fhd,70*4);  /*Read Sac Float Field*/
read(fd1,ihd,40*4);  /*Read Sac Int   Field*/
read(fd1,chd,24*8);  /*Read Sac Char. Field*/
npts=ihd[9];
dt=fhd[0];
read(fd2,fhd,70*4);  /*Read Sac Float Field*/
read(fd2,ihd,40*4);  /*Read Sac Int   Field*/
read(fd2,chd,24*8);  /*Read Sac Char. Field*/
if(ihd[9] != npts)
    {
    fprintf(stderr,"SAC2HELM ERROR npts not equal\n");
    exit(-1);
    }
if(fhd[0] != dt)
    {
    fprintf(stderr,"SAC2HELM ERROR dt not equal\n");
    exit(-1);
    }
read(fd3,fhd,70*4);  /*Read Sac Float Field*/
read(fd3,ihd,40*4);  /*Read Sac Int   Field*/
read(fd3,chd,24*8);  /*Read Sac Char. Field*/
if(ihd[9] != npts)
    {
    fprintf(stderr,"SAC2HELM ERROR npts not equal\n");
    exit(-1);
    }
if(fhd[0] != dt)
    {
    fprintf(stderr,"SAC2HELM ERROR dt not equal\n");
    exit(-1);
    }
fprintf(stderr,"npts=%d   dt=%f\n",npts,dt);

sprintf(line1,"     0.0000e+00     0.0000e+00      0  0  0.00");
sprintf(line2,"%8d  %8.5f %11.4e",npts,dt,0.0);

tr=(float *)malloc(sizeof(float)*npts);
read(fd1,tr,npts*sizeof(float));
fprintf(fd4,"       3\n");
fprintf(fd4,"(7e14.5)\n");
fprintf(fd4,"%s\n",line1);
fprintf(fd4,"%s\n",line2);
for(i=1; i <= npts; i++)
   if(i % 7 == 0 && i != 1 && i != npts)
     fprintf(fd4,"  %12.5e\n",tr[i-1]);
   else
     fprintf(fd4,"  %12.5e",tr[i-1]);

read(fd2,tr,npts*sizeof(float));
fprintf(fd4,"\n%s\n",line1);
fprintf(fd4,"%s\n",line2);
for(i=1; i <= npts; i++)
   if(i % 7 == 0 && i != 1 && i != npts)
     fprintf(fd4,"  %12.5e\n",tr[i-1]);
   else
     fprintf(fd4,"  %12.5e",tr[i-1]);

read(fd3,tr,npts*sizeof(float));
fprintf(fd4,"\n%s\n",line1);
fprintf(fd4,"%s\n",line2);
for(i=1; i <= npts; i++)
   if(i % 7 == 0 && i != 1 && i != npts)
     fprintf(fd4,"  %12.5e\n",tr[i-1]);
   else
     fprintf(fd4,"  %12.5e",tr[i-1]);


}
