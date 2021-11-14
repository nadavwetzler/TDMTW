#include<sys/file.h>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <fcntl.h>

void convlv();


main(ac,av)
 int ac;
 char **av;
 {
 int i, npts1, npts2, N, test, ihd[40], ihdsave[40],fd1, fd2, fd3, conv=-1;
 float dt,dt1, dt2, fhd[70], fhdsave[70], h2olev=0.01;
 float tau=0.0, *tmp,*outtrace;
 double *tr1, *tr2, *ans;
 char chd[8][24], chdsave[8][24], in1[100], in2[100], out[100];

setpar(ac,av);
mstpar("in1","s",in1);
mstpar("in2","s",in2);
mstpar("out","s",out);
getpar("h2olev","f",&h2olev);
getpar("conv","d",&conv);
getpar("tau","f",&tau);
endpar();

fd1=open(in1,O_RDONLY,0644);
fd2=open(in2,O_RDONLY,0644);
fd3=open(out,O_CREAT | O_TRUNC | O_WRONLY,0644);

read(fd1,fhdsave,70*4);  /*Read Sac Float Field*/
read(fd1,ihdsave,40*4);  /*Read Sac Int   Field*/
read(fd1,chdsave,24*8);  /*Read Sac Char. Field*/
npts1=ihdsave[9];
dt1=fhdsave[0];
fprintf(stderr,"npts=%d dt=%f\n",npts1,dt1);

read(fd2,fhd,70*4);  /*Read Sac Float Field*/
read(fd2,ihd,40*4);  /*Read Sac Int   Field*/
read(fd2,chd,24*8);  /*Read Sac Char. Field*/
npts2=ihd[9];
if((npts2 % 2) == 0)
   {
   npts2 -= 1;
   fprintf(stderr,"Resp points changed to %d\n",npts2);
   }
dt2=fhd[0];
fprintf(stderr,"npts=%d dt=%f\n",npts2,dt2);

if(dt1 != dt2)
    {
    fprintf(stderr,"DECON ERROR dt not equal\n");
    exit(-1);
    }
    dt=dt1;

if (2*npts2 > npts1)
   test=2*npts2; 
else
   test=npts1;

N=(int)(log10((float)test)/log10(2.0) + 2.0);
N=pow(2.0,(float)N);
/*
tau=(float)npts1/2.0*dt1;
*/
fprintf(stderr,"npts1=%d 2*npts2=%d N=%d\n",npts1,2*npts2, N);


tr1=(double *)malloc(sizeof(double)*(N+1));
tr2=(double *)malloc(sizeof(double)*(N+1));
tmp=(float *)malloc(sizeof(float)*(N+1));
for(i=0; i < N; i++)
  tr1[i]=tr2[i]=tmp[i]=0.0;
ans=(double *)malloc(sizeof(double)*2*(N+1));
outtrace=(float *)malloc(sizeof(float)*2*(N+1));
read(fd1,tmp,npts1*sizeof(float));
for(i=0; i<N;i++)
   {
   tr1[i]=(double)tmp[i];
   tmp[i]=0.0;
   }
read(fd2,tmp,npts2*sizeof(float));
for(i=0; i<N;i++)
   tr2[i]=(double)tmp[i];

free(tmp);

convlv(tr1-1,N,tr2-1,2*npts2 - 1,conv,ans-1,h2olev,tau,dt);

/*ihd[9]=npts2;*/
ihdsave[9]=npts1;
if(conv == 1) ihdsave[9]=npts1 + npts2 -1;
write(fd3,fhdsave,70*4);  /*Read Sac Float Field*/
write(fd3,ihdsave,40*4);  /*Read Sac Int   Field*/
write(fd3,chdsave,24*8);  /*Read Sac Char. Field*/
for(i=0; i<2*N+1; i++)
  outtrace[i]=(float)ans[i];
write(fd3,outtrace,ihdsave[9]*sizeof(float));


}
