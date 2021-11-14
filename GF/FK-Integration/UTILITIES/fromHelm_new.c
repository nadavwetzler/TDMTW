#include<stdio.h>
#include<stdlib.h>
#include<string.h>
/*This program converts Helmberger Ascii to headerless binary. It assumes all
vectors have the same length*/

main()
 {
 int i, j, ntr, n;
 float *tr, *vec, dt;
 char line[120];

 fgets(line,120,stdin); /*number of traces line*/
 sscanf(line,"%d",&ntr);
 fprintf(stderr,"ntr=%d\n",ntr);
 fgets(line,120,stdin); /*format line*/
 fgets(line,120,stdin); /*variables line*/
 fgets(line,120,stdin); /*npts and dt line*/
 sscanf(line,"%d %f",&n,&dt);
 fprintf(stderr,"n=%d  dt=%g\n",n,dt);

 tr=(float *)malloc(sizeof(float)*n);
 vec=(float *)malloc(sizeof(float)*n*ntr);
 for(i=0; i<ntr; i++)
   {
   fprintf(stderr,"trace=%d\n",i);
   for(j=0; j<n; j++)
     {
     if(j == (n-1))
       fscanf(stdin,"%f\n",&tr[j]); /*add return for incomplete line*/
     else
       fscanf(stdin,"%f",&tr[j]);

     vec[i*n+j]=tr[j];
     }
   if(i < (ntr-1))
     {
     fgets(line,120,stdin); /*variables line*/
     fprintf(stderr,"%s",line);
     fgets(line,120,stdin); /*npts and dt line*/
     fprintf(stderr,"%s",line);
     }
   }
 write(1,vec,sizeof(float)*n*ntr);
 }
