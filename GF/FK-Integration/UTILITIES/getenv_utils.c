#ifndef lint
static char sccsid[] = "%W% %G% %U%";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/************************************************************************/
/* getenv_lf:								*/
/*	Get the value of an environment variable as a double.		*/
/*	Return 1 on success, 0 on failure.				*/
/************************************************************************/
int getenv_lf (var, pval)
    char *var;
    double *pval;
{
    char *str;
    if (str = getenv(var)) {
	char *p;
	double val;
	val = strtod (str, &p);
	if (*p == '\0') {
	    *pval = val;
	    return(1);
	}
    }
    return (0);
}

/************************************************************************/
/* getenv_f:								*/
/*	Get the value of an environment variable as a float.		*/
/*	Return 1 on success, 0 on failure.				*/
/************************************************************************/
int getenv_f (var, pval)
    char *var;
    float *pval;
{
    char *str;
    if (str = getenv(var)) {
	char *p;
	double val;
	val = (float)strtod (str, &p);
	if (*p == '\0') {
	    *pval = val;
	    return(1);
	}
    }
    return (0);
}

/************************************************************************/
/* getenv_d:								*/
/*	Get the value of an environment variable as a int.		*/
/*	Return 1 on success, 0 on failure.				*/
/************************************************************************/
int getenv_d (var, pval)
    char *var;
    int *pval;
{
    char *str;
    if (str = getenv(var)) {
	char *p;
	int val;
	val = strtol (str, &p, 10);
	if (*p == '\0') {
	    *pval = val;
	    return(1);
	}
    }
    return (0);
}

/************************************************************************/
/* getenv_s:								*/
/*	Get the value of an environment variable as a string.		*/
/*	Return 1 on success, 0 on failure.				*/
/************************************************************************/
int getenv_s (var, pval)
    char *var;
    char **pval;
{
    char *str;
    if (str = getenv(var)) {
	*pval = str;
	return(1);
    }
    return (0);
}

