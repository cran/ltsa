/**************************************************************************************
DLSim
input parameters:
	error - vector of length n, innovation sequence
	nR - length of time series
	c - vector of length n containing the autocovariances
	EPSL - machine epsilon in double precision
	
	output parameters:
	z - vector of length n, simulated time series
	ifault -  0, OK ;  1, c is not p.d.    
***************************************************************************************/
  
#include "nrutil.h"
#include <math.h>

void DLSim(double *z, double *error, int *nR, double *c, double *EPSL, int *fault)  
{
	double sum;
	int i,j,k,n;
	double EPS;

	VECTOR v,phi,phiki,phikj;
    
	*fault = 0;
	n = *nR;
	EPS = *EPSL;
	if (n < 1) *fault = 1;
	v = Vector(n);
	phi = Vector(n);
	phiki = Vector(n);
	phikj = Vector(n);
	
	v[0] = c[0];
	z[0] = error[0]*sqrt(v[0]);
	
	if (c[0] <= EPS) *fault = 1;
	phi[0] = c[1] / c[0];
	phiki[0] = phi[0];
	v[1] = v[0] * (1.0 - phi[0]*phi[0]);
	z[1] = sqrt(v[1])*error[1] + phi[0] * z[0];
	
	if (v[1] <= EPS) *fault = 1;
	
	for (k = 2; k < n; k++)
	{
		sum = 0.0;
		for (i = 1; i < k; i++)
			sum += phiki[i - 1] * c[k - i];	         
		phi[k - 1] = (c[k] - sum) / v[k - 1]; 
		for (j = 1; j< k; j++)
		{
			phikj[j - 1] = phiki[j - 1] - phi[k - 1] * phiki[k - j - 1];
		}
		phikj[k - 1] = phi[k - 1];
		sum = 0.0;
		for (j = 1; j<=k; j++)
		{
			sum += phikj[j - 1] * z[k - j];
			phiki[j - 1] = phikj[j - 1];
		}
		v[k] = v[k - 1] * (1.0 - phi[k - 1]*phi[k - 1]);
		if (v[k] <= EPS) *fault = 1;
		error[k] *= sqrt(v[k]); 
		z[k] = error[k] + sum;
		}
	
	free_vector(v);
	free_vector(phi);     
	free_vector(phiki);   
	free_vector(phikj);	  
  
	return;    
}
