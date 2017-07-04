#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <omp.h>

double euc(double *x, double *y, int len) {
	double sum = 0.0;
	int count = 0;
	for (int i=0; i<len; ++i) {
	  if ((x[i] != -1) && (y[i] != -1)) {
	    sum += (x[i]-y[i])*(x[i]-y[i]);
	    count += 1;
	  }
	}	

	if (count == 0) {
		return -1;
	}
	double res = sqrt(sum/count);
	return res;
	
}

double na(double *x, double *y, int len) {
  int count = 0;
  for (int i=0; i<len; ++i) {
    if ((x[i] != -1) && (y[i] != -1)) {
      count += 1;
    }
  }

  return count;

}


double man(double *x, double *y, int len) {
  double sum = 0.0;
  int count = 0;
  for (int i=0; i<len; ++i) {
    if ((x[i] != -1) && (y[i] != -1)) {
      sum += abs(x[i]-y[i]);
      count += 1;
    }
  }

  if (count == 0) {
    return -1;
  }
  double res = sum/count;
  return res;

}

double allel(double *x, double *y, int len, int diff) {
  double sum = 0.0;
  int count = 0;
  for (int i=0; i<len; ++i) {
    if ((x[i] != -1) && (y[i] != -1)) {
      sum += (abs(x[i]-y[i])>diff)?1:0;
      count += 1;
    }
  }

  if (count == 0) {
    return -1;
  }
  double res = sum/count;
  return res;

}

void euc_dist(double *a, int *dim, double* res) {
	int nrow,ncol;
	nrow = dim[0];
	ncol = dim[1];

	#pragma omp parallel for schedule(dynamic, 5)
	for (int i=0; i<ncol;++i) {
		for( int j=0; j<i; ++j) {
			res[i*ncol+j] = euc(a+i*nrow,a+j*nrow,nrow);
		}
	}
}

void manhattan_dist(double *a, int *dim, double* res) {
  int nrow,ncol;
  nrow = dim[0];
  ncol = dim[1];

#pragma omp parallel for schedule(dynamic, 5)
  for (int i=0; i<ncol;++i) {
    for( int j=0; j<i; ++j) {
      res[i*ncol+j] = man(a+i*nrow,a+j*nrow,nrow);
    }
  }
}


void NA_dist(double *a, int *dim, double* res) {
  int nrow,ncol;
  nrow = dim[0];
  ncol = dim[1];

#pragma omp parallel for schedule(dynamic, 5)
  for (int i=0; i<ncol;++i) {
    for( int j=0; j<i; ++j) {
      res[i*ncol+j] = na(a+i*nrow,a+j*nrow,nrow);
    }
  }
}


void allel_dist(double *a, int *dim, double* res, int *diff) {
  int nrow,ncol;
  nrow = dim[0];
  ncol = dim[1];

#pragma omp parallel for schedule(dynamic, 5)
  for (int i=0; i<ncol;++i) {
    for( int j=0; j<i; ++j) {
      res[i*ncol+j] = allel(a+i*nrow,a+j*nrow,nrow,*diff);
    }
  }
}
