#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

double bathtubC(double age) {
	double p = 0.6*exp(-age/4.)+(-1.+exp(age*log(2.)/20.));
	if(p>1.)
		{
			p=1.;
		}
	return p;
}
