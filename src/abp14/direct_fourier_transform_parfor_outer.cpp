#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>

#include "tbb/parallel_for.h"

namespace hpce
{
namespace abp14
{
class direct_fourier_transform_parfor_outer
	: public fourier_transform
{

typedef tbb::blocked_range<unsigned> my_range_t;
char *v=getenv("HPCE_DIRECT_INNER_K");
unsigned K = v==NULL ? 16 : atoi(v);

protected:
	/*! We can do any size transform */
	virtual size_t calc_padded_size(size_t n) const
	{
		return n;
	}

	virtual void forwards_impl(
		size_t n,	const complex_t &wn,
		const complex_t *pIn,
		complex_t *pOut
	) const
	{
		assert(n>0);
				
		my_range_t range(0,n,K);
		const real_t PI2=2*3.1415926535897932384626433832795;

		// = -i*2*pi / n
		complex_t neg_im_2pi_n=-complex_t(0, 1)*PI2 / (real_t)n;
		
		// Fill the output array with zero
		std::fill(pOut, pOut+n, 0);

		auto f=[&](const my_range_t &chunk){
	        	for(unsigned kk=chunk.begin(); kk!=chunk.end(); kk++){ 
				for(size_t ii=0; ii<n; ii++){
					pOut[kk] += pIn[ii] * exp( neg_im_2pi_n * (double)kk * (double)ii );
				}
			}
		};
		tbb::parallel_for(range, f, tbb::simple_partitioner());
	}	
	virtual void backwards_impl(
		size_t n,	const complex_t &/*wn*/,
		const complex_t *pIn,
		complex_t *pOut
	) const
	{
		assert(n>0);
		my_range_t range(0,n,K);
	
		const real_t PI2=2*3.1415926535897932384626433832795;

		// = i*2*pi / n
		complex_t im_2pi_n=complex_t(0, 1)*PI2 / (real_t)n;
		const real_t scale=real_t(1)/n;

		// Fill the output array with zero
		std::fill(pOut, pOut+n, 0);


		
		auto f=[&](const my_range_t &chunk){
	        	for(unsigned kk=chunk.begin(); kk!=chunk.end(); kk++){ 
				for(size_t ii=0; ii<n; ii++){
					pOut[kk] += pIn[ii] * exp( im_2pi_n * (double)kk * (double)ii ) * scale;
				}
			}
		};
		tbb::parallel_for(range, f, tbb::simple_partitioner());
			
	}

public:
	virtual std::string name() const
	{ return "hpce.abp14.direct_fourier_transform"; }

	virtual bool is_quadratic() const
	{ return true; }
};

std::shared_ptr<fourier_transform> Create_direct_fourier_transform_parfor_outer()
{
	return std::make_shared<direct_fourier_transform_parfor_outer>();
}

}; // namespace abp14

}; // namespace hpce
