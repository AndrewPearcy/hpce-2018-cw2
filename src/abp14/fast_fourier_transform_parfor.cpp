#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>
#include "tbb/parallel_for.h"


namespace hpce
{
namespace abp14
{
class fast_fourier_transform_parfor
	: public fourier_transform
{

typedef tbb::blocked_range<unsigned> my_range_t;

protected:
	char *v = getenv("HPCE_FTT_LOOP_K");
	size_t K = v==NULL ? 2 : atoi(v);
	/* Standard radix-2 FFT only supports binary power lengths */
	virtual size_t calc_padded_size(size_t n) const
	{
		assert(n!=0);

		size_t ret=1;
		while(ret<n){
			ret<<=1;
		}

		return ret;
	}

	virtual void recurse(
		size_t n,	const complex_t &wn,
		const complex_t *pIn, size_t sIn,
		complex_t *pOut, size_t sOut
	) const
	{
		assert(n>0);

		if (n == 1){
			pOut[0] = pIn[0];
		}else if (n == 2){
			pOut[0] = pIn[0]+pIn[sIn];
			pOut[sOut] = pIn[0]-pIn[sIn];
		}else{
			size_t m = n/2;

			recurse(m,wn*wn,pIn,2*sIn,pOut,sOut);
			recurse(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut);

			complex_t w=complex_t(1, 0);
			
			if(m <= K) {
			  for (size_t j=0;j<m;j++){
				  complex_t t1 = w*pOut[m+j];
				  complex_t t2 = pOut[j]-t1;
				  pOut[j] = pOut[j]+t1;                
				  pOut[j+m] = t2;                          
				  w = w*wn;
			  }
			} else {	
				tbb::parallel_for(tbb::blocked_range<unsigned>(0,m,K), [&](const tbb::blocked_range<unsigned> &chunk){
    	w=pow(wn,chunk.begin());
	for(unsigned i=chunk.begin(); i!=chunk.end(); i++){
					  complex_t t1 = w*pOut[m+i];
					  complex_t t2 = pOut[i]-t1;
					  pOut[i] = pOut[i]+t1;                
					  pOut[m+i] = t2;                          
					  w = w*wn;
				  }	
				}, tbb::simple_partitioner());
			}
/*
			if(m <= K) {
			  for (size_t j=0;j<m;j++){
				  complex_t t1 = w*pOut[m+j];
				  complex_t t2 = pOut[j]-t1;
				  pOut[j] = pOut[j]+t1;                
				  pOut[j+m] = t2;                          
				  w = w*wn;
			  }
			} else {	
				my_range_t range(0,m,K);	
				auto f=[&](const my_range_t &chunk){
						
				  for (unsigned j=chunk.begin();j!=chunk.end();j++){
				    w = pow(wn,j);
				    complex_t t1 = w*pOut[m+j];
				    complex_t t2 = pOut[j]-t1;
				    pOut[j] = pOut[j]+t1;                 
				    pOut[j+m] = t2;                       
				  }
				};
				tbb::parallel_for(range, f, tbb::simple_partitioner());
				//f(range);
			}
*/		}
	}

	virtual void forwards_impl(
		size_t n,	const complex_t &wn,
		const complex_t *pIn,
		complex_t *pOut
	) const
	{
		assert(n>0);

		recurse(n,wn,pIn,1,pOut,1);
	}

	virtual void backwards_impl(
		size_t n,	const complex_t &wn,
		const complex_t *pIn,
		complex_t *pOut
	) const
	{
		complex_t reverse_wn=real_t(1)/wn;
		recurse(n, reverse_wn, pIn, 1, pOut, 1);

		real_t scale=real_t(1)/n;
		for(size_t i=0;i<n;i++){
			pOut[i]=pOut[i]*scale;
		}
	}

public:
	virtual std::string name() const
	{ return "hpce.abp.fast_fourier_transform_parfor"; }

	virtual bool is_quadratic() const
	{ return false; }
};

std::shared_ptr<fourier_transform> Create_fast_fourier_transform_parfor()
{
	return std::make_shared<fast_fourier_transform_parfor>();
}

}; // namespace abp14

}; // namespace hpce
