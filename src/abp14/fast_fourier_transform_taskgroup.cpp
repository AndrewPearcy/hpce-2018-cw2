#include "fourier_transform.hpp"

#include <cmath>
#include <cassert>
#include "tbb/task_group.h"

namespace hpce
{
namespace abp14
{
class fast_fourier_transform_taskgroup
	: public fourier_transform
{
protected:

	char *v = getenv("HPCE_FTT_RECURSION_K");
	size_t grain = v==NULL ? 16 : atoi(v);
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
		tbb::task_group group;
		assert(n>0);

		if (n == 1){
			pOut[0] = pIn[0];
		}else if (n == 2){
			pOut[0] = pIn[0]+pIn[sIn];
			pOut[sOut] = pIn[0]-pIn[sIn];
		}else{

			size_t m = n/2;
			if(n > grain) {
				auto left=[&](){ recurse(m,wn*wn,pIn,2*sIn,pOut,sOut); };
				auto right=[&]() { recurse(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut); };
				group.run(left);
				group.run(right);
				group.wait();


			} else {
				recurse(m,wn*wn,pIn,2*sIn,pOut,sOut);
				recurse(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut);
				
			}
			complex_t w=complex_t(1, 0);

			for (size_t j=0;j<m;j++){
			  complex_t t1 = w*pOut[m+j];
			  complex_t t2 = pOut[j]-t1;
			  pOut[j] = pOut[j]+t1;                 /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
			  pOut[j+m] = t2;                          /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
			  w = w*wn;
			}
		}
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
	{ return "hpce.abp.fast_fourier_transform_taskgroup"; }

	virtual bool is_quadratic() const
	{ return false; }
};

std::shared_ptr<fourier_transform> Create_fast_fourier_transform_taskgroup()
{
	return std::make_shared<fast_fourier_transform_taskgroup>();
}

}; // namespace abp14
}; // namespace hpce