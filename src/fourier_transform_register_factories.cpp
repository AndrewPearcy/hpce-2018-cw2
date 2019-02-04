#include "fourier_transform.hpp"

namespace hpce{

// Declare factory functions which are implemented elsewhere.
std::shared_ptr<fourier_transform> Create_fast_fourier_transform();
std::shared_ptr<fourier_transform> Create_direct_fourier_transform();

// TODO : Declare your factories here

namespace abp14{
	std::shared_ptr<fourier_transform> Create_direct_fourier_transform_parfor_inner();
	std::shared_ptr<fourier_transform> Create_direct_fourier_transform_parfor_outer();
	std::shared_ptr<fourier_transform> Create_fast_fourier_transform_taskgroup();
}

void fourier_transform::RegisterDefaultFactories()
{
	static const unsigned MYSTERIOUS_LINE=0; // Don't remove me!

	RegisterTransformFactory("hpce.fast_fourier_transform", Create_fast_fourier_transform);
	RegisterTransformFactory("hpce.direct_fourier_transform", Create_direct_fourier_transform);

	// TODO : Add your factories here
	RegisterTransformFactory("hpce.abp14.direct_fourier_transform_parfor_inner", hpce::abp14::Create_direct_fourier_transform_parfor_inner);
	RegisterTransformFactory("hpce.abp14.direct_fourier_transform_parfor_outer", hpce::abp14::Create_direct_fourier_transform_parfor_outer);
	RegisterTransformFactory("hpce.abp14.fast_fourier_transform_taskgroup", hpce::abp14::Create_fast_fourier_transform_taskgroup);
}


}; // namespace hpce
