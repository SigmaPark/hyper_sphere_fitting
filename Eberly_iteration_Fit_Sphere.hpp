#pragma once
#ifndef _PRAC_EBERLY_FIT_SPHERE_
#define _PRAC_EBERLY_FIT_SPHERE_

#include "Sphere_Fitting_interface.hpp"
#include "SGM/Abbreviable/Nullable.hpp"

namespace prac
{

	template<class T, std::size_t D, class CON>
	static auto Eberly_fit_sphere(CON const& con, std::size_t const max_iteration)
	->	sgm::Nullable< sgm::Array<T, D+1> >;

}


template<class T, std::size_t D, class CON>
auto prac::Eberly_fit_sphere(CON const& con, std::size_t const max_iteration)
->	sgm::Nullable< sgm::Array<T, D+1> >
{
	using Vec_t = s3d::Vector<T, D>;

	auto const xbar
	= [&con]()-> Vec_t
	{
		Vec_t sx = decltype(sx)::Zero();
		std::size_t n = 0;

		for(auto const& pt : con)
			sx += sph_fit::To_Vector<T, D>(pt),
			++n;

		return sx / n;
	}();

	auto c = xbar;

	for(auto k = max_iteration;  k-->0;)
	{
		std::size_t n = 0;
		T a_norm_bar = 0;
		Vec_t a_hat_bar = decltype(a_hat_bar)::Zero();

		for(auto const& pt : con)
		{
			Vec_t const a = c - sph_fit::To_Vector<T, D>(pt);
			T const a_norm = a.norm();

			a_norm_bar += a_norm,
			a_hat_bar += a/a_norm,
			++n;
		}

		a_norm_bar /= n,
		a_hat_bar /= n;

		Vec_t const c_new = xbar + a_norm_bar*a_hat_bar;
		
		if( (c - c_new).sqr_norm() < T(1e-6) )
			return sph_fit::As_Array<T, D>(c, a_norm_bar);
		
		c = c_new;
	}

	return sgm::Null_t{};
}


#endif // #ifndef _PRAC_EBERLY_FIT_SPHERE_