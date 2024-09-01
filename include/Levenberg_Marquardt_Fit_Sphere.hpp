/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_LEVENBERG_MARQUARDT_FIT_SPHERE_
#define _PRAC_LEVENBERG_MARQUARDT_FIT_SPHERE_

#include <limits>
#include "Sphere_Fitting_interface.hpp"
#include "Abbreviable/Nullable.hpp"


namespace prac
{

	template<class T, std::size_t D, class CON>
	static auto Levenberg_Marquardt_fit_sphere
	(	CON const& con, sgm::Array<T, D+1> const& x0, T const lambda0, T const mu
	,	T const stable_ratio_diff, std::size_t const max_iteration
	)->	sgm::Nullable< sgm::Array<T, D+1> >;

	template<class T, std::size_t D, class CON>
	class Levenberg_Marquardt_fit_sphere_impl;

}
//========//========//========//========//=======#//========//========//========//========//=======#


template<class T, std::size_t D, class CON>
class prac::Levenberg_Marquardt_fit_sphere_impl
{
private: 
	using _point_t = s3d::Vector<T, D>;
	using _vec_t = s3d::Vector<T, D+1>;
	using _mat_t = s3d::Matrix<T, D+1, D+1>;

public:
	static auto Calc
	(	CON const& con, sgm::Array<T, D+1> const& x0, T const lambda0, T const mu
	,	T const stable_ratio_diff, std::size_t const max_iteration
	)->	sgm::Nullable< sgm::Array<T, D+1> >;

private:
	static auto _Diag_mat(_mat_t const& M) noexcept-> _mat_t;

	static auto _err(_point_t const& v, _point_t const& c, T const r) noexcept-> T;
	static auto _del_err(_point_t const& v, _point_t const& c) noexcept-> _vec_t;

	static auto _Tau_and_s
	(	_point_t const& c, T const r, CON const& con, T const lambda
	) 	noexcept-> s3d::Family<_vec_t, T>;
};


template<class T, std::size_t D, class CON>
auto prac::Levenberg_Marquardt_fit_sphere_impl<T, D, CON>
::	Calc
(	CON const& con, sgm::Array<T, D+1> const& x0, T const lambda0, T const mu
,	T const stable_ratio_diff, std::size_t const max_iteration
)->	sgm::Nullable< sgm::Array<T, D+1> >
{
	for
	(	auto [x, s, L, k] 
		= 	sgm::Make_Family
			( 	*reinterpret_cast<_vec_t const*>(&x0)
			, 	std::numeric_limits<T>::max() 
			,	lambda0
			,	max_iteration
			)
	;  	k-->0
	;
	)
	{
		_point_t const& c = *reinterpret_cast<_point_t const*>(&x);
		T const& r = reinterpret_cast<T const*>(&x)[D];
		auto const [tau, s_new] = _Tau_and_s(c, r, con, L);
		T const d = 1 - s_new/s;

		x -= tau;
		s = s_new;

		if(0 <= d && d < stable_ratio_diff)
			return sph_fit::As_Array<T, D>(c, r);
	
		L = d > 0 ? L/mu : L*mu;
	}

	return sgm::Null_t{};			
}


template<class T, std::size_t D, class CON>
auto prac::Levenberg_Marquardt_fit_sphere_impl<T, D, CON>
::	_Diag_mat(_mat_t const& M) noexcept-> _mat_t
{
	_mat_t res = _mat_t::Zero();

	for(std::size_t i = 0;  i < D+1;  ++i)
		res(i, i) = M(i, i);

	return res;		
}


template<class T, std::size_t D, class CON>
auto prac::Levenberg_Marquardt_fit_sphere_impl<T, D, CON>
::	_err(_point_t const& v, _point_t const& c, T const r) noexcept-> T
{
	return r - (v - c).norm();  
}


template<class T, std::size_t D, class CON>
auto prac::Levenberg_Marquardt_fit_sphere_impl<T, D, CON>
::	_del_err(_point_t const& v, _point_t const& c) noexcept-> _vec_t
{
	_vec_t res;

	res.head(D) = (v - c).normalized();
	res(D) = 1;

	return res;
}


template<class T, std::size_t D, class CON>
auto prac::Levenberg_Marquardt_fit_sphere_impl<T, D, CON>
::	_Tau_and_s
(	_point_t const& c, T const r, CON const& con, T const lambda
) 	noexcept-> s3d::Family<_vec_t, T>
{
	_mat_t JJt = _mat_t::Zero();
	_vec_t Jf = _vec_t::Zero();
	T s = 0;

	for(auto const& pt : con)
	{
		auto const& v = sph_fit::To_Vector<T, D>(pt);
		auto const e = _err(v, c, r);
		auto const del_e = _del_err(v, c);

		s += e*e;
		Jf += e*del_e;
		JJt += del_e.dyadic(del_e);
	}

	return 
	{	( JJt + lambda*_Diag_mat(JJt) ).inv()*Jf
	, 	s
	};
}
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T, std::size_t D, class CON>
auto prac::Levenberg_Marquardt_fit_sphere
(	CON const& con, sgm::Array<T, D+1> const& x0, T const lambda0, T const mu
,	T const stable_ratio_diff, std::size_t const max_iteration
)->	sgm::Nullable< sgm::Array<T, D+1> >
{
	return 
	Levenberg_Marquardt_fit_sphere_impl<T, D, CON>::Calc
	(	con, x0, lambda0, mu, stable_ratio_diff, max_iteration
	);
}


#endif // #ifndef _PRAC_LEVENBERG_MARQUARDT_FIT_SPHERE_