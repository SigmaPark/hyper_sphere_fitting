/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_TEST_SPHERE_FITTING_
#define _PRAC_TEST_SPHERE_FITTING_

#include "Euclid/Euclid.hpp"
#include "Functor/Functor.hpp"
#include <cstdint>
#include <cassert>
#include <random>
#include <iostream>


namespace prac::test
{

    template<class T, std::size_t D>
    struct Hyper_Sphere;


	template<class T, std::size_t D>
	class Random_hyper_sphere;

	template<class T, std::size_t D>
	class Test_Sphere_Fitting;

}
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T, std::size_t D>
struct prac::test::Hyper_Sphere
{
	using elem_t = T;
	static auto constexpr Dim = D;

	s3d::Vector<T, D> center;
	T radius;
};
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T, std::size_t D>
class prac::test::Random_hyper_sphere
{
public:
	Random_hyper_sphere(T const min_x, T const max_x, T const mu_r, T const sigma_r) noexcept;

	auto get(std::mt19937_64& ran_eng) const noexcept-> Hyper_Sphere<T, D>;

private:
	T _min_x, _max_x, _mu_r, _sigma_r;
};


template<class T, std::size_t D>
prac::test::Random_hyper_sphere<T, D>
::	Random_hyper_sphere(T const min_x, T const max_x, T const mu_r, T const sigma_r) noexcept
:	_min_x(min_x), _max_x(max_x), _mu_r(mu_r), _sigma_r(sigma_r)
{}


template<class T, std::size_t D>
auto prac::test::Random_hyper_sphere<T, D>
::	get(std::mt19937_64& ran_eng) const noexcept-> Hyper_Sphere<T, D>
{
	auto const center
	= [min_x = _min_x, max_x = _max_x, &ran_eng]()-> s3d::Vector<T, D>
	{
		std::uniform_real_distribution<T> urd(min_x, max_x);
		s3d::Vector<T, D> res;

		for( std::size_t k = 0;  k < D;  res(k++) = urd(ran_eng) );

		return res;
	}();

	auto const radius
	=	_sigma_r == 0 ? _mu_r
	:	std::normal_distribution<T>(_mu_r, _sigma_r)(ran_eng);

	return {center, radius};
}
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T, std::size_t D>
class prac::test::Test_Sphere_Fitting
{
public:
	static auto constexpr Dim = D;

	using elem_t = T;
	using Vec_t = s3d::Vector<T, D>;
	using Sphere_t = Hyper_Sphere<T, D>;

	Test_Sphere_Fitting
	(	Sphere_t const& sphere, std::size_t const nof_points
	,	T const gaussian_noise_sigma, T const partial_rate, std::mt19937_64& ran_engine
	);

	auto reset
	(	Sphere_t const& sphere, std::size_t const nof_points
	,	T const gaussian_noise_sigma, T const partial_rate, std::mt19937_64& ran_engine
	)->	void;

	Test_Sphere_Fitting(Test_Sphere_Fitting const&) = delete;
	auto operator=(Test_Sphere_Fitting const&)-> Test_Sphere_Fitting& = delete;

	auto sample_points() const noexcept-> sgm::Array<Vec_t> const&;
	auto original_sphere() const noexcept-> Sphere_t const&;
	
	auto dc_dr_rmse(sgm::Array<T, D+1> const& fit_sphere) const noexcept-> sgm::Array<T, 3>;

private:
	sgm::Array<Vec_t> _points;
	Sphere_t _answer_sphere;


	static auto _Generate_samples
	(	std::size_t const nof_points, Sphere_t const& answer_sphere
	,	T const gaussian_noise_sigma, T const partial_rate
	,	std::mt19937_64& ran_engine
	)->	sgm::Array<Vec_t>;
};


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>
::	reset
(	Sphere_t const& sphere, std::size_t const nof_points
, 	T const gaussian_noise_sigma, T const partial_rate, std::mt19937_64& ran_engine
)->	void
{
	_points = _Generate_samples(nof_points, sphere, gaussian_noise_sigma, partial_rate, ran_engine);
	_answer_sphere = sphere;
}


template<class T, std::size_t D>
prac::test::Test_Sphere_Fitting<T, D>::Test_Sphere_Fitting
(	Sphere_t const& sphere, std::size_t const nof_points
, 	T const gaussian_noise_sigma, T const partial_rate, std::mt19937_64& ran_engine
) :	_points(), _answer_sphere()
{
	reset(sphere, nof_points, gaussian_noise_sigma, partial_rate, ran_engine);  
}


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>::sample_points() const noexcept
->	sgm::Array<Vec_t> const&{  return _points;  }


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>::original_sphere() const noexcept
->	Sphere_t const&{  return _answer_sphere;  }


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>
::	dc_dr_rmse(sgm::Array<T, D+1> const& fit_sphere) const noexcept-> sgm::Array<T, 3>
{
	auto const fit_sph
	= [&fit_sphere]()-> Sphere_t
	{
		Vec_t c;

		for(std::size_t k = 0;  k < D;  ++k)
			c(k) = fit_sphere[k];

		return {c, fit_sphere[D]};
	}();

	auto const rmse_dist
	= [&fit_sph, &points = sample_points()]
	{
		auto sqr_f = [](T const t){  return t*t;  };

		T sqrsum = 0;

		for(auto const& pt : points)
			sqrsum += sqr_f( (pt - fit_sph.center).norm() - fit_sph.radius );

		return std::sqrt(sqrsum/points.size());
	}();

	T const 
		dc = s3d::Distance(fit_sph.center, original_sphere().center),
		dr = fit_sph.radius - original_sphere().radius;

	return {dc, dr, rmse_dist};
}


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>::_Generate_samples
(	std::size_t const nof_points, Sphere_t const& answer_sphere
,	T const gaussian_noise_sigma, T const partial_rate
,	std::mt19937_64& ran_engine
)->	sgm::Array<Vec_t>
{
	assert(nof_points >= D + 1);
	assert(gaussian_noise_sigma >= 0);
	assert(0 <= partial_rate && partial_rate <= 1);

	using UnitVec_t = s3d::UnitVec<T, D>;

	auto uniform_random_direction_f
	= [&rng = ran_engine]()-> UnitVec_t
	{
		std::uniform_real_distribution<T> urd(-1, 1);
		Vec_t pt;

		do
			for( std::size_t k = 0;  k < D;  pt(k++) = urd(rng) );
		while( pt.sqr_norm() < T(1e-6) );

		return pt;
	};

	auto uniform_sphere_f
	= [&urdirec_f = uniform_random_direction_f](std::size_t const n)-> sgm::Array<Vec_t>
	{
		sgm::Array<Vec_t> res(n);

		while(res.size() != n)
			res >> urdirec_f().vec();

		return res;
	};

	auto points_on_sphere_f
	= [&rng = ran_engine, &urdirec_f = uniform_random_direction_f, sph = answer_sphere]
	(sgm::Array<Vec_t> points, T const sigma)-> decltype(points)
	{
		std::normal_distribution<T> nd(0, sigma);

		for(auto& pt : points)
			pt = sph.center + sph.radius*pt + nd(rng)*urdirec_f();

		return points;
	};

	auto partial_sphere_f
	= [cut_direction = uniform_random_direction_f(), sph = answer_sphere]
	(sgm::Array<Vec_t> points, T const rate)-> decltype(points)
	{
		if(rate >= 1)
			return points;

		s3d::Plane<T, D> const cut_plane
		(	sph.center + sph.radius*(2*rate - 1)*cut_direction
		,	cut_direction
		);

		return sgm::fp::Filter_f(  points, SGM_LAMBDA( cut_plane.signed_dist_to(_0) < 0 )  );
	};

	return
	partial_sphere_f
	(	points_on_sphere_f
		(	uniform_sphere_f(nof_points)
		,	gaussian_noise_sigma
		)
	,	partial_rate
	);
}

#endif // #ifndef _PRAC_TEST_SPHERE_FITTING_