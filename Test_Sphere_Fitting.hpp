/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_TEST_SPHERE_FITTING_
#define _PRAC_TEST_SPHERE_FITTING_

#include "S3D/Euclid/Euclid.hpp"
#include "SGM/Functor/Functor.hpp"
#include <cstdint>
#include <cassert>
#include <random>
#include <iostream>


namespace prac::test
{

    template<class T, std::size_t D>
    struct Hyper_Sphere;


	template<class T, std::size_t D>
	class Test_Sphere_Fitting;


    template<class T, std::size_t D>
    static auto Sample_points
    (	std::size_t const nof_points, Hyper_Sphere<T, D> const& answer_sphere
	,	T const gaussian_noise_sigma, T const partial_rate
	,	std::uint64_t const ran_seed 
    )->	sgm::Array< s3d::Vector<T, D> >;

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
	,	T const gaussian_noise_sigma, T const partial_rate, std::uint64_t const ran_seed
	);

	Test_Sphere_Fitting(Test_Sphere_Fitting const&) = delete;
	auto operator=(Test_Sphere_Fitting const&)-> Test_Sphere_Fitting& = delete;

	auto sample_points() const noexcept-> sgm::Array<Vec_t> const&;
	auto original_sphere() const noexcept-> Sphere_t const&;
	auto test(Sphere_t const& fit_sphere) const noexcept-> sgm::Array<T, 2>;
	auto test(sgm::Array<T, D+1> const& fit_sphere) const noexcept-> sgm::Array<T, 2>;

private:
	sgm::Array<Vec_t> _points;
	Sphere_t _answer_sphere;


	static auto _Generate_samples
	(	std::size_t const nof_points, Sphere_t const& answer_sphere
	,	T const gaussian_noise_sigma, T const partial_rate
	,	std::uint64_t const ran_seed
	)->	sgm::Array<Vec_t>;
};


template<class T, std::size_t D>
prac::test::Test_Sphere_Fitting<T, D>::Test_Sphere_Fitting
(	Sphere_t const& sphere, std::size_t const nof_points
, 	T const gaussian_noise_sigma, T const partial_rate, std::uint64_t const ran_seed
)
:	_points( _Generate_samples(nof_points, sphere, gaussian_noise_sigma, partial_rate, ran_seed) )
,	_answer_sphere(sphere)
{
	using std::cout, std::endl;

	cout 
	<<	"sphere dimension = " << D << endl
	<<	"max #of sample points = " << nof_points << endl
	<<	"gaussian noise sigma = " << gaussian_noise_sigma << endl
	<<	"sphere partial rate = " << 100*partial_rate << " %\n";

	cout << "original sphere center = [";

	for(std::size_t k = 0;  k < D;  cout << ", ",  ++k)
		cout << sphere.center(k);

	cout << "]\n";

	cout << "original sphere radius = " << sphere.radius << endl;	
}


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>::sample_points() const noexcept
->	sgm::Array<Vec_t> const&{  return _points;  }


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>::original_sphere() const noexcept
->	Sphere_t const&{  return _answer_sphere;  }


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>
::	test(Sphere_t const& fit_sphere) const noexcept-> sgm::Array<T, 2>
{
	return
	{	s3d::Distance(fit_sphere.center, original_sphere().center)
	,	fit_sphere.radius - original_sphere().radius
	};
}

template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>
::	test(sgm::Array<T, D+1> const& fit_sphere) const noexcept-> sgm::Array<T, 2>
{
	Vec_t c;

	for(std::size_t k = 0;  k < D;  ++k)
		c(k) = fit_sphere[k];

	return test(Sphere_t{c, fit_sphere[D]});
}


template<class T, std::size_t D>
auto prac::test::Test_Sphere_Fitting<T, D>::_Generate_samples
(	std::size_t const nof_points, Sphere_t const& answer_sphere
,	T const gaussian_noise_sigma, T const partial_rate
,	std::uint64_t const ran_seed
)->	sgm::Array<Vec_t>
{
	assert(nof_points >= D + 1);
	assert(gaussian_noise_sigma >= 0);
	assert(0 <= partial_rate && partial_rate <= 1);

	using UnitVec_t = s3d::UnitVec<T, D>;

	std::mt19937_64 ran_engine{ran_seed};

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