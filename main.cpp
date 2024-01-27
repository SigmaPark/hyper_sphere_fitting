/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#include "Test_Sphere_Fitting.hpp"
#include "Least_Square_Fit_Sphere.hpp"
#include "Sumith_YD_Fit_Sphere.hpp"
#include "Eberly_iteration_Fit_Sphere.hpp"
#include "MSAC_Fit_Sphere.hpp"
#include "Levenberg_Marquardt_Fit_Sphere.hpp"
#include <iostream>


static auto Testbed()
{
	using real_t = double;
	auto constexpr dimension = std::size_t(3);

	using Vec_t = s3d::Vector<real_t, dimension>;

	auto constexpr nof_points = std::size_t(10'000);
	auto constexpr gaussian_noise_sigma = real_t(10);
	auto constexpr partial_rate = real_t(0.6);
	auto constexpr ran_seed = std::uint64_t(777);

	static_assert(nof_points >= dimension + 1);
	static_assert(gaussian_noise_sigma >= 0);
	static_assert(0 <= partial_rate && partial_rate <= 1);

	auto random_vec_f
	= [dimension, rng = std::mt19937_64{ran_seed}]() mutable-> Vec_t
	{
		std::uniform_real_distribution<real_t> urd(-500, 500);
		Vec_t res;

		for( std::size_t k = 0;  k < dimension;  res(k++) = urd(rng) );

		return res;
	};

	prac::test::Hyper_Sphere<real_t, dimension> const answer_sphere{random_vec_f(), real_t(1'000)};

	return
	prac::test::Test_Sphere_Fitting<real_t, dimension>
	(	answer_sphere, nof_points, gaussian_noise_sigma, partial_rate, ran_seed
	);
}


int main()
{
	prac::test::Test_Sphere_Fitting const testbed = ::Testbed();

	using real_t = typename decltype(testbed)::elem_t;
	static auto constexpr dim = decltype(testbed)::Dim;

	auto const& sample_points = testbed.sample_points();

	auto report_f
	= [&testbed](sgm::Array<real_t, dim+1> const& result_sphere, char const* title)
	{
		auto const [dc, dr] = +testbed.test(result_sphere);

		auto const rmse
		= 	[	&c = *reinterpret_cast< s3d::Vector<real_t, dim> const* >(result_sphere.data())
			,	r = result_sphere.back()
			, 	&samples = testbed.sample_points()
			]
		{
			real_t s = 0;

			for(auto const& v :  samples)
				s += std::pow( (v - c).norm() - r, 2 );

			return std::sqrt(s / samples.size());
		}();

		std::cout 
		<<	title 
		<<	" : dc = " << dc << ", dr = " << dr << ", rmse = " << rmse << std::endl;
	};

	std::cout << std::endl;

	report_f
	( 	*reinterpret_cast< sgm::Array<real_t, dim+1> const* >(&testbed.original_sphere())
	, 	"Original Sphere" 
	);
	
	report_f( prac::Least_square_fit_sphere<real_t, dim>(sample_points), "Least Square" );
	report_f( prac::Sumith_YD_fit_sphere<real_t, dim>(sample_points), "Sumith YD" );

	if
	(	auto const result = prac::Eberly_fit_sphere<real_t, dim>(sample_points, 200)
	;	result.has_value()
	)
		report_f(result.v(), "Eberly's iteration");
	else
		std::cerr << "Eberly's iteration fitting Failed." << std::endl;
	
	
	report_f( prac::MSAC_fit_sphere<real_t, dim>( sample_points, 1'000, real_t(20) ), "MSAC" );

	{
		auto constexpr stable_ratio_diff = real_t(1e-4);
		auto constexpr max_iteration = std::size_t(20);
		auto const mu = real_t( std::exp(1) );

		auto initial_sphere 
		=	*reinterpret_cast< sgm::Array<real_t, dim+1> const* >(&testbed.original_sphere());

		std::mt19937_64 ran_eng{0};
		std::normal_distribution<real_t> nd(0, initial_sphere[dim]/10);

		for(auto& t : initial_sphere)
			t += nd(ran_eng);

		auto const result
		=	prac::Levenberg_Marquardt_fit_sphere<real_t, dim>
			(	sample_points
			,	initial_sphere
			,	real_t(1), mu, stable_ratio_diff, max_iteration
			);

		if(result.has_value())
			report_f(result, "Levenberg-Marquardt");
		else
			std::cerr << "Levenberg-Marquardt fitting Failed." << std::endl;
	}

    return 0;
}