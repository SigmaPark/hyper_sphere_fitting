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
#include <string>
#include <iostream>
#include <chrono>
#include <map>
#include <cassert>


template<class T>
class Rms_and_stdev
{
public:
	Rms_and_stdev() noexcept : _sum(0), _sqr_sum(0), _nof_data(0){}

	auto reset() noexcept-> void{  _sum = _sqr_sum = 0;  _nof_data = 0;  }

	auto record(T const val) noexcept-> void
	{
		_sum += val,
		_sqr_sum += val*val,
		++_nof_data;
	}

	auto nof_data() const noexcept-> std::size_t{  return _nof_data;  }
	
	auto get() const noexcept-> sgm::Family<T, T>
	{
		assert(_nof_data != 0);

		auto const avg_of_sqrs = _sqr_sum / _nof_data;
		auto const avg = _sum / _nof_data;
		auto const variance = avg_of_sqrs - avg*avg;
		
		auto const 
			rms = static_cast<T>( std::sqrt(avg_of_sqrs) ),
			stdev = static_cast<T>( std::sqrt(variance) );
		
		return {rms, stdev};
	}

private:
	T _sum, _sqr_sum;
	std::size_t _nof_data;
};
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T>
class Test_Records
{
public:
	auto record_result(T const dc, T const dr, T const rmse) noexcept-> void
	{
		_dc.record(dc);
		_dr.record(dr);
		_rmse.record(rmse);
	}

	auto record_time(std::chrono::duration<T> const& time) noexcept-> void
	{
		using time_unit_t = std::chrono::microseconds;

		_time.record
		(	static_cast<T>
			(	std::chrono::duration_cast<time_unit_t>(time).count() 
			)
		);
	}

	auto record_as_failure() noexcept-> void
	{
		++_nof_failed;
	}

	auto report(std::string const& title) const-> void
	{
		auto const nof_success = _rmse.nof_data();

		if(nof_success == 0)
		{
			std::cerr << title << " has no record.\n";

			return;
		}

		auto const [rms_dc, sigma_dc] = _dc.get();
		auto const [rms_dr, sigma_dr] = _dr.get();
		auto const [rms_rmse, sigma_rmse] = _rmse.get();
		auto const [rms_time, sigma_time] = _time.get();

		assert(_rmse.nof_data() == _time.nof_data());


		std::cout 
		<<	title << std::endl
		<<	'\t' <<	nof_success << " pass out of " << nof_success + _nof_failed << std::endl
		<<	'\t' << "dc = " << rms_dc << " +- " << sigma_dc << std::endl
		<<	'\t' << "dr = " << rms_dr << " +- " << sigma_dr << std::endl
		<<	'\t' << "rmse = " << rms_rmse << " +- " << sigma_rmse << std::endl
		<<	'\t' << "time = " << rms_time << " +- " << sigma_time << " us\n";
	}


private:
	Rms_and_stdev<T> _dc{}, _dr{}, _rmse{}, _time{};
	std::size_t _nof_failed = 0;
};
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class TREC, class MEASURE, class FN, class...ARGS>
static auto Add_test_result(TREC& rec, MEASURE&& measure, FN&& fn, ARGS&&...args)-> void
{
	auto const start_tp = std::chrono::system_clock::now();
	auto const result = fn( sgm::Forward<ARGS>(args)... );
	auto const time = std::chrono::system_clock::now() - start_tp;

	if constexpr( sgm::is_Nullable<decltype(result)>::value )
		if(!result.has_value())
		{	
			rec.record_as_failure();

			return;
		}

	auto const [dc, dr, rmse] = +measure(result);

	rec.record_result(dc, dr, rmse);
	rec.record_time(time);
}


int main()
{
	using std::size_t;

	static auto constexpr Dim = size_t(3);
	
	using real_t = double;
	using Vec_t = s3d::Vector<real_t, Dim>;

	auto const ran_sphere
	=	[]()-> prac::test::Random_hyper_sphere<real_t, Dim>
		{
			auto constexpr max_position = real_t(500), sphere_radius = real_t(1'000);

			return {-max_position, max_position, sphere_radius, 0};
		}();
	
	auto constexpr nof_points_per_test = size_t(1'000);
	auto constexpr gaussian_noise_sigma = real_t(10);
	auto constexpr partial_rate = real_t(0.6);

	auto constexpr nof_tests = size_t(1000);
	auto constexpr ran_seed = std::uint64_t(777);
	
	std::mt19937_64 ran_engine{ran_seed};
	std::map< std::string, Test_Records<real_t> > fitting_method;

	for(auto count = nof_tests;  count-->0;)
	{
		prac::test::Test_Sphere_Fitting const testbed
		(	ran_sphere.get(ran_engine), nof_points_per_test, gaussian_noise_sigma, partial_rate
		, 	ran_engine
		);

		auto const measure_f
		=	sgm::Memfunc(testbed, &prac::test::Test_Sphere_Fitting<real_t, Dim>::dc_dr_rmse);

		auto const& sample_points = testbed.sample_points();


		::Add_test_result
		(	fitting_method["Least Square"], measure_f
		, 	SGM_1st_Class_Citizen(prac::Least_square_fit_sphere, <real_t, Dim>)
		,	sample_points
		);

		::Add_test_result
		(	fitting_method["Sumith YD"], measure_f
		, 	SGM_1st_Class_Citizen(prac::Sumith_YD_fit_sphere, <real_t, Dim>)
		,	sample_points
		);

		{
			size_t constexpr max_iteration = 200;

			::Add_test_result
			(	fitting_method["Eberly's iteration"], measure_f
			, 	SGM_1st_Class_Citizen(prac::Eberly_fit_sphere, <real_t, Dim>)
			,	sample_points, max_iteration
			);
		}
		{
			auto constexpr stable_ratio_diff = real_t(1e-4);
			auto constexpr max_iteration = size_t(20);
			auto const mu = real_t( std::exp(1) );

			auto initial_sphere 
			=	*reinterpret_cast< sgm::Array<real_t, Dim+1> const* >
				(	&testbed.original_sphere()
				);

			std::mt19937_64 ran_eng{0};
			std::normal_distribution<real_t> nd(0, initial_sphere[Dim]/10);

			for(auto& t : initial_sphere)
				t += nd(ran_eng);

			::Add_test_result
			(	fitting_method["Levenberg-Marquardt"], measure_f
			, 	SGM_1st_Class_Citizen(prac::Levenberg_Marquardt_fit_sphere, <real_t, Dim>)
			,	sample_points, initial_sphere
			,	real_t(1), mu, stable_ratio_diff, max_iteration
			);
		}
		{
			size_t constexpr nof_trials = 1000;
			real_t constexpr thres_val = 20;

			::Add_test_result
			(	fitting_method["MSAC"], measure_f
			, 	SGM_1st_Class_Citizen(prac::MSAC_fit_sphere, <real_t, Dim>)
			,	sample_points, nof_trials, thres_val
			);
		}
	}

	for(auto const& [method, record] : fitting_method)
		record.report(method);

    return 0;
}
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#