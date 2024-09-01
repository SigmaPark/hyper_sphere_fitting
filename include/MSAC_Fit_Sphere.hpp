/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_MSAC_FIT_SPHERE_
#define _PRAC_MSAC_FIT_SPHERE_

#include <cassert>
#include <random>
#include <set>
#include "Sphere_Fitting_interface.hpp"


namespace prac
{

	template<class T, std::size_t D, class CON>
	static auto MSAC_fit_sphere(CON const& con, std::size_t const trials, T const threshold)
	->	sgm::Array<T, D+1>;

}


namespace prac::_msac_impl
{

	template<class T, std::size_t D>
	class Random_sphere_generation;

}
//========//========//========//========//=======#//========//========//========//========//=======#


template<class T, std::size_t D>
class prac::_msac_impl::Random_sphere_generation
{
public:
	using Vec_t = s3d::Vector<T, D>;
	
	Random_sphere_generation
	(	sgm::Array<Vec_t const*>&& ptrs
	,	T const threshold, std::uint64_t const ran_seed = 0
	);

	auto operator()() noexcept-> sgm::Family<Vec_t, T, T>;

private:
	sgm::Array<Vec_t const*> _ptrs;
	std::mt19937_64 _ran_eng;
	T _threshold;

	auto _sphere(std::set<std::size_t> const& indices) const;
	auto _score(Vec_t const& c, T const r) const-> T;
};


template<class T, std::size_t D>
prac::_msac_impl::Random_sphere_generation<T, D>::Random_sphere_generation
(	sgm::Array<Vec_t const*>&& ptrs
, 	T const threshold, std::uint64_t const ran_seed
) :	_ptrs( sgm::Move(ptrs) ), _ran_eng{ran_seed}, _threshold(threshold)
{}


template<class T, std::size_t D>
auto prac::_msac_impl::Random_sphere_generation<T, D>
::	_sphere(std::set<std::size_t> const& indices) const
{
	s3d::Matrix<T, D+1, D+1> A;
	s3d::Vector<T, D+1> b;
	auto itr = indices.cbegin();

	for(std::size_t k = 0;  k < D+1;  ++k)
	{
		auto const idx = *itr++;
		Vec_t const& vk = *_ptrs[idx];

		A.row(k).head(D) = 2*vk.transpose();
		A.row(k)(D) = 1;
		b(k) = vk.sqr_norm();
	}

	s3d::Vector<T, D+1> const x = A.inv()*b;
	Vec_t const c = x.head(D);
	T const r = std::sqrt( x(D) + c.sqr_norm() );

	return sgm::Make_Family(c, r);
}


template<class T, std::size_t D>
auto prac::_msac_impl::Random_sphere_generation<T, D>::_score(Vec_t const& c, T const r) const-> T
{
	auto const sqr_thres = std::pow(_threshold, 2);
	T score = 0;

	for(auto const p : _ptrs)
		if(  T const ee = std::pow( (*p - c).norm() - r, 2 );  ee < sqr_thres  )
			score += sqr_thres - ee;

	return score;
}


template<class T, std::size_t D>
auto prac::_msac_impl::Random_sphere_generation<T, D>
::	operator()() noexcept-> sgm::Family<Vec_t, T, T>
{
	auto const indices
	= [&ran_eng = _ran_eng, max_num =_ptrs.size() - 1]
	{
		std::uniform_int_distribution<std::size_t> uid(0, max_num);
		std::set<std::size_t> idx_set;

		while(idx_set.size() != D+1)
			idx_set.emplace( uid(ran_eng) );

		return idx_set;
	}();

	auto const [sph_c, sph_r] = _sphere(indices);
	auto const score = _score(sph_c, sph_r);

	return {sph_c, sph_r, score};
}
//--------//--------//--------//--------//-------#//--------//--------//--------//--------//-------#


template<class T, std::size_t D, class CON>
auto prac::MSAC_fit_sphere(CON const& con, std::size_t const trials, T const threshold)
->	sgm::Array<T, D+1>
{
	using Vec_t = typename _msac_impl::Random_sphere_generation<T, D>::Vec_t;

	auto ran_sph_gen
	= [threshold, &con]()-> _msac_impl::Random_sphere_generation<T, D>
	{
		sgm::Array<Vec_t const*> sample_ptrs( sgm::Size(con) );

		for(auto const& pt : con)
			sample_ptrs >> &sph_fit::To_Vector<T, D>(pt);

		assert(sample_ptrs.size() > D+1);

		return {sgm::Move(sample_ptrs), threshold};
	}();

	Vec_t c;
	T r;
	T s = 0;

	for(auto d = trials;  d-->0;)
		if(auto const [nc, nr, ns] = ran_sph_gen();  ns > s)
			c = nc,  r = nr,  s = ns;

	return sph_fit::As_Array<T, D>(c, r);
}


#endif // #ifndef _PRAC_MSAC_FIT_SPHERE_