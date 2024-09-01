/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_SPHERE_FITTING_INTERFACE_
#define _PRAC_SPHERE_FITTING_INTERFACE_

#include "Array/Array.hpp"
#include "Euclid/Euclid.hpp"


namespace prac::sph_fit
{

	template<class T, std::size_t D, class Q>
	inline auto To_Vector(Q const& q) noexcept-> s3d::Vector<T, D> const&;


	template<class T, std::size_t D>
	inline auto As_Array(s3d::Vector<T, D> const& c, T const r) noexcept
	->	sgm::Array<T, D+1>;

}


template<class T, std::size_t D, class Q>
auto prac::sph_fit::To_Vector(Q const& q) noexcept-> s3d::Vector<T, D> const&
{
	return *reinterpret_cast< s3d::Vector<T, D> const* >( sgm::Address_of(q) );
}


template<class T, std::size_t D>
auto prac::sph_fit::As_Array(s3d::Vector<T, D> const& c, T const r) noexcept
->	sgm::Array<T, D+1>
{
	sgm::Array<T, D + 1> res;

	for(std::size_t i = 0;  i < D;  ++i)
		res[i] = c(i);

	res[D] = r;
		
	return res;	
}

#endif // #ifndef _PRAC_SPHERE_FITTING_INTERFACE_