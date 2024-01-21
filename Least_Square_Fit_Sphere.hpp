/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_LEAST_SQUARE_FIT_SPHERE_
#define _PRAC_LEAST_SQUARE_FIT_SPHERE_

#include "Sphere_Fitting_interface.hpp"


namespace prac
{

	template<class T, std::size_t D, class CON>
	static auto Least_square_fit_sphere(CON const& con)-> sgm::Array<T, D+1>;

}


template<class T, std::size_t D, class CON>
auto prac::Least_square_fit_sphere(CON const& con)-> sgm::Array<T, D+1>
{
	using std::size_t;

	auto const [M, g]
	= [&con]
	{
		s3d::Matrix<T, D + 1, D + 1> M = decltype(M)::Zero();
		s3d::Vector<T, D + 1> g = decltype(g)::Zero();

		for(auto const& p : con)
		{
			auto const& v = sph_fit::To_Vector<T, D>(p);
			auto const vv = v.sqr_norm();
	
			for(size_t i = 0;  i < D;  ++i)
			{
				g(i) += 2*v(i)*vv;
		
				for(size_t j = i;  j < D;  ++j)
					M(i, j) += 4*v(i)*v(j);
			}

			g(D) += vv;
		
			for(size_t i = 0;  i < D;  ++i)
				M(i, D) += 2*v(i);
			
			M(D, D) += 1;
		}
		
		for(size_t i = 0;  i < D + 1;  ++i)
			for(size_t j = 0;  j < i;  ++j)
				M(i, j) = M(j, i);

		return sgm::Make_Family(M, g);
	}();

	
	s3d::Vector<T, D + 1> const x = M.inv()*g;
	s3d::Vector<T, D> const c = x.head(D);
	T const r = std::sqrt( x(D) + c.sqr_norm() );

	return sph_fit::As_Array<T, D>(c, r);
}


#endif // #ifndef _PRAC_LEAST_SQUARE_FIT_SPHERE_