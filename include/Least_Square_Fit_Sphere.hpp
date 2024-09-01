/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_LEAST_SQUARE_FIT_SPHERE_
#define _PRAC_LEAST_SQUARE_FIT_SPHERE_

#include "Sphere_Fitting_interface.hpp"
#include "_Matrix_operations.hpp"


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
	=	[&con]
		{
			s3d::Matrix<T, D+1, D+1> M = decltype(M)::Zero();
			s3d::Vector<T, D+1> g = decltype(g)::Zero();

			for(auto const& p : con)
			{
				auto const& v = sph_fit::To_Vector<T, D>(p);
				auto const vv = v.sqr_norm();
				
				auto const a 
				=	[&v]() noexcept-> s3d::Vector<T, D+1>
					{
						s3d::Vector<T, D+1> res;

						res.head(D) = 2*v;
						res(D) = 1;

						return res;
					}();

				M += mat_op_helper::Self_dyadic_upper<T, D+1>(a);
				g += vv*a;
			}

			mat_op_helper::Copy_upper_to_lower<T, D+1>(M);

			return sgm::Make_Family(M, g);
		}();

	
	s3d::Vector<T, D+1> const x = M.inv()*g;
	s3d::Vector<T, D> const c = x.head(D);
	T const r = std::sqrt( x(D) + c.sqr_norm() );

	return sph_fit::As_Array<T, D>(c, r);
}


#endif // #ifndef _PRAC_LEAST_SQUARE_FIT_SPHERE_