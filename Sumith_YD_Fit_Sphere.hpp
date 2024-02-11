/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_SUMITH_YD_FIT_SPHERE_
#define _PRAC_SUMITH_YD_FIT_SPHERE_

#include "Sphere_Fitting_interface.hpp"
#include "_Matrix_operations.hpp"


namespace prac
{

	template<class T, std::size_t D, class CON>
	static auto Sumith_YD_fit_sphere(CON const& con)-> sgm::Array<T, D+1>;

}


template<class T, std::size_t D, class CON>
auto prac::Sumith_YD_fit_sphere(CON const& con)-> sgm::Array<T, D+1>
{
	using std::size_t;

	auto diag_f
	= [](s3d::Matrix<T, D, D> const& M)-> s3d::Vector<T, D>
	{
		s3d::Vector<T, D> res;

		for(size_t i = 0;  i < D;  ++i)
			res(i) = M(i, i);

		return res;
	};

	auto const [N, S1, S2, S3]
	= [diag_f, &con]
	{
		size_t N = 0;
		s3d::Vector<T, D> S1 = decltype(S1)::Zero();
		s3d::Matrix<T, D, D> S2 = decltype(S2)::Zero(), S3 = S2;

		for(auto const& p : con)
		{
			auto const& v = sph_fit::To_Vector<T, D>(p);
			auto const uvdv = mat_op_helper::Self_dyadic_upper<T, D>(v);

			++N;
			S1 += v;
			S2 += uvdv;
			S3 += diag_f(uvdv).dyadic(v);
		}

		mat_op_helper::Copy_upper_to_lower<T, D>(S2);

		return sgm::Make_Family(N, S1, S2, S3);
	}();

	auto const tr_S2
	= [&S2=S2]
	{
		T res = 0;

		for(size_t k = 0;  k < D;  ++k)
			res += S2(k, k);

		return res;
	}();

	auto const S3_colsum
	= [&S3=S3]
	{
		s3d::Vector<T, D> res = decltype(res)::Zero();

		for(size_t i = 0;  i < D;  ++i)
			for(size_t j = 0;  j < D;  ++j)
				res(j) += S3(i, j);

		return res;
	}();

	s3d::Matrix<T, D, D> const A = 2*( mat_op_helper::Self_dyadic<T, D>(S1) - N*S2 );
	s3d::Vector<T, D> const b = tr_S2*S1 - N*S3_colsum;
	
	s3d::Vector<T, D> const c = A.inv()*b;
	T const r = std::sqrt(  ( tr_S2 - 2*c.dot(S1) )/N + c.sqr_norm()  );

	return sph_fit::As_Array<T, D>(c, r);	
}


#endif // #ifndef _PRAC_SUMITH_YD_FIT_SPHERE_