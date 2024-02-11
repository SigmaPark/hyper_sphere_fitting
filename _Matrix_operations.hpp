/*  SPDX-FileCopyrightText: (c) 2024 Jin-Eon Park <greengb@naver.com> <sigma@gm.gist.ac.kr>
*   SPDX-License-Identifier: MIT License
*/
//========//========//========//========//=======#//========//========//========//========//=======#


#pragma once
#ifndef _PRAC_MATRIX_OPERATIONS_
#define _PRAC_MATRIX_OPERATIONS_

#include "S3D/Hamilton/Hamilton.hpp"


namespace prac::mat_op_helper
{

	template<class T, std::size_t D>
	static auto Self_dyadic(s3d::Vector<T, D> const& v) noexcept-> s3d::Matrix<T, D, D>;	

	template<class T, std::size_t D>
	static auto Self_dyadic_upper(s3d::Vector<T, D> const& v) noexcept-> s3d::Matrix<T, D, D>;

	template<class T, std::size_t D>
	static auto Copy_upper_to_lower(s3d::Matrix<T, D, D>& M) noexcept-> void;

}


template<class T, std::size_t D>
auto prac::mat_op_helper::Self_dyadic(s3d::Vector<T, D> const& v) noexcept-> s3d::Matrix<T, D, D>
{
	auto M = Self_dyadic_upper<T, D>(v);

	Copy_upper_to_lower<T, D>(M);

	return M;
}


template<class T, std::size_t D>
auto prac::mat_op_helper::Self_dyadic_upper(s3d::Vector<T, D> const& v) noexcept
-> 	s3d::Matrix<T, D, D>
{
	s3d::Matrix<T, D, D> res = decltype(res)::Zero();

	for(std::size_t i = 0;  i < D;  ++i)
		for(std::size_t j = i;  j < D;  ++j)
			res(i, j) = v(i)*v(j);

	return res;
}


template<class T, std::size_t D>
auto prac::mat_op_helper::Copy_upper_to_lower(s3d::Matrix<T, D, D>& M) noexcept-> void
{
	for(std::size_t i = 0;  i < D;  ++i)
		for(std::size_t j = 0;  j < i;  ++j)
			M(i, j) = M(j, i);
}


#endif // #ifndef _PRAC_MATRIX_OPERATIONS_