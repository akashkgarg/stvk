//------------------------------------------------------------------------------
// types.h
//------------------------------------------------------------------------------

#ifndef __LOOMAI_TYPES_H__
#define __LOOMAI_TYPES_H__

//-----------------------------------------------------------------------------
// includes
//-----------------------------------------------------------------------------

#include <memory>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

//-----------------------------------------------------------------------------
// typedefs
//-----------------------------------------------------------------------------

//- signed integer -------------------------------------------------------------

using i8           = int8_t;
using i16          = int16_t;
using i32          = int32_t;
using i64          = int64_t;

//- unsigned integer -----------------------------------------------------------

using u8           = uint8_t;
using u16          = uint16_t;
using u32          = uint32_t;
using u64          = uint64_t;

//- float ----------------------------------------------------------------------

using f32          = float;
using f64          = double;
using Scalar       = f64;
using Real         = f64;

//- offsets --------------------------------------------------------------------

using Offset       = u32;
using Index        = u32;

//- vectors --------------------------------------------------------------------

using Vec4f        = Eigen::Vector4f;
using Vec4d        = Eigen::Vector4d;
using Vec4i        = Eigen::Vector4i;
using Vec4uc       = Eigen::Matrix<unsigned char, 4, 1>;
using Vec4         = Vec4d;

using Vec3f        = Eigen::Vector3f;
using Vec3d        = Eigen::Vector3d;
using Vec3i        = Eigen::Vector3i;
using Vec3uc       = Eigen::Matrix<unsigned char, 3, 1>;
using Vec3u        = Eigen::Matrix<unsigned, 3, 1>;
using Vec3         = Vec3d;

using Vec2f        = Eigen::Vector2f;
using Vec2d        = Eigen::Vector2d;
using Vec2i        = Eigen::Vector2i;
using Vec2uc       = Eigen::Matrix<unsigned char, 2, 1>;
using Vec2         = Vec2d;

using VecXf        = Eigen::VectorXf;
using VecXd        = Eigen::VectorXd;
using VecXi        = Eigen::VectorXi;
using VecXb        = Eigen::Matrix<bool, Eigen::Dynamic, 1>;
using VecXuc       = Eigen::Matrix<unsigned char, Eigen::Dynamic, 1>;
using VecXu        = Eigen::Matrix<unsigned, Eigen::Dynamic, 1>;
using VecX         = VecXd;

//- transformations ------------------------------------------------------------

using Quaterniond  = Eigen::Quaterniond;
using Quaternionf  = Eigen::Quaternionf;
using Quaternion   = Quaterniond;

using AngleAxisd   = Eigen::AngleAxisd;
using AngleAxisf   = Eigen::AngleAxisf;

//- matrices -------------------------------------------------------------------

using Matrix4f     = Eigen::Matrix4f;
using Matrix4d     = Eigen::Matrix4d;
using Matrix4      = Matrix4d;

using MatrixXd     = Eigen::MatrixXd;
using MatrixXf     = Eigen::MatrixXf;
using MatrixXu     = Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixXuc    = Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>;
using MatrixX      = MatrixXd;

using Matrix3f     = Eigen::Matrix3f;
using Matrix3d     = Eigen::Matrix3d;
using Matrix3      = Matrix3d;

using MatrixX3f    = Eigen::MatrixX3f;
using MatrixX3d    = Eigen::MatrixX3d;
using MatrixX3i    = Eigen::MatrixX3i;
using MatrixX3     = Eigen::MatrixX3d;

using MatrixX2f    = Eigen::MatrixX2f;
using MatrixX2d    = Eigen::MatrixX2d;
using MatrixX3i    = Eigen::MatrixX3i;
using MatrixX3u    = Eigen::Matrix<unsigned, Eigen::Dynamic, 3>;
using MatrixX3uc   = Eigen::Matrix<unsigned char, Eigen::Dynamic, 3>;
using MatrixX2i    = Eigen::MatrixX2i;
using MatrixX2u    = Eigen::Matrix<unsigned, Eigen::Dynamic, 2>;
using MatrixX2     = Eigen::MatrixX2d;

using SparseMatrix     = Eigen::SparseMatrix<Scalar>;
using AdjacencyMatrix  = Eigen::SparseMatrix<unsigned int>;
using NNZTriplet       = Eigen::Triplet<Scalar>;

using IdxList = std::vector<Index>;
using IdxSet =  std::unordered_set<Index>;

//------------------------------------------------------------------------------
// Useful conversion macros.
//------------------------------------------------------------------------------

#define ARR_TO_VEC3(p) Vec3(p[0], p[1], p[2])

//-----------------------------------------------------------------------------

#endif // __LOOMAI_TYPES_H__
