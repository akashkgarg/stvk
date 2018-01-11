//------------------------------------------------------------------------------
// tri_mesh.h
//------------------------------------------------------------------------------

#ifndef __LOOMAI_GEOMETRY_TRIMESH_H__
#define __LOOMAI_GEOMETRY_TRIMESH_H__

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------

#include <types.h>
#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>

//------------------------------------------------------------------------------
// class declarations
//------------------------------------------------------------------------------

struct TriMesh {

    //----------------------------------------------------------------------
    // typedefs
    //----------------------------------------------------------------------

    using  ptr       = std::shared_ptr<TriMesh>;
    using  const_ptr = std::shared_ptr<const TriMesh>;

    //----------------------------------------------------------------------
    // construction
    //----------------------------------------------------------------------

    TriMesh();
    TriMesh(const std::string &name);

    TriMesh(const TriMesh& other);
    TriMesh(const_ptr other);

    //----------------------------------------------------------------------
    // members
    //----------------------------------------------------------------------

    void Scale (double factor);
    void ComputeVertexNormals ();
    void FlipVertexNormals ();

    //----------------------------------------------------------------------
    // static methods.
    //----------------------------------------------------------------------

    static MatrixX3d ComputeVertexNormals(const MatrixX3u& faces,
                                          const MatrixX3d& verts);

    static void ComputeAdjacencyMatrix(const MatrixX3u &F,
                                       AdjacencyMatrix &A);

    static void ComputeEdges(const MatrixX3u &F, MatrixX2u &E);

    //----------------------------------------------------------------------
    // data members
    //----------------------------------------------------------------------

    std::string name;
    MatrixX3d  v;
    MatrixX3uc vc;
    MatrixX3d  n;
    MatrixX2d  vt;
    MatrixX2d  vt2; // 2nd uv set if applicable
    MatrixX3u  f;
    MatrixX3u  ft;
};

//------------------------------------------------------------------------------

#endif // __LOOMAI_GEOMETRY_TRIMESH_H__
