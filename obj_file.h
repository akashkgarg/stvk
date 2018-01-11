//------------------------------------------------------------------------------
// obj_file.h
//------------------------------------------------------------------------------

#ifndef __LOOMAI_IO_OBJFILE_H__
#define __LOOMAI_IO_OBJFILE_H__

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------

#include <types.h>
#include <tri_mesh.h>
#include <string>
#include <Eigen/Core>
#include <Eigen/Geometry>

//------------------------------------------------------------------------------
// class declarations
//------------------------------------------------------------------------------

class ObjFile {
public:
    static TriMesh::ptr read(const std::string &filepath, double scale = 1.0);

    static void read(const std::string &filepath,
                     TriMesh::ptr mesh,
                     double scale = 1.0);

    static void write(TriMesh::const_ptr mesh,
                      const std::string &filepath);

private:
    static void read(std::istream &inputstream, TriMesh::ptr mesh,
                     double scale);

    static void write(TriMesh::const_ptr mesh, FILE *file);

    template <typename T, typename EigenMatrix>
        static void ListToMatrix(const std::vector<std::vector<T> > &vals,
                                 EigenMatrix &M);
};

//------------------------------------------------------------------------------

#endif // __LOOMAI_IO_OBJFILE_H__
