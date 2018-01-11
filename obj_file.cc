//------------------------------------------------------------------------------
// obj_file.cc
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cstdio>
#include <fstream>
#include <exception>
#include <vector>

#include "obj_file.h"
#include "tri_mesh.h"

//------------------------------------------------------------------------------
// class implementation
//------------------------------------------------------------------------------

#define MSG(c) (std::stringstream() << c).str()

//------------------------------------------------------------------------------
// Read methods.
//------------------------------------------------------------------------------

TriMesh::ptr
ObjFile::read(const std::string &filepath, double scale)
{
    auto mesh = std::make_shared<TriMesh>(filepath);
    read(filepath, mesh, scale);
    return mesh;
}

//------------------------------------------------------------------------------

void
ObjFile::read(const std::string &filepath, TriMesh::ptr mesh,
                double scale)
{
#if defined(ANDROID) && !defined(NO_JNI)
    // Get memory stream.
    AssetManager::DataType data;
    bool b = AssetManager::instance().readAsset(filepath, data);
    if (!b) {
        throw std::invalid_argument(MSG("Could not open obj file: " << filepath));
    }
    memstream ms(&data[0], data.size());
    read(ms, mesh, scale);
#else
    // Open file, and check for error
    std::ifstream file(filepath.c_str());
    if (!file.good()) {
        throw std::invalid_argument(MSG("Could not open obj file: " << filepath));
    }
    read(file, mesh, scale);
    file.close();
#endif // ANDROID
}

//------------------------------------------------------------------------------

void
ObjFile::read(std::istream &inputstream, TriMesh::ptr mesh, double scale)
{
    // File open was succesfull so clear outputs
    std::vector<std::vector<Scalar> > V;
    std::vector<std::vector<Scalar> > VT;
    std::vector<std::vector<Scalar> > VN;
    std::vector<std::vector<Index> >  F;
    std::vector<std::vector<Index> >  FT;
    std::vector<std::vector<Index> >  FN;

    const std::string v("v");
    const std::string vn("vn");
    const std::string vt("vt");
    const std::string f("f");
    const std::string comment("#");

    const int max_line_bytes = 2048;

    int line_no = 1;
    std::string readline;
    while (std::getline(inputstream, readline)) {
        const char *line = readline.c_str();
        char type[max_line_bytes];
        // Read first word containing type
        if(sscanf(line, "%s", type) == 1)
        {
            // Get pointer to rest of line right after type
            const char * l = &line[strlen(type)];
            if (type == v)
            {
                double x[4];
                int count =
                    sscanf(l,"%lf %lf %lf %lf\n",&x[0],&x[1],&x[2],&x[3]);
                if (count != 3 && count != 4) {
                    throw std::invalid_argument(
                        MSG("vertex on line " << line_no <<
                            " should have 3 or 4 coordinates"));
                }
                std::vector<Scalar> vertex(count);
                for (int i = 0;i<count;i++) {
                    vertex[i] = x[i] * scale;
                }
                V.push_back(vertex);
            } else if(type == vn) {
                double x[3];
                int count = sscanf(l,"%lf %lf %lf\n",&x[0],&x[1],&x[2]);
                if(count != 3) {
                    throw std::invalid_argument(
                        MSG("normal on line " << line_no <<
                            " should have 3 or 4 coordinates"));
                }
                std::vector<Scalar > normal(count);
                for (int i = 0;i<count;i++) {
                    normal[i] = x[i];
                }
                VN.push_back(normal);
            } else if(type == vt) {
                double x[3];
                int count = sscanf(l,"%lf %lf %lf\n",&x[0],&x[1],&x[2]);
                if (count != 2 && count != 3) {
                    throw std::invalid_argument(
                        MSG("texture coords on line  " << line_no <<
                            "%d should have 2 or 3 coordinates"));
                }
                std::vector<Scalar > tex(count);
                for (int i = 0;i<count;i++) {
                    tex[i] = x[i];
                }
                VT.push_back(tex);
            } else if(type == f) {
                const auto & shift = [&V](const int i)->int
                    {
                        return i<0 ? i+V.size() : i-1;
                    };
                const auto & shift_t = [&VT](const int i)->int
                    {
                        return i<0 ? i+VT.size() : i-1;
                    };
                const auto & shift_n = [&VN](const int i)->int
                    {
                        return i<0 ? i+VN.size() : i-1;
                    };

                std::vector<Index> f;
                std::vector<Index> ft;
                std::vector<Index> fn;
                // Read each "word" after type
                char word[max_line_bytes];
                int offset;
                while (sscanf(l,"%s%n",word,&offset) == 1)
                {
                    // adjust offset
                    l += offset;
                    // Process word
                    long int i,it,in;
                    if(sscanf(word,"%ld/%ld/%ld",&i,&it,&in) == 3) {
                        f.push_back(shift(i));
                        ft.push_back(shift_t(it));
                        fn.push_back(shift_n(in));
                    } else if(sscanf(word,"%ld/%ld",&i,&it) == 2) {
                        f.push_back(shift(i));
                        ft.push_back(shift_t(it));
                    } else if(sscanf(word,"%ld//%ld",&i,&in) == 2) {
                        f.push_back(shift(i));
                        fn.push_back(shift_n(in));
                    } else if(sscanf(word,"%ld",&i) == 1) {
                        f.push_back(shift(i));
                    } else {
                        throw std::invalid_argument(
                            MSG("face on line " << line_no <<
                                " has invalid element format"));
                    }
                }

                const auto& addFace = [&F,&FT,&FN,&line_no](
                                      const std::vector<Index>& f,
                                      const std::vector<Index>& ft,
                                      const std::vector<Index>& fn) {
                    if ((f.size() > 0 &&
                         fn.size() == 0 && ft.size() == 0) ||
                        (f.size() > 0 &&
                         fn.size() == f.size() && ft.size() == 0) ||
                        (f.size() > 0 &&
                         fn.size() == 0 && ft.size() == f.size()) ||
                        (f.size() > 0 &&
                         fn.size() == f.size() && ft.size() == f.size())) {
                        // No matter what add each type to lists so that
                        // lists are the correct lengths
                        F .push_back(f);
                        FT.push_back(ft);
                        FN.push_back(fn);
                    } else {
                        throw std::invalid_argument (
                            MSG("face on line " << line_no <<
                                " has invalid format"));
                    }
                };

                // Ensure only triangle or quad meshes.
                // Quad meshes will be triangulated.
                if (f.size() > 4) {
                    throw std::invalid_argument (
                        MSG("face on line " << line_no <<
                            " is not a triangle or quad."));
                } else if (f.size() == 4) {
                    std::vector<Index> t1_f;
                    std::vector<Index> t2_f;
                    t1_f.push_back(f[0]);
                    t1_f.push_back(f[1]);
                    t1_f.push_back(f[2]);
                    t2_f.push_back(f[0]);
                    t2_f.push_back(f[2]);
                    t2_f.push_back(f[3]);

                    std::vector<Index> t1_ft;
                    std::vector<Index> t2_ft;
                    if (ft.size() == 4) {
                        t1_ft.push_back(ft[0]);
                        t1_ft.push_back(ft[1]);
                        t1_ft.push_back(ft[2]);
                        t2_ft.push_back(ft[0]);
                        t2_ft.push_back(ft[2]);
                        t2_ft.push_back(ft[3]);
                    }

                    std::vector<Index> t1_fn;
                    std::vector<Index> t2_fn;
                    if (fn.size() == 4) {
                        t1_fn.push_back(fn[0]);
                        t1_fn.push_back(fn[1]);
                        t1_fn.push_back(fn[2]);
                        t2_fn.push_back(fn[0]);
                        t2_fn.push_back(fn[2]);
                        t2_fn.push_back(fn[3]);
                    }

                    addFace (t1_f, t1_ft, t1_fn);
                    addFace (t2_f, t2_ft, t2_fn);
                } else {
                    addFace (f, ft, fn);
                }
            } else if(strlen(type) >= 1 && (type[0] == '#' ||
                                            type[0] == 'g'  ||
                                            type[0] == 's'  ||
                                            strcmp("usemtl",type)==0 ||
                                            strcmp("mtllib",type)==0)) {
                //ignore comments or other shit
            } else {
                //ignore any other lines
            }
        } else {
            // ignore empty line
        }
        line_no++;
    }

    assert(F.size() == FN.size());
    assert(F.size() == FT.size());

    if (VN.size() > 0 && VN.size() > V.size()) {
        // When normals are present in the incoming OBJ, they may not
        //  match the indices of vertices (which can be shared between faces).
        // When this is the case and since we don't support multiple normals per
        //  vertex then we have to split the vertices so they match the normals
        // In the Face vertices
        std::vector<std::vector<Scalar> > V_split = VN;

        for (size_t faceI = 0; faceI < F.size(); faceI++) {
            for (size_t vtxI = 0; vtxI < F[faceI].size(); vtxI++) {
                V_split[FN[faceI][vtxI]] = V[F[faceI][vtxI]];
            }
        }

        V = V_split;
        F = FN;
    }

    ListToMatrix(V, mesh->v);
    ListToMatrix(VN, mesh->n);
    ListToMatrix(VT, mesh->vt);
    ListToMatrix(F, mesh->f);
    ListToMatrix(FT, mesh->ft);
}

//------------------------------------------------------------------------------

template <typename T, typename EigenMatrix>
void
ObjFile::ListToMatrix(const std::vector<std::vector<T> > &vals,
                          EigenMatrix &M)
{
    M.resize(vals.size(), M.cols());

    for (unsigned i = 0; i < M.rows(); ++i) {
        for (unsigned j = 0; j < M.cols(); ++j) {
            M(i, j) = vals[i][j];
        }
    }
}

//------------------------------------------------------------------------------
// Write methods.
//------------------------------------------------------------------------------

void
ObjFile::write(TriMesh::const_ptr mesh, const std::string &filepath)
{
    FILE * file = fopen(filepath.c_str(), "w");
    if (NULL == file) {
        fprintf(stderr,"IOError: %s could not be opened...\n",
                filepath.c_str());
        throw std::invalid_argument("Could not open obj file for writing.");
    }
    write(mesh, file);
}

//------------------------------------------------------------------------------

void
ObjFile::write(TriMesh::const_ptr mesh, FILE *file)
{
    // write vertices
    for (unsigned int i = 0; i < mesh->v.rows(); ++i) {
        fprintf(file, "v %0.17g %0.17g %0.17g\n",
                mesh->v(i,0), mesh->v(i, 1), mesh->v(i, 2));
    }

    // write normals, if mesh has them.
    if (mesh->n.rows() > 0) {
        for (unsigned i = 0; i < mesh->n.rows(); ++i) {
            fprintf(file, "vn %0.17g %0.17g %0.17g\n",
                    mesh->n(i, 0), mesh->n(i, 1), mesh->n(i, 2));
        }
    }

    // write texture coords.
    if (mesh->vt.rows() > 0) {
        for (unsigned i = 0; i < mesh->vt.rows(); ++i) {
            fprintf(file, "vt %0.17g %0.17g\n",
                    mesh->vt(i, 0), mesh->vt(i, 1));
        }
    }

    if (mesh->f.cols() != 3 || mesh->ft.cols() != 3) {
        throw std::invalid_argument("Mesh has non triangular faces");
    }

    // write faces
    for (unsigned i = 0; i < mesh->f.rows(); ++i) {
        const auto& shift = [](const int i)->int { return i + 1; };
        if (mesh->ft.rows() > 0) {
            fprintf(file, "f %u/%u %u/%u %u/%u\n",
                    shift(mesh->f(i, 0)), shift(mesh->ft(i, 0)),
                    shift(mesh->f(i, 1)), shift(mesh->ft(i, 1)),
                    shift(mesh->f(i, 2)), shift(mesh->ft(i, 2)));
        } else {
            fprintf(file, "f %u %u %u\n",
                    shift(mesh->f(i, 0)),
                    shift(mesh->f(i, 1)),
                    shift(mesh->f(i, 2)));
        }
    }

    fclose(file);
}

//------------------------------------------------------------------------------
// Explicit template instantiations.
//------------------------------------------------------------------------------

template
void ObjFile::ListToMatrix<Scalar, MatrixX3d>
(const std::vector<std::vector<Scalar> >&, MatrixX3d &);

template
void ObjFile::ListToMatrix<Scalar, MatrixX2d>
(const std::vector<std::vector<Scalar> >&, MatrixX2d &);

template
void ObjFile::ListToMatrix<Index, MatrixX3i>
(const std::vector<std::vector<Index> >&, MatrixX3i &);

//------------------------------------------------------------------------------
