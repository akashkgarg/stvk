//------------------------------------------------------------------------------
// tri_mesh.cc
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// includes
//------------------------------------------------------------------------------

#include <tri_mesh.h>

//------------------------------------------------------------------------------
// class implementation
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

TriMesh::TriMesh() :
    name("")
{
}

//------------------------------------------------------------------------------

TriMesh::TriMesh(const std::string &name) :
    name(name)
{
}

//------------------------------------------------------------------------------

TriMesh::TriMesh(const TriMesh& other) :
    name(other.name), v(other.v), n(other.n), vt(other.vt),
    f(other.f), ft(other.ft)
{
}

TriMesh::TriMesh(const_ptr other) :
    name(other->name), v(other->v), n(other->n), vt(other->vt),
    f(other->f), ft(other->ft)
{
}

//------------------------------------------------------------------------------

void
TriMesh::Scale (double factor)
{
    if (factor != 1.0) {
        v *= factor;
    }
}

//------------------------------------------------------------------------------

MatrixX3d
TriMesh::ComputeVertexNormals(const MatrixX3u& f, const MatrixX3d& v)
{
    using namespace std;

    const int num_v = v.rows();
    const int num_f = f.rows();

    // Resize for output
    MatrixX3d n(num_v, 3);
    n.setZero();

    MatrixX3d W(num_f, 3);
    W.setConstant(1.); // constant weight. probably should do something smarter.

    // compute per face normals first.
    MatrixX3d FN;
    FN.resize(num_f, 3);
    for (int fi = 0; fi < num_f; ++fi) {
        auto p0 = v.row( f(fi,0) );
        auto p1 = v.row( f(fi,1) );
        auto p2 = v.row( f(fi,2) );
        auto n0 = (p1 - p0).cross(p2 - p0);
        auto n1 = (p2 - p1).cross(p0 - p1);
        auto n2 = (p0 - p2).cross(p1 - p2);
        // careful sum
        for(int d = 0; d < 3; ++d) {
            // This is a little _silly_ in terms of complexity, but its recursive
            // implementation is clean looking...
            const function<double(double,double,double)> sum3 =
                [&sum3](double a, double b, double c)->double {
                if(fabs(c)>fabs(a)) {
                    return sum3(c,b,a);
                }
                // c < a
                if(fabs(c)>fabs(b)) {
                    return sum3(a,c,b);
                }
                // c < a, c < b
                if(fabs(b)>fabs(a)) {
                    return sum3(b,a,c);
                }
                return (a+b)+c;
            };
            FN(fi,d) = sum3(n0(d),n1(d),n2(d));
        }
        // sum better not be sure, or else NaN
        FN.row(fi) /= FN.row(fi).norm();
    }

    // loop over faces
    for(int i = 0; i < num_f; i++) {
        // throw normal at each corner
        for(int j = 0; j < 3; j++) {
            n.row( f(i,j) ) += W(i,j) * FN.row(i);
        }
    }

    // take average via normalization
    n.rowwise().normalize();

    return n;
}

//------------------------------------------------------------------------------

void
TriMesh::ComputeVertexNormals()
{
    n = ComputeVertexNormals (f, v);
}

//------------------------------------------------------------------------------

void
TriMesh::FlipVertexNormals()
{
    for (int i = 0; i < n.rows(); ++i) {
        n(i,0) *= -1.0;
        n(i,1) *= -1.0;
        n(i,2) *= -1.0;
    }
}

//------------------------------------------------------------------------------

void
TriMesh::ComputeAdjacencyMatrix(const MatrixX3u &F,
                                AdjacencyMatrix &A)
{
    // Type store to store the non-zero values in our adjacency matrix.
    using Triplet = Eigen::Triplet<unsigned int>;

    // the list of non-zeros.
    std::vector<Triplet> nnz;

    // approximate about 2x edges than faces.
    nnz.reserve(F.size()*2);

    for (unsigned i = 0; i < F.rows(); ++i) {
        for (unsigned j = 0; j < F.cols(); ++j) {
            auto from_idx = F(i, j);
            auto to_idx   = F(i, (j+1) % F.cols());
            nnz.push_back(Triplet(from_idx, to_idx, 1));
            nnz.push_back(Triplet(to_idx, from_idx, 1));
        }
    }

    auto n = F.maxCoeff()+1;
    A.resize(n, n);
    A.setZero();

    // predict proper number of nonzeros. For tri meshes, assume valence 6.
    A.reserve(6*n);

    // Always just one for any value in nnz.
    A.setFromTriplets(nnz.begin(), nnz.end(),
                      [](const AdjacencyMatrix::Scalar&,
                         const AdjacencyMatrix::Scalar&) {return 1;});
}

//------------------------------------------------------------------------------

void
TriMesh::ComputeEdges(const MatrixX3u &F, MatrixX2u &E)
{
    AdjacencyMatrix A;
    ComputeAdjacencyMatrix(F, A);

    E.resize(A.nonZeros()/2, 2);

    // Iterate non-zero entries.
    int e_row = 0;
    using InnerIterator = AdjacencyMatrix::InnerIterator;
    for (unsigned i = 0; i < A.outerSize(); ++i) {
        for (InnerIterator it(A, i); it; ++it) {
            // only add half the edges.
            if (it.row() < it.col()) {
                E(e_row, 0) = it.row();
                E(e_row, 1) = it.col();
                ++e_row;
            }
        }
    }
}

//------------------------------------------------------------------------------
