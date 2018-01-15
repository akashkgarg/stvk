
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <array>
#include <tri_mesh.h>
#include <random>


using CClothModel = TriMesh;
typedef Eigen::Matrix<double, 3, 2> Matrix3x2d;
typedef Eigen::Matrix<double, 3, 9> Matrix3x9d;
typedef Eigen::Matrix<double, 1, 3> RowVec3d;
typedef Eigen::Matrix<double, 9, 1> Vec9d;
typedef Eigen::Matrix<double, 9, 9> Mat9d;
using CSparseMat = MatrixXd;
using CVector3d = Vec3d;
using Matrix3x3 = Matrix3d;

const double _mu = 500;
const double _poisson = 1.0; 

struct PrecomputedTriangleElement
{
	double A;
	Eigen::Matrix3d DmInverse;
	Matrix3x2d dFdu[3][3];
};

double GetAreaUV(const CClothModel *cloth,
                 const Vec3u &f)
{
    auto a = cloth->vt.row(f[0]);
    auto b = cloth->vt.row(f[1]);
    auto c = cloth->vt.row(f[2]);

    Eigen::Matrix2d D;
    D.col(0) = (b - a);
    D.col(1) = (c - a);

    return fabs(D.determinant() / 2.0);
}

class StVenantKirchoffForce {
protected:
	std::vector<PrecomputedTriangleElement> m_elements;

public:
	void precompute(const CClothModel* pCloth)
	{
		m_elements.resize(pCloth->f.rows());

		for (int i = 0; i < pCloth->f.rows(); i++)
		{
            auto f = pCloth->f.row(i);

            // material space
            Eigen::Vector2d x[3];
			for (int j = 0; j < 3; j++)
			{
				double u, v;
                u = pCloth->vt.row(f[j])[0];
                v = pCloth->vt.row(f[j])[1];
				x[j] << u, v;
			}

            Vec3d d0, d1, d2;
            d0 << x[1] - x[0], 0.0;
            d1 << x[2] - x[0], 0.0;
            d2 = d0.cross(d1);
            double dn = d2.norm();
            Matrix3d Dm3;
            Dm3.col(0) = d0;
            Dm3.col(1) = d1;
            Dm3.col(2) = d2/dn;
            double a = 0.5 * dn;

            Matrix3d invDm = Dm3.inverse();

            m_elements[i].A = a;
            m_elements[i].DmInverse = invDm;
		}
	}

    Matrix3d deformationGradient(CClothModel *cloth, const size_t fidx)
    {
        auto f = cloth->f.row(fidx);
        Vec3d v0 = cloth->v.row(f[0]);
        Vec3d v1 = cloth->v.row(f[1]);
        Vec3d v2 = cloth->v.row(f[2]);
        Vec3d n  = ((v1 - v0).cross(v2 - v0)).normalized();

        // Shape matrix.
        Matrix3d Ds;
        Ds.col(0) = v1 - v0;
        Ds.col(1) = v2 - v0;
        Ds.col(2) = n;

        // Deformation gradient.
        Matrix3d F = Ds * m_elements[fidx].DmInverse;

        return F;
    }

    Matrix3d greenStrain(const Matrix3d &F)
    {
        // Green strain tensor
        Matrix3d G = (F.transpose() * F - Matrix3d::Identity()) * 0.5;
        return G;
    }

    double computeElementEnergy(CClothModel *cloth, const size_t fidx)
	{
        Matrix3d F = deformationGradient(cloth, fidx);

        // Green strain tensor
        Matrix3d G = greenStrain(F);

        // sigma
        Matrix3d sigma = _mu * (1.0 - _poisson) * G + Matrix3d::Identity() * (_mu * _poisson * G.trace());

        auto inner = [](const Matrix3d &a, const Matrix3d &b) -> double {
            double r=0;
            for (int j=0; j<3; j++)
                for (int i=0; i<3; i++)
                    r+=a.col(j)[i]*b.col(j)[i];
            return r;
        };

        double W = m_elements[fidx].A * 0.5 * inner(sigma, G);

        return W;
	}

    Matrix3x9d kronecker(const RowVec3d &A, const Matrix3d &B)
    {
        Matrix3x9d C;
        for (int i = 0; i < 1; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        C(i*3+k,j*3+l) = A(i,j)*B(k,l);
        return C;
        
    }

    void computeForce(CClothModel *cloth, const size_t fidx,
                      Matrix3d &forces, CSparseMat &J)
    {
        auto face = cloth->f.row(fidx);
        const Vec3d x[3] = {cloth->v.row(face[0]), cloth->v.row(face[1]), cloth->v.row(face[2])};

        Matrix3d F = deformationGradient(cloth, fidx);
        Matrix3d G = greenStrain(F);
        Matrix3d Y = m_elements[fidx].DmInverse;
        Matrix3d D;
        D.row(0) = -Y.row(0) - Y.row(1);
        D.row(1) = Y.row(0);
        D.row(2) = Y.row(1);

        Matrix3x9d DD[3] = { kronecker(D.col(0).transpose(), Matrix3d::Identity()),
                             kronecker(D.col(1).transpose(), Matrix3d::Identity()),
                             kronecker(D.col(2).transpose(), Matrix3d::Identity()) };

        Vec9d X;
        for (int i=0; i<9; i++)
            X[i] = x[i/3][i%3];
        Vec3d f[3] = {DD[0]*X,DD[1]*X,DD[2]*X};

        Vec9d grad_f = Vec9d::Zero();
        Mat9d hess_f = Mat9d::Zero();
        double gf = -m_elements[fidx].A * _mu; // negate for forces/jacobian
        //Mat3x3 Gc = max(G,Mat3x3(0)); // posdef
        Matrix3d Gc = G.array().unaryExpr([](double x) { return std::max(x, 0.0); }); // posdef
        for (int i=0; i<3; i++) {
            for(int j=0; j<=i; j++) {
                Vec9d dG = 0.5 * (DD[i].transpose()*f[j] + DD[j].transpose()*f[i]);
                Mat9d dG2 = 0.5 * (DD[i].transpose()*DD[j] + DD[j].transpose()*DD[i]);
                Vec9d d_trace_c = 0.5*dG; // posdef

                if (i==j)
                    grad_f += gf * (1.0 - _poisson) * G(i,j) * dG +
                        gf * _poisson * G.trace() * dG;
                else
                    grad_f += 2.0 * gf * (1.0 - _poisson) * G(i,j) * dG;
                
                if (i==j)
                    hess_f += gf * (1.0 - _poisson) * (dG*dG.transpose() + Gc(i,j) * dG2) +
                        gf * _poisson * (Gc.trace() * dG2 + (d_trace_c*dG.transpose()));
                else
                    hess_f += 2.0 * gf * (1.0 - _poisson) * ((dG*dG.transpose()) +\
                                                             Mat9d::Zero()); // posdef
            }
        }
        std::cout << "---------------Forces-----------------------------" << std::endl;
        std::cout << grad_f.transpose() << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;

        std::cout << "---------------Hessian-----------------------------" << std::endl;
        std::cout << hess_f.transpose() << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
    }
};
