
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
using CSparseMat = MatrixXd;
using CVector3d = Vec3d;
using Matrix3x3 = Matrix3d;

const double _mu = 500;

struct PrecomputedTriangleElement
{
	double A;
	Eigen::Matrix2d DmInverse;
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

			Eigen::Vector2d x[3];
			for (int j = 0; j < 3; j++)
			{
				double u, v;
                u = pCloth->vt.row(f[j])[0];
                v = pCloth->vt.row(f[j])[1];
				x[j] << u, v;
			}

			Eigen::Matrix2d Dm;
			Dm.col(0) = x[1] - x[0];
			Dm.col(1) = x[2] - x[0];
			m_elements[i].DmInverse = Dm.inverse();

			//m_elements[i].A = 0.5 * std::sqrt(std::fabs(Dm.determinant()));
			m_elements[i].A = GetAreaUV(pCloth, f);

			// Derivative of deformation gradient F w.r.t. to each displacement u
			Eigen::Matrix2d& DmI = m_elements[i].DmInverse;
			Eigen::MatrixXd dFdu = Eigen::MatrixXd(6, 9);
			dFdu << -DmI(0, 0) - DmI(1, 0), 0.0, 0.0, DmI(0, 0), 0.0, 0.0, DmI(1, 0), 0.0, 0.0,
					-DmI(0, 1) - DmI(1, 1), 0.0, 0.0, DmI(0, 1), 0.0, 0.0, DmI(1, 1), 0.0, 0.0,
					0.0, -DmI(0, 0) - DmI(1, 0), 0.0, 0.0, DmI(0, 0), 0.0, 0.0, DmI(1, 0), 0.0,
					0.0, -DmI(0, 1) - DmI(1, 1), 0.0, 0.0, DmI(0, 1), 0.0, 0.0, DmI(1, 1), 0.0,
					0.0, 0.0, -DmI(0, 0) - DmI(1, 0), 0.0, 0.0, DmI(0, 0), 0.0, 0.0, DmI(1, 0),
					0.0, 0.0, -DmI(0, 1) - DmI(1, 1), 0.0, 0.0, DmI(0, 1), 0.0, 0.0, DmI(1, 1);

			// Store as a 3x2 matrix for each vertex and each coordinate
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					m_elements[i].dFdu[j][k] = dFdu.block<2, 3>(2 * j, 3 * k).transpose();
		}
	}

	double computeElementEnergy(const size_t idx,
                                const Eigen::Vector3d& x1,
                                const Eigen::Vector3d& x2,
                                const Eigen::Vector3d& x3)
	{
		// Shape coordinate basis
		Matrix3x2d Ds;
		Ds.col(0) = x2 - x1;
		Ds.col(1) = x3 - x1;

		// Deformation gradient F = Ds * Dm^-1
		Matrix3x2d F = Ds * m_elements[idx].DmInverse;

		Eigen::Matrix2d E = 0.5 * (F.transpose() * F - Eigen::Matrix2d::Identity());
		//E = clampMatrixEigenvalues(E);

		double W = m_elements[idx].A * _mu * (E.transpose() * E).trace();
		return W;
	}

	template<typename Func>
	Vec3d computeFDForce(const size_t idx, const double fx,
                         const std::array<Eigen::Vector3d, 3>& x,
                         Func energy)
	{
		const double eps = 1.0e-6;

		std::array<Eigen::Vector3d, 3> x_eps;
		for (size_t i = 0; i < 3; i++)
			x_eps[i] = x[i];

		Vec3d force;
		for (size_t k = 0; k < 3; k++)
		{
			x_eps[idx][k] += eps;
			force[k] = (energy(x_eps[0], x_eps[1], x_eps[2]) - fx) / eps; // finite-diff force
			x_eps[idx][k] = x[idx][k]; // Reset the vertex
		}

        // negate. force is negative grad of energy
		return -force;
	}

	template<typename Func>
	Matrix3d computeFDHessian(const size_t idx1, const size_t idx2,
                              const double fx,
                              const std::array<Eigen::Vector3d, 3>& x,
                              Func energy)
	{
		const double eps = 1.0e-6;

		std::array<Eigen::Vector3d, 3> x_eps;
		for (size_t i = 0; i < 3; i++)
			x_eps[i] = x[i];

		Matrix3d H;
		for (size_t k = 0; k < 3; k++)
		{
			for (size_t l = 0; l < 3; l++)
			{
                // diagonal entries
				if (idx1 == idx2 && k == l)
				{
					x_eps[idx1][k] += eps;
					const double f1 = energy(x_eps[0], x_eps[1], x_eps[2]);

					x_eps[idx1][k] = x[idx1][k] - eps;
					const double f2 = energy(x_eps[0], x_eps[1], x_eps[2]);

					x_eps[idx1][k] = x[idx1][k];
                    // FIXME: storage order of underlying matrix matters here
                    // when using [] operator.
					H.data()[3*k+l] = (f1 - 2.0 * fx + f2) / (eps * eps);
				}
				else
				{
					x_eps[idx1][k] += eps;
					x_eps[idx2][l] += eps;
					const double f1 = energy(x_eps[0], x_eps[1], x_eps[2]);

					x_eps[idx2][l] -= 2.0 * eps;
					const double f2 = energy(x_eps[0], x_eps[1], x_eps[2]);

					x_eps[idx1][k] -= 2.0 * eps;
					const double f4 = energy(x_eps[0], x_eps[1], x_eps[2]);

					x_eps[idx2][l] += 2.0 * eps;
					const double f3 = energy(x_eps[0], x_eps[1], x_eps[2]);

					x_eps[idx1][k] = x[idx1][k];
					x_eps[idx2][l] = x[idx2][l];

                    // FIXME: storage order of underlying matrix matters here
                    // when using [] operator.
					H.data()[3*k+l] = (f1 - f2 - f3 + f4) / (4.0 * eps * eps);
				}
			}
		}

		return H;
	}

    void gradEnergy(CClothModel *pCloth, Matrix3d &grad)
    {
        // ONLY A SINGLE ELEMENT IN OUR EXPERIMENTS.
        int i = 0;
        grad.setZero();

        //CFace3d* f = pCloth->GetFace(i);
        auto f = pCloth->f.row(i);

        std::array<Eigen::Vector3d, 3> x;
        for (int j = 0; j < 3; j++) {
            //x[j] << f->v(j)->m_pPos->x(), f->v(j)->m_pPos->y(), f->v(j)->m_pPos->z();
            auto pos = pCloth->v.row(f[j]);
            x[j] << pos[0], pos[1], pos[2];
        }

        // Shape coordinate basis
        Matrix3x2d Ds;
        for (size_t j = 0; j < 2; ++j)
            Ds.col(j) = x[j + 1] - x[0];

        // Deformation gradient F = Ds * Dm^-1
        Matrix3x2d F = Ds * m_elements[i].DmInverse;

        // 1st Piola stress tensor, S = k * E, where E = 0.5 * (F^T F - I) (Green-Lagrange strain)
        // Remove negative eigenvalues so that it can always recover triangle shape and is always P.D.
        Eigen::Matrix2d S = _mu * (F.transpose() * F - Eigen::Matrix2d::Identity());
        //S = clampMatrixEigenvalues(S);

        // 2nd Piola Stress tensor, P = F S, distribute to 2nd and 3rd node to get forces
        Matrix3x2d P2 = m_elements[i].A * F * S * m_elements[i].DmInverse.transpose();

        CVector3d f1(P2.col(0).data());
        CVector3d f2(P2.col(1).data());

        // Forces should sum to zero, so first vertex force is -(f_1 + f_2)
        grad.col(0) -= (f1 + f2);
        grad.col(1) += f1;
        grad.col(2) += f2;
    }

    void gradEnergyFD(CClothModel *pCloth, Matrix3d &grad)
    {
        double dx = 1.e-6;
        grad.setZero();

        auto v_old = pCloth->v;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                // compute central difference.

                // fij_forward
                pCloth->v(i,j) += dx;
                auto fij_forward = computeElementEnergy(0,
                                                        pCloth->v.row(0),
                                                        pCloth->v.row(1),
                                                        pCloth->v.row(2));
                pCloth->v = v_old;

                // fij_backward
                pCloth->v(i,j) -= dx;
                auto fij_backward = computeElementEnergy(0,
                                                         pCloth->v.row(0),
                                                         pCloth->v.row(1),
                                                         pCloth->v.row(2));
                pCloth->v = v_old;

                grad(i, j) = (fij_forward - fij_backward) / (2*dx);
            }
        }

        pCloth->v = v_old;
    }

    void hessianFD(CClothModel *pCloth, CSparseMat *H)
    {
        double dx = 1.e-6;
        H->setZero();

        // compute using finite differences.
        auto v_old = pCloth->v;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    Matrix3d grad;

                    // compute central difference.

                    // x_jk + dx
                    pCloth->v(j, k) += dx;
                    gradEnergy(pCloth, grad);
                    auto fi_forward = grad.col(i);
                    pCloth->v = v_old;

                    // x_jk - dx
                    pCloth->v(j, k) -= dx;
                    gradEnergy(pCloth, grad);
                    auto fi_backward = grad.col(i);
                    pCloth->v = v_old;

                    H->block<3, 1>(i*3, j*3 + k) = (fi_forward - fi_backward) / (2*dx);
                }
            }
        }

        // restore.
        pCloth->v = v_old;
    }

    void hessianAnalytical(CClothModel *pCloth, CSparseMat *dFdx)
    {
        // ONLY A SINGLE ELEMENT
        int i = 0;

        //CFace3d* f = pCloth->GetFace(i);
        auto f = pCloth->f.row(i);

        std::array<Eigen::Vector3d, 3> x;
        for (int j = 0; j < 3; j++) {
            //x[j] << f->v(j)->m_pPos->x(), f->v(j)->m_pPos->y(), f->v(j)->m_pPos->z();
            auto pos = pCloth->v.row(f[j]);
            x[j] << pos[0], pos[1], pos[2];
        }

        // Shape coordinate basis
        Matrix3x2d Ds;
        for (size_t j = 0; j < 2; ++j)
            Ds.col(j) = x[j + 1] - x[0];

        // Deformation gradient F = Ds * Dm^-1
        Matrix3x2d F = Ds * m_elements[i].DmInverse;

        // 1st Piola stress tensor, S = k * E, where E = 0.5 * (F^T F - I) (Green-Lagrange strain)
        // Remove negative eigenvalues so that it can always recover triangle shape and is always P.D.
        Eigen::Matrix2d S = _mu * (F.transpose() * F - Eigen::Matrix2d::Identity());
        //S = clampMatrixEigenvalues(S);

        // Take derivative of each force f_j with respect to each vertex k
        for (int j = 0; j < 3; j++)
        {
            // j-th force is A * F * S * Bm_j, where Bm_j is appropriate column(s) from Dm^-1
            Eigen::Vector2d Bm_j;
            if (j == 0)
                Bm_j = (-m_elements[i].DmInverse.row(0) - m_elements[i].DmInverse.row(1)).transpose();
            else
                Bm_j = m_elements[i].DmInverse.row(j-1).transpose();

            for (int k = j; k < 3; k++)
            {
                // Each column of 3x3 stiffness matrix is derivative of force f_j w.r.t. vertex k's l'th coordinate
                Eigen::Matrix3d K_jk;
                for (int l = 0; l < 3; l++)
                {
                    // Get dFdu for k-th vertex, l-th coordinate
                    Matrix3x2d& dFdu = m_elements[i].dFdu[k][l];
                    K_jk.col(l) = (dFdu * S + _mu * F * (F.transpose() * dFdu + dFdu.transpose() * F)) * Bm_j * m_elements[i].A;
                }

                Matrix3x3 K;
                K = K_jk;
                //K.SetCol(CVector3d(K_jk.col(0).data()), CVector3d(K_jk.col(1).data()), CVector3d(K_jk.col(2).data()));
                //dFdx->AddBlockMatrix(f->v(j)->GetIndex(), f->v(k)->GetIndex(), K);
                dFdx->block(f[j], f[k], 3, 3) = K;
            }
        }
    }
    
    void computeForces(CClothModel* pCloth, CSparseMat* dFdx,
                       Matrix3d &out_forces,
                       Matrix3d &fd_forces,
                       Matrix3d &F_approx)
    {
        // Compute forces.
        gradEnergy(pCloth, out_forces);
        // negate the grad, since forces = -\grad E
        out_forces *= -1;

        // negate the grad, since forces = -\grad E
        gradEnergyFD(pCloth, fd_forces);
        fd_forces *= -1;

        // Compute hessian
        //hessianFD(pCloth, dFdx);
        hessianAnalytical(pCloth, dFdx);
        // negate the hessian, since jacobian = -\grad^2 E
        dFdx->operator*=(-1);

        // Validate hessian using F(x) = F(x + dx) - dF. dF = dFdx * dx
        const double eps = 1.e-6;
        // make dx some random small perturbation
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        MatrixXd dx(9, 1);
        for (int i = 0; i < 9; ++i) {
            dx(i, 0) = eps * dist(mt);
        }

        auto dF = (*dFdx) * dx;

        // compute forces at x + dx
        auto v_old = pCloth->v;
        for (int i = 0; i < 3; ++i) {
            pCloth->v.row(i) += dx.block<3,1>(i*3, 0).transpose();
        }
        Matrix3d Fxdx;
        gradEnergy(pCloth, Fxdx);
        Fxdx *= -1; // negate for forces

        for (int i = 0; i < 3; ++i) {
            F_approx.col(i) = Fxdx.col(i) - dF.block<3, 1>(i*3, 0);
        }

        // restore
        pCloth->v = v_old;
    }

#if 0
    void computeForces(CClothModel* pCloth, CSparseMat* dFdx,
                       Matrix3d &out_forces,
                       Matrix3d &fd_forces, CSparseMat *fd_dFdx)
	{
		std::cout << "Computing forces" << std::endl;
		for (int i = 0; i < pCloth->f.rows(); i++)
		{
			//CFace3d* f = pCloth->GetFace(i);
            auto f = pCloth->f.row(i);

			std::array<Eigen::Vector3d, 3> x;
			for (int j = 0; j < 3; j++) {
				//x[j] << f->v(j)->m_pPos->x(), f->v(j)->m_pPos->y(), f->v(j)->m_pPos->z();
                auto pos = pCloth->v.row(f[j]);
				x[j] << pos[0], pos[1], pos[2];
            }

			// Shape coordinate basis
			Matrix3x2d Ds;
			for (size_t j = 0; j < 2; ++j)
				Ds.col(j) = x[j + 1] - x[0];

			// Deformation gradient F = Ds * Dm^-1
			Matrix3x2d F = Ds * m_elements[i].DmInverse;

			// 1st Piola stress tensor, S = k * E, where E = 0.5 * (F^T F - I) (Green-Lagrange strain)
			// Remove negative eigenvalues so that it can always recover triangle shape and is always P.D.
			Eigen::Matrix2d S = _mu * (F.transpose() * F - Eigen::Matrix2d::Identity());
			// S = clampMatrixEigenvalues(S);

			// 2nd Piola Stress tensor, P = F S, distribute to 2nd and 3rd node to get forces
			Matrix3x2d forces = -m_elements[i].A * F * S * m_elements[i].DmInverse.transpose();

			CVector3d f1(forces.col(0).data());
			CVector3d f2(forces.col(1).data());

			// Forces should sum to zero, so first vertex force is -(f_1 + f_2)
            // TODO: IMPL THIS:
            out_forces.col(0) -= (f1 + f2);
            out_forces.col(1) += f1;
            out_forces.col(2) += f2;
			/* *pCloth->m_PhysWorld.m_Particles[f->v1()->GetIndex()]->m_pForceVec -= f1 + f2; */
			/* *pCloth->m_PhysWorld.m_Particles[f->v2()->GetIndex()]->m_pForceVec += f1; */
			/* *pCloth->m_PhysWorld.m_Particles[f->v3()->GetIndex()]->m_pForceVec += f2; */

			// Compare forces to FD results
			const double fx = computeElementEnergy(i, x[0], x[1], x[2]);

			for (size_t j = 0; j < 3; j++)
			{
				CVector3d force = computeFDForce(j, fx, x,
					[&](Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3)
                                                 { return computeElementEnergy(i, x1, x2, x3); });
                fd_forces.col(j) = force;

				for (size_t k = j; k < 3; k++)
				{
                    fd_dFdx->block(f[j], f[k], 3, 3) = computeFDHessian(j, k, fx, x,
						[&](Eigen::Vector3d& x1, Eigen::Vector3d& x2, Eigen::Vector3d& x3) { return computeElementEnergy(i, x1, x2, x3); });
				}
			}

			// Take derivative of each force f_j with respect to each vertex k
			for (int j = 0; j < 3; j++)
			{
				// j-th force is A * F * S * Bm_j, where Bm_j is appropriate column(s) from Dm^-1
				Eigen::Vector2d Bm_j;
				if (j == 0)
					Bm_j = (-m_elements[i].DmInverse.row(0) - m_elements[i].DmInverse.row(1)).transpose();
				else
					Bm_j = m_elements[i].DmInverse.row(j-1).transpose();

				for (int k = j; k < 3; k++)
				{
					// Each column of 3x3 stiffness matrix is derivative of force f_j w.r.t. vertex k's l'th coordinate
					Eigen::Matrix3d K_jk;
					for (int l = 0; l < 3; l++)
					{
						// Get dFdu for k-th vertex, l-th coordinate
						Matrix3x2d& dFdu = m_elements[i].dFdu[k][l];
						K_jk.col(l) = (dFdu * S + _mu * F * (F.transpose() * dFdu + dFdu.transpose() * F)) * Bm_j * m_elements[i].A;
					}

					Matrix3x3 K;
                    K = K_jk;
					//K.SetCol(CVector3d(K_jk.col(0).data()), CVector3d(K_jk.col(1).data()), CVector3d(K_jk.col(2).data()));
					//dFdx->AddBlockMatrix(f->v(j)->GetIndex(), f->v(k)->GetIndex(), K);
                    dFdx->block(f[j], f[k], 3, 3) = K;
				}
			}
		}
	}
 #endif

	Eigen::Matrix2d clampMatrixEigenvalues(Eigen::Matrix2d& input)
	{
		static const double eps = 1.0e-6;

		Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> ev(input);	
		assert(ev.info() == Eigen::Success);

		Eigen::Vector2d evs;
		evs << std::max(ev.eigenvalues()[0], eps), std::max(ev.eigenvalues()[1], eps);
		return (ev.eigenvectors() * evs.asDiagonal() * ev.eigenvectors().transpose());
	}
};
