#include <StVK.h>
#include <obj_file.h>
#include <fstream>

void
compute(TriMesh::ptr mesh,
        const Vec3d &u0, const Vec3d &u1, const Vec3d &u2,
        std::ofstream &fout)
{
    // save
    auto v_old = mesh->v;

    // deform.
    mesh->v.row(0) += u0.transpose();
    mesh->v.row(1) += u1.transpose();
    mesh->v.row(2) += u2.transpose();

    // std::cout << "-----------------Configuration--------------------" << std::endl;
    // std::cout << mesh->v.transpose() << std::endl;
    // std::cout << "--------------------------------------------------" << std::endl;

    std::cout << "-----------------Deformation----------------------" << std::endl;
    std::cout << (mesh->v - v_old).transpose() << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    StVenantKirchoffForce force;

    force.precompute(mesh.get());

    double W = force.computeElementEnergy(0, mesh->v.row(0), mesh->v.row(1), mesh->v.row(2));

    std::cout << "Energy = " << W << std::endl;

    MatrixXd dFdx(9, 9);
    dFdx.setZero();

    Matrix3d forces = Matrix3d::Zero();
    Matrix3d fd_forces = Matrix3d::Zero();
    Matrix3d F_approx = Matrix3d::Zero();

    force.computeForces(mesh.get(), &dFdx, forces, fd_forces, F_approx);

    std::cout << "-----------------Forces---------------------------" << std::endl;
    std::cout << forces << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    std::cout << "-----------------F approx-------------------------" << std::endl;
    std::cout << F_approx << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    std::cout << "-----------------F finite-diff--------------------" << std::endl;
    std::cout << F_approx << std::endl;
    std::cout << "--------------------------------------------------" << std::endl;

    // // Coerce small effectively zero values to zero.
    // forces = forces.unaryExpr([](double x){return (fabs(x) < 1.e-10) ? 0.0 : x;});
    // dFdx = dFdx.unaryExpr([](double x){return (fabs(x) < 1.e-10) ? 0.0 : x;});
    // fd_forces = fd_forces.unaryExpr([](double x){return (fabs(x) < 1.e-10) ? 0.0 : x;});
    // F_approx = F_approx.unaryExpr([](double x){return (fabs(x) < 1.e-10) ? 0.0 : x;});

    // // Filter out nans and infs
    // // auto forces_rel_diff = ((forces - fd_forces).array() / forces.array()) \
    // //     .unaryExpr([](double x) { return (std::isinf(x) || std::isnan(x)) ? 0.0 : x; });
    // auto forces_abs_diff = (forces - fd_forces).array();
    // double forces_max_abs_err = forces_abs_diff.abs().maxCoeff();
    // std::cout << forces_max_abs_err << std::endl;

    // // auto Fapprox_rel_diff = ((forces - F_approx).array() / forces.array()) \
    // //     ; //.unaryExpr([](double x) { return (std::isinf(x) || std::isnan(x)) ? 0.0 : x; });
    // // double Fapprox_max_rel_err = Fapprox_rel_diff.abs().maxCoeff();
    // // std::cout << Fapprox_max_rel_err << std::endl;
    // auto Fapprox_abs_diff = (forces - F_approx).array();
    // double Fapprox_max_abs_err = Fapprox_abs_diff.abs().maxCoeff();
    // std::cout << Fapprox_max_abs_err << std::endl;

    // // output force magnitudes to file.
    // double disp = std::max(std::max(u0.norm(), u1.norm()), u2.norm());
    // fout << disp << ",";
    // fout << forces.col(0).norm() << ",";
    // fout << forces.col(1).norm() << ",";
    // fout << forces.col(2).norm() << ",";
    // fout << forces_max_abs_err << ",";
    // fout << Fapprox_max_abs_err << std::endl;

    // restore
    mesh->v = v_old;
}

void
deform(TriMesh::ptr mesh,
       const Vec3d &d0, const Vec3d &d1, const Vec3d &d2,
       const std::string &filepath)
{
    std::ofstream ofs(filepath, std::ofstream::out);

    const double step = 0.1;

    for (int i = 0; i < 2; ++i) {
        compute(mesh, i*step * d0, i*step * d1, i*step * d2, ofs);
    }

    ofs.close();
}

void
make_equilateral(TriMesh::ptr mesh)
{
    mesh->v.row(0) = Vec3d(0, 1, 0);
    mesh->vt.row(0) = Vec2d(0, 1);

    mesh->v.row(1) = Vec3d(-sqrt(3)/2.0, -0.5, 0);
    mesh->vt.row(1) = Vec2d(-sqrt(3)/2.0, -0.5);

    mesh->v.row(2) = Vec3d(sqrt(3)/2.0, -0.5, 0);
    mesh->vt.row(2) = Vec2d(sqrt(3)/2.0, -0.5);
}

int main(void)
{
    TriMesh::ptr mesh = ObjFile::read("tri.obj");
    // 0.
    // deform(mesh, Vec3d::Zero(), Vec3d::Zero(), Vec3d::Zero(), "0.csv");
    // deform(mesh, Vec3d(0, 1, 0), Vec3d(0, 1, 0), Vec3d(0, 1, 0), "0.csv");

    // // 1.
    // deform(mesh, Vec3d(0, 1, 0), Vec3d::Zero(), Vec3d::Zero(), "1.csv");

    // // 2.
    // deform(mesh, Vec3d(0, -1, 0), Vec3d::Zero(), Vec3d::Zero(), "2.csv");

    
    // deform(mesh, Vec3d::Zero(), Vec3d(1, 0, 0), Vec3d::Zero(), "3.csv");

    // // 4. 
    // deform(mesh, Vec3d::Zero(), Vec3d(-1, 0, 0), Vec3d::Zero(), "4.csv");

    make_equilateral(mesh);

    // 5.
    //deform(mesh, Vec3d(0, 1, 0).normalized(), Vec3d::Zero(), Vec3d::Zero(), "5.csv");

    // // 6.
    // deform(mesh, Vec3d(0, 1, 0).normalized(), Vec3d::Zero(), Vec3d::Zero(), "6.csv");

    // 7.
    deform(mesh,
          (-mesh->v.row(0)).normalized(),
          (-mesh->v.row(1)).normalized(),
          (-mesh->v.row(2)).normalized(), "7.csv");
    
    // // 8.
    // deform(mesh,
    //       (mesh->v.row(0)).normalized(),
    //       (mesh->v.row(1)).normalized(),
    //       (mesh->v.row(2)).normalized(), "8.csv");
    return 0;
}
