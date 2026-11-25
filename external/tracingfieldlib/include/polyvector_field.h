/*
 * //Non-Commercial License Agreement
 * //Copyright (c) 2025, Nico Pietroni
 * //
 * //1. Permission is hereby granted to use, copy, modify, and distribute
 * //   this software for academic research, teaching, and personal
 * //   non-commercial purposes, provided that this copyright notice is
 * //   included in all copies or substantial portions of the software.
 * //
 * //2. Commercial use of this software, including but not limited to use in
 * //   commercial products, services, tools, or consulting, is strictly
 * //   prohibited without prior written permission from the copyright holder.
 * //
 * //3. This software is provided "as is", without warranty of any kind,
 * //   express or implied, including but not limited to the warranties of
 * //   merchantability, fitness for a particular purpose, and
 * //   noninfringement. In no event shall the authors or copyright holders
 * //   be liable for any claim, damages, or other liability, arising from,
 * //   out of, or in connection with the software or the use or other
 * //   dealings in the software.
 * //
 * //4. For commercial licensing, please contact:
 * //   Nico Pietroni
 * //   Email: nico.pietroni@gmail.com
 */

//#ifndef POLYVECTOR_FIELD
//#define POLYVECTOR_FIELD
// Guard renamed: "POLYVECTOR_FIELD" clashes with directional::fieldTypeEnum::POLYVECTOR_FIELD
#ifndef POLYVECTOR_FIELD_H
#define POLYVECTOR_FIELD_H

//#include <Eigen/Geometry>
//#include <Eigen/SparseCholesky>
//#include <Eigen/Eigenvalues>
//#include <unsupported/Eigen/Polynomials>
//#include <igl/triangle_triangle_adjacency.h>
//#include <igl/doublearea.h>
//#include <igl/speye.h>
// #include <igl/PI.h>
// #include <iostream>


#include <Eigen/Sparse>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <igl/eigs.h>

class PolyVectorField
{
#define UNIFORM_WEIGHTS 0
#define BARYCENTRIC_WEIGHTS 1
#define INV_COT_WEIGHTS 2

    inline static void multiply_polynomials(Eigen::RowVectorXcd &p1,
                                            const Eigen::RowVectorXcd &p2)
    {
        Eigen::RowVectorXcd newp = Eigen::RowVectorXcd::Zero(p1.size() + p2.size());
        for (int i = 0; i < p1.size(); i++)
            for (int j = 0; j < p2.size(); j++)
                newp(i + j) += p1(i) * p2(j);

        p1 = newp;
    }

    inline static void multiply_polynomials(Eigen::MatrixXcd &p1,
                                            const Eigen::MatrixXcd &p2)
    {
        assert(p1.rows() == p2.rows());
        Eigen::MatrixXcd newp = Eigen::MatrixXcd::Zero(p1.rows(), p1.cols() + p2.cols());
        for (int i = 0; i < p1.cols(); i++)
            for (int j = 0; j < p2.cols(); j++)
                newp.col(i + j).array() += p1.col(i).array() * p2.col(j).array();

        p1 = newp;
    }

    inline static void polynomial_eval_from_roots(const Eigen::RowVectorXcd &roots,
                                                  const std::complex<double> &evalPoint,
                                                  std::complex<double> &polyValue)
    {
        polyValue = std::complex<double>(1.0, 0.0);
        for (int i = 0; i < roots.size(); i++)
            polyValue *= (evalPoint - roots(i));
    }

    inline static void polynomial_eval_from_roots(const Eigen::MatrixXcd &roots,
                                                  const Eigen::VectorXcd &evalPoints,
                                                  Eigen::VectorXcd &polyValues)
    {
        assert(evalPoints.rows() == roots.rows());
        polyValues = Eigen::VectorXcd::Constant(evalPoints.rows(), std::complex<double>(1.0, 0.0));
        for (int i = 0; i < roots.cols(); i++)
            polyValues.array() *= (evalPoints - roots.col(i)).array(); // rowwise-product faster?;
    }

    inline static void polynomial_eval(const Eigen::RowVectorXcd &coeffs,
                                       const std::complex<double> &evalPoint,
                                       std::complex<double> &polyValue)
    {
        polyValue = std::complex<double>(0.0, 0.0);
        std::complex<double> z = 1;
        for (int i = 0; i < coeffs.size(); i++)
        {
            polyValue += coeffs(i) * z;
            z *= evalPoint;
        }
        std::cout << "z: " << z << std::endl;
        polyValue += z; // the biggest and unit power
    }

    inline static void polynomial_eval(const Eigen::MatrixXcd &coeffs,
                                       const Eigen::VectorXcd &evalPoints,
                                       Eigen::VectorXcd &polyValues)
    {
        polyValues = Eigen::VectorXcd::Zero(coeffs.rows());
        Eigen::VectorXcd z = Eigen::VectorXcd::Constant(coeffs.rows(), 1.0);
        for (int i = 0; i < coeffs.cols(); i++)
        {
            polyValues.array() += coeffs.col(i).array() * z.array();
            z.array() *= evalPoints.array();
        }
        // std::cout<<"z: "<<z<<std::endl;
        polyValues.array() += z.array(); // the biggest and unit power
    }

    inline static void complex_eigs(const Eigen::SparseMatrix<std::complex<double>> &Q,
                                    const Eigen::SparseMatrix<std::complex<double>> M,
                                    const int numEigs,
                                    Eigen::MatrixXcd &U,
                                    Eigen::VectorXcd &S)
    {
        // TODO: should insist on real matrix M also nominally
        Eigen::SparseMatrix<double> QReal(2 * Q.rows(), 2 * Q.rows());
        std::vector<Eigen::Triplet<double>> QRealTriplets, MRealTriplets;
        for (int k = 0; k < Q.outerSize(); ++k)
            for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(Q, k); it; ++it)
            {
                QRealTriplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value().real()));
                QRealTriplets.push_back(Eigen::Triplet<double>(it.row(), Q.cols() + it.col(), -it.value().imag()));
                QRealTriplets.push_back(Eigen::Triplet<double>(Q.rows() + it.row(), it.col(), it.value().imag()));
                QRealTriplets.push_back(Eigen::Triplet<double>(Q.rows() + it.row(), Q.cols() + it.col(), it.value().real()));
            }
        QReal.setFromTriplets(QRealTriplets.begin(), QRealTriplets.end());

        Eigen::SparseMatrix<double> MReal(2 * M.rows(), 2 * M.rows());
        for (int k = 0; k < M.outerSize(); ++k)
            for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(M, k); it; ++it)
            {
                MRealTriplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value().real()));
                MRealTriplets.push_back(Eigen::Triplet<double>(it.row() + M.rows(), it.col() + M.cols(), it.value().real()));
            }

        MReal.setFromTriplets(MRealTriplets.begin(), MRealTriplets.end());
        Eigen::MatrixXd UReal;
        Eigen::VectorXd SReal;

        igl::eigs(QReal, MReal, numEigs, igl::EIGS_TYPE_SM, UReal, SReal);

        U.resize(UReal.rows() / 2, UReal.cols());
        S = SReal;

        U = UReal.block(0, 0, UReal.rows() / 2, UReal.cols()).cast<std::complex<double>>().array() * std::complex<double>(1, 0) +
            UReal.block(UReal.rows() / 2, 0, UReal.rows() / 2, UReal.cols()).cast<std::complex<double>>().array() * std::complex<double>(0, 1);
    }

    inline static void circumcircle(const Eigen::MatrixXd &V,
                                    const Eigen::MatrixXi &F,
                                    Eigen::MatrixXd &circumcenter,
                                    Eigen::VectorXd &circumradius)
    {
        // Using the formulation from https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_from_cross-_and_dot-products
        using namespace Eigen;

        circumcenter.resize(F.rows(), 3);
        circumradius.resize(F.rows());
        for (int i = 0; i < F.rows(); i++)
        {
            RowVector3d P12 = V.row(F(i, 1)) - V.row(F(i, 0));
            RowVector3d P23 = V.row(F(i, 2)) - V.row(F(i, 1));
            RowVector3d P31 = V.row(F(i, 0)) - V.row(F(i, 2));

            double denominator = 2.0 * P12.cross(P23).squaredNorm();
            double alpha = -P12.dot(P31) * P23.squaredNorm();
            double beta = -P12.dot(P23) * P31.squaredNorm();
            double gamma = -P31.dot(P23) * P12.squaredNorm();
            circumcenter.row(i) = V.row(F(i, 0)).array() * alpha + V.row(F(i, 1)).array() * beta + V.row(F(i, 2)).array() * gamma;

            // testing to see if face center is correct
            double radius1 = (circumcenter.row(i) - V.row(F(i, 0))).squaredNorm();
            double radius2 = (circumcenter.row(i) - V.row(F(i, 1))).squaredNorm();
            double radius3 = (circumcenter.row(i) - V.row(F(i, 2))).squaredNorm();

            std::cout << "radii1,2,3: " << radius1 << " " << radius2 << " " << radius3 << std::endl;

            circumradius(i) = (circumcenter.row(i) - V.row(F(i, 0))).norm();
        }
    }

    struct PolyVectorData
    {
    public:
        // User parameters
        Eigen::VectorXi constFaces;   // list of faces where there are (partial) constraints. The faces can repeat to constrain more vectors
        Eigen::MatrixXd constVectors; // corresponding to constFaces.

        int N;                      // Degree of field
        int sizeF;                  // #faces
        bool signSymmetry;          // Whether field enforces a ssign symmetry (only when N is even, otherwise by default set to false)
        bool perfectRoSy;           // Whether the field must be perfect rotationally-symmetric (but not unit).
        double wSmooth;             // Weight of smoothness
        double wRoSy;               // Weight of rotational-symmetry. "-1" means a perfect RoSy field (power field)
        Eigen::VectorXd wAlignment; // Weight of alignment per each of the constfaces. "-1" means a fixed vector

        int lapType; // Choice of weights (from UNIFORM_WEIGHTS,BARYCENTRIC_WEIGHTS or INV_COT_WEIGHTS)

        Eigen::SparseMatrix<std::complex<double>> smoothMat; // Smoothness energy
        Eigen::SparseMatrix<std::complex<double>> roSyMat;   // Rotational-symmetry energy
        Eigen::SparseMatrix<std::complex<double>> alignMat;  //(soft) alignment energy.
        Eigen::SparseMatrix<std::complex<double>> reducMat;  // reducing the fixed dofs (for instance with sign symmetry or fixed partial constraints)
        Eigen::VectorXcd reducRhs;                           // The uncompressed PV coeffs are reducMat*true_dofs+reducRhs
        Eigen::VectorXcd alignRhs;                           // encoding the soft constraints

        // Mass and stiffness matrices
        Eigen::SparseMatrix<std::complex<double>> WSmooth, WAlign, WRoSy, M;
        double totalRoSyWeight, totalConstrainedWeight, totalSmoothWeight; // for co-scaling energies

        PolyVectorData() : signSymmetry(true), lapType(BARYCENTRIC_WEIGHTS), wSmooth(1.0), wRoSy(0.0)
        {
            wAlignment.resize(0);
            constFaces.resize(0);
            constVectors.resize(0, 3);
        }
        ~PolyVectorData() {}
    };

    // Precalculate the operators according to the user-prescribed parameters. Must be called whenever any of them changes
    // Inputs:
    //  V:      #V by 3 vertex coordinates.
    //  F:      #F by 3 face vertex indices.
    //  EV:     #E by 2 matrix of edges (vertex indices)
    //  EF:     #E by 2 matrix of oriented adjacent faces
    //  B1, B2: #F by 3 matrices representing the local base of each face.
    //  PolyVectorData (must fill non-default values in advance)
    // Outputs:
    //  PolyVectorData:       Updated structure with all operators
    inline static void Precompute(const Eigen::MatrixXd &V,
                                  const Eigen::MatrixXi &F,
                                  const Eigen::MatrixXi &EV,
                                  const Eigen::MatrixXi &EF,
                                  const Eigen::MatrixXd &B1,
                                  const Eigen::MatrixXd &B2,
                                  const int N,
                                  PolyVectorData &pvData)
    {

        using namespace std;
        using namespace Eigen;

        // Building the smoothness matrices, with an energy term for each inner edge and degree
        int rowCounter = 0;
        std::vector<Triplet<complex<double>>> dTriplets, WTriplets;
        pvData.N = N;
        pvData.sizeF = F.rows();
        if (pvData.N % 2 != 0)
            pvData.signSymmetry = false; // it has to be for odd N

        /************Smoothness matrices****************/
        VectorXd stiffnessWeights = VectorXd::Zero(EF.rows());

        if (pvData.lapType == UNIFORM_WEIGHTS)
        {
            stiffnessWeights.setOnes();
        }
        else
        {
            MatrixXd faceCenter;
            VectorXd radius; // stub
            if (pvData.lapType == BARYCENTRIC_WEIGHTS)
                igl::barycenter(V, F, faceCenter);
            if (pvData.lapType == INV_COT_WEIGHTS)
                circumcircle(V, F, faceCenter, radius);
            for (int i = 0; i < EF.rows(); i++)
            {
                if ((EF(i, 0) == -1) || (EF(i, 1) == -1))
                    continue; // boundary edge

                RowVector3d midEdge = (V.row(EV(i, 0)) + V.row(EV(i, 1))) / 2.0;
                double dualLength = (faceCenter.row(EF(i, 0)) - midEdge).norm() + (faceCenter.row(EF(i, 1)) - midEdge).norm(); // TODO: that's only an approximation of barycentric height...
                double primalLength = (V.row(EV(i, 0)) - V.row(EV(i, 1))).norm();
                if (dualLength > 10e-9) // smaller than 10e-9 might happen for inv cot weights
                    stiffnessWeights(i) = primalLength / dualLength;
            }
        }

        pvData.totalSmoothWeight = stiffnessWeights.sum();

        vector<Triplet<complex<double>>> WSmoothTriplets, MTriplets;
        VectorXd doubleAreas;
        igl::doublearea(V, F, doubleAreas);
        for (int n = 0; n < pvData.N; n++)
        {
            for (int i = 0; i < EF.rows(); i++)
            {
                if ((EF(i, 0) == -1) || (EF(i, 1) == -1))
                    continue; // boundary edge

                // Compute the complex representation of the common edge
                RowVector3d e = V.row(EV(i, 1)) - V.row(EV(i, 0));
                RowVector2d vef = Vector2d(e.dot(B1.row(EF(i, 0))), e.dot(B2.row(EF(i, 0)))).normalized();
                complex<double> ef(vef(0), vef(1));
                Vector2d veg = Vector2d(e.dot(B1.row(EF(i, 1))), e.dot(B2.row(EF(i, 1)))).normalized();
                complex<double> eg(veg(0), veg(1));

                // Add the term conj(f)^n*ui - conj(g)^n*uj to the differential matrix
                dTriplets.push_back(Triplet<complex<double>>(rowCounter, n * F.rows() + EF(i, 0), pow(conj(ef), pvData.N - n)));
                dTriplets.push_back(Triplet<complex<double>>(rowCounter, n * F.rows() + EF(i, 1), -1. * pow(conj(eg), pvData.N - n)));

                // stiffness weights
                WSmoothTriplets.push_back(Triplet<complex<double>>(rowCounter, rowCounter, stiffnessWeights(i)));
                rowCounter++;
            }

            for (int i = 0; i < F.rows(); i++)
                MTriplets.push_back(Triplet<complex<double>>(n * F.rows() + i, n * F.rows() + i, doubleAreas(i) / 2.0));
        }

        pvData.smoothMat.resize(rowCounter, pvData.N * F.rows());
        pvData.smoothMat.setFromTriplets(dTriplets.begin(), dTriplets.end());

        pvData.WSmooth.resize(rowCounter, rowCounter);
        pvData.WSmooth.setFromTriplets(WSmoothTriplets.begin(), WSmoothTriplets.end());

        pvData.M.resize(pvData.N * F.rows(), pvData.N * F.rows());
        pvData.M.setFromTriplets(MTriplets.begin(), MTriplets.end());

        // creating reduction transformation
        VectorXi numFaceConstraints = Eigen::VectorXi::Zero(F.rows());
        int realN = (pvData.signSymmetry ? N / 2 : N);
        realN = (pvData.wRoSy < 0.0 ? 1 : realN);
        // MatrixXcd faceConstraints(F.rows(),realN);
        std::vector<MatrixXcd> localFaceReducMats;
        localFaceReducMats.resize(F.rows());
        std::vector<VectorXcd> localFaceReducRhs;
        localFaceReducRhs.resize(F.rows());

        for (int i = 0; i < F.rows(); i++)
        {
            localFaceReducMats[i] = MatrixXcd::Identity(realN, realN);
            localFaceReducRhs[i] = VectorXcd::Zero(realN);
        }

        /*************Hard-constraint reduction matrices******************/
        for (int i = 0; i < pvData.constFaces.size(); i++)
        {
            if (pvData.wAlignment(i) >= 0.0)
                continue; // here we only handle the reduction caused by a fixed dof
            if (numFaceConstraints(pvData.constFaces(i)) == realN)
                continue; // overconstrained; we ignore any further constraints on that face

            int currFaceNumDof = numFaceConstraints(pvData.constFaces(i));

            complex<double> constVectorComplexSingle = std::complex<double>(pvData.constVectors.row(i).dot(B1.row(pvData.constFaces(i))), pvData.constVectors.row(i).dot(B2.row(pvData.constFaces(i))));
            complex<double> constVectorComplex = (pvData.signSymmetry ? constVectorComplexSingle * constVectorComplexSingle : constVectorComplexSingle);
            constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorComplexSingle, N) : constVectorComplex);

            MatrixXcd singleReducMat = MatrixXcd::Zero(realN - currFaceNumDof, realN - currFaceNumDof - 1);
            VectorXcd singleReducRhs = VectorXcd::Zero(realN - currFaceNumDof);
            for (int j = 0; j < realN - currFaceNumDof - 1; j++)
            {
                singleReducMat(j, j) = -constVectorComplex;
                singleReducMat(j + 1, j) = complex<double>(1.0, 0.0);
            }
            singleReducRhs(realN - currFaceNumDof - 1) = -constVectorComplex;

            localFaceReducRhs[pvData.constFaces(i)] = localFaceReducMats[pvData.constFaces(i)] * singleReducRhs + localFaceReducRhs[pvData.constFaces(i)];
            localFaceReducMats[pvData.constFaces(i)] = localFaceReducMats[pvData.constFaces(i)] * singleReducMat;

            numFaceConstraints(pvData.constFaces(i))++;

            // faceConstraints(pvData.constFaces(i), numFaceConstraints(pvData.constFaces(i))++) = constVectorComplex;
        }

        // creating the global reduction matrices
        double colCounter = 0;
        pvData.reducRhs = VectorXcd::Zero(pvData.N * F.rows());
        vector<Triplet<complex<double>>> reducMatTriplets;
        int jump = (pvData.signSymmetry ? 2 : 1);
        jump = (pvData.wRoSy < 0.0 ? pvData.N : jump);
        for (int i = 0; i < F.rows(); i++)
        {
            // std::cout<<"localFaceReducMats[i]: "<<localFaceReducMats[i]<<std::endl;
            for (int j = 0; j < pvData.N; j += jump)
            {
                for (int k = 0; k < localFaceReducMats[i].cols(); k++)
                    reducMatTriplets.push_back(Triplet<complex<double>>(j * F.rows() + i, colCounter + k, localFaceReducMats[i](j / jump, k)));

                pvData.reducRhs(j * F.rows() + i) = localFaceReducRhs[i](j / jump);
            }

            colCounter += localFaceReducMats[i].cols();
        }

        pvData.reducMat.resize(pvData.N * F.rows(), colCounter);
        pvData.reducMat.setFromTriplets(reducMatTriplets.begin(), reducMatTriplets.end());

        /****************rotational-symmetry matrices********************/

        if (pvData.wRoSy >= 0.0)
        { // this is anyhow enforced, this matrix is unnecessary)
            vector<Triplet<complex<double>>> roSyTriplets, WRoSyTriplets;
            for (int i = F.rows(); i < pvData.N * F.rows(); i++)
            {
                roSyTriplets.push_back(Triplet<complex<double>>(i, i, 1.0));
                WRoSyTriplets.push_back(Triplet<complex<double>>(i, i, doubleAreas(i % F.rows()) / 2.0));
            }

            pvData.roSyMat.resize(N * F.rows(), N * F.rows());
            pvData.roSyMat.setFromTriplets(roSyTriplets.begin(), roSyTriplets.end());

            pvData.WRoSy.resize(N * F.rows(), N * F.rows());
            pvData.WRoSy.setFromTriplets(WRoSyTriplets.begin(), WRoSyTriplets.end());

            pvData.totalRoSyWeight = ((double)pvData.N) * doubleAreas.sum() / 2.0;
        }
        else
        {
            pvData.roSyMat.resize(0, N * F.rows());
            pvData.WRoSy.resize(0, 0); // Eeven necessary?
            pvData.totalRoSyWeight = 1.0;
        }

        /*****************Soft alignment matrices*******************/
        rowCounter = 0;
        vector<Triplet<complex<double>>> alignTriplets;
        vector<VectorXcd> alignRhsList;
        vector<Triplet<complex<double>>> WAlignTriplets;
        pvData.totalConstrainedWeight = 0.0;
        bool noSoftAlignment = true;
        for (int i = 0; i < pvData.constFaces.size(); i++)
        {
            if (pvData.wAlignment(i) < 0.0)
                continue; // here we only handle soft alignments

            noSoftAlignment = false;

            complex<double> constVectorComplexSingle = std::complex<double>(pvData.constVectors.row(i).dot(B1.row(pvData.constFaces(i))), pvData.constVectors.row(i).dot(B2.row(pvData.constFaces(i))));
            complex<double> constVectorComplex = (pvData.signSymmetry ? constVectorComplexSingle * constVectorComplexSingle : constVectorComplexSingle);
            constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorComplexSingle, N) : constVectorComplex);
            numFaceConstraints(pvData.constFaces(i))++;
            // faceConstraints(pvData.constFaces(i), numFaceConstraints(pvData.constFaces(i))++) = constVectorComplex;

            MatrixXcd singleReducMat = MatrixXcd::Zero(realN, realN - 1);
            VectorXcd singleReducRhs = VectorXcd::Zero(realN);
            for (int j = 0; j < realN - 1; j++)
            {
                singleReducMat(j, j) = -constVectorComplex;
                singleReducMat(j + 1, j) = complex<double>(1.0, 0.0);
            }
            singleReducRhs(realN - 1) = -constVectorComplex;

            MatrixXcd IAiA;
            if (realN > 1)
            {
                MatrixXcd invSingleReducMat = singleReducMat.completeOrthogonalDecomposition().pseudoInverse();
                IAiA = MatrixXcd::Identity(realN, realN) - singleReducMat * invSingleReducMat;
            }
            else
                IAiA = MatrixXcd::Ones(realN, realN);
            singleReducRhs = IAiA * singleReducRhs;
            for (int j = 0; j < IAiA.rows(); j++)
                for (int k = 0; k < IAiA.cols(); k++)
                    alignTriplets.push_back(Triplet<complex<double>>(rowCounter + j, k * jump * F.rows() + pvData.constFaces(i), IAiA(j, k)));

            alignRhsList.push_back(singleReducRhs);
            for (int j = 0; j < singleReducRhs.size(); j++)
            {
                WAlignTriplets.push_back(Triplet<complex<double>>(rowCounter + j, rowCounter + j, pvData.wAlignment(i) * doubleAreas(pvData.constFaces(i)) / 2.0));
                pvData.totalConstrainedWeight += doubleAreas(pvData.constFaces(i)) / 2.0;
            }
            rowCounter += realN;
        }

        if (noSoftAlignment)
            pvData.totalConstrainedWeight = 1.0; // it wouldn't be used, except just to avoid a division by zero in the energy formulation

        pvData.alignRhs.resize(rowCounter);
        for (int i = 0; i < alignRhsList.size(); i++)
            pvData.alignRhs.segment(i * realN, realN) = alignRhsList[i];

        pvData.alignMat.resize(rowCounter, N * F.rows());
        pvData.alignMat.setFromTriplets(alignTriplets.begin(), alignTriplets.end());

        pvData.WAlign.resize(rowCounter, rowCounter);
        pvData.WAlign.setFromTriplets(WAlignTriplets.begin(), WAlignTriplets.end());
    }

public:
    // Computes a polyvector on the entire mesh
    // Inputs:
    //  PolyVectorData: The data structure which should have been initialized with polyvector_precompute()
    // Outputs:
    //  polyVectorField: #F by N The output interpolated field, in polyvector (complex polynomial) format.
    static inline void Construct(const PolyVectorData &pvData,
                                 Eigen::MatrixXcd &polyVectorField)
    {
        using namespace std;
        using namespace Eigen;

        // forming total energy matrix;
        SparseMatrix<complex<double>> totalUnreducedLhs = pvData.wSmooth * (pvData.smoothMat.adjoint() * pvData.WSmooth * pvData.smoothMat) / pvData.totalSmoothWeight + (pvData.wRoSy * pvData.roSyMat.adjoint() * pvData.WRoSy * pvData.roSyMat) / pvData.totalRoSyWeight + (pvData.alignMat.adjoint() * pvData.WAlign * pvData.alignMat) / pvData.totalConstrainedWeight;
        VectorXcd totalUnreducedRhs = (pvData.alignMat.adjoint() * pvData.WAlign * pvData.alignRhs) / pvData.totalConstrainedWeight;

        SparseMatrix<complex<double>> totalLhs = pvData.reducMat.adjoint() * totalUnreducedLhs * pvData.reducMat;
        VectorXcd totalRhs = pvData.reducMat.adjoint() * (totalUnreducedRhs - totalUnreducedLhs * pvData.reducRhs);
        polyVectorField = MatrixXcd::Zero(pvData.sizeF, pvData.N);
        if (pvData.constFaces.size() == 0) // alignmat should be empty and the reduction matrix should be only sign symmetry, if applicable
        {
            // using a matrix with only the first sizeF x sizeF block
            vector<Triplet<complex<double>>> X0LhsTriplets, X0MTriplets;
            SparseMatrix<complex<double>> X0Lhs, X0M;
            for (int k = 0; k < totalUnreducedLhs.outerSize(); ++k)
                for (SparseMatrix<std::complex<double>>::InnerIterator it(totalUnreducedLhs, k); it; ++it)
                    if ((it.row() < pvData.sizeF) && (it.col() < pvData.sizeF))
                        X0LhsTriplets.push_back(Triplet<complex<double>>(it.row(), it.col(), it.value() * pvData.totalSmoothWeight)); // to bypass the early convergence of igl eigenvalues...

            X0Lhs.resize(pvData.sizeF, pvData.sizeF);
            X0Lhs.setFromTriplets(X0LhsTriplets.begin(), X0LhsTriplets.end());

            for (int k = 0; k < pvData.M.outerSize(); ++k)
                for (SparseMatrix<std::complex<double>>::InnerIterator it(pvData.M, k); it; ++it)
                    if ((it.row() < pvData.sizeF) && (it.col() < pvData.sizeF))
                        X0MTriplets.push_back(Triplet<complex<double>>(it.row(), it.col(), it.value() * pvData.totalSmoothWeight));

            X0M.resize(pvData.sizeF, pvData.sizeF);
            X0M.setFromTriplets(X0MTriplets.begin(), X0MTriplets.end());

            // Extracting first eigenvector
            Eigen::MatrixXcd U;
            Eigen::VectorXcd S;
            complex_eigs(X0Lhs, X0M, 10, U, S);
            int smallestIndex;
            S.cwiseAbs().minCoeff(&smallestIndex);

            polyVectorField.col(0) = U.col(smallestIndex);
        }
        else
        { // just solving the system
            SimplicialLDLT<SparseMatrix<complex<double>>> solver;
            // solver.analyzePattern(totalLhs);   // for this step the numerical values of A are not used
            solver.compute(totalLhs);
            VectorXcd reducedDofs = solver.solve(totalRhs);
            assert(solver.info() == Success);
            VectorXcd fullDofs = pvData.reducMat * reducedDofs + pvData.reducRhs;
            for (int i = 0; i < pvData.N; i++)
                polyVectorField.col(i) = fullDofs.segment(i * pvData.sizeF, pvData.sizeF);

            // std::cout<<"Smoothness energy: "<<pvData.wSmooth * (fullDofs.adjoint()*pvData.smoothMat.adjoint()*pvData.WSmooth*pvData.smoothMat*fullDofs)/pvData.totalSmoothWeight<<std::endl;
            // std::cout<<"RoSy Energy: "<<fullDofs.adjoint()* ((pvData.wRoSy*pvData.roSyMat.adjoint()*pvData.WRoSy*pvData.roSyMat)/pvData.totalRoSyWeight)*fullDofs<<std::endl;
            // TODO: measure energy in the correct way
            /*std::cout<<"Alignment energy: "<< fullDofs.adjoint()*((pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignMat*fullDofs - pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignRhs)/pvData.totalConstrainedWeight)<<std::endl;*/
        }
    }
    // minimal version without auxiliary data
    inline static void Construct(const Eigen::MatrixXd &V,
                                 const Eigen::MatrixXi &F,
                                 const Eigen::VectorXi &constFaces,
                                 const Eigen::MatrixXd &constVectors,
                                 const double smoothWeight,
                                 const double roSyWeight,
                                 const Eigen::VectorXd &alignWeights,
                                 const int N,
                                 Eigen::MatrixXcd &polyVectorField)
    {
        Eigen::MatrixXi EV, xi, EF;
        Eigen::MatrixXd B1, B2, xd;
        igl::local_basis(V, F, B1, B2, xd);
        PolyVectorData pvData;
        pvData.constFaces = constFaces;
        pvData.constVectors = constVectors;
        pvData.wAlignment = alignWeights;
        pvData.wSmooth = smoothWeight;
        pvData.wRoSy = roSyWeight;
        igl::edge_topology(V, F, EV, xi, EF);
        Precompute(V, F, EV, EF, B1, B2, N, pvData);
        Construct(pvData, polyVectorField);
    }

    // Converts a field in PolyVector representation to raw represenation.
    // Inputs:
    //  B1, B2:           #F by 3 matrices representing the local base of each face.
    //  polyVectorField:  #F by N complex PolyVectors
    //  N:                The degree of the field.
    // Outputs:
    //  raw:              #F by 3*N matrix with all N explicit vectors of each directional in raw format xyzxyz
    //    returns true if succeeded
    inline static bool PolyvectorToField(const Eigen::MatrixXd &B1,
                                         const Eigen::MatrixXd &B2,
                                         const Eigen::MatrixXcd &pvField,
                                         const int N,
                                         Eigen::MatrixXd &rawField,
                                         bool signSymmetry = true,
                                         const double rootTolerance = 1e-8)
    {
        using namespace std;
        using namespace Eigen;
        rawField.resize(B1.rows(), 3 * N);
        if (N % 2 != 0)
            signSymmetry = false; // by definition
        MatrixXcd actualPVField;
        int actualN;
        if (signSymmetry)
        {
            actualPVField.resize(pvField.rows(), N / 2);
            for (int i = 0; i < N; i += 2)
                actualPVField.col(i / 2) = pvField.col(i);
            actualN = N / 2;
        }
        else
        {
            actualPVField = pvField;
            actualN = N;
        }

        MatrixXcd roots(actualPVField.rows(), actualN);
        roots.col(0).array() = (-actualPVField.col(0).array()).pow(1.0 / (double)actualN);
        for (int i = 1; i < actualN; i++)
            roots.col(i).array() = roots.col(i - 1).array() * std::exp(std::complex<double>(0, 2.0 * M_PI / (double)actualN));

        MatrixXd rootError = MatrixXd::Constant(actualPVField.rows(), actualN, 1000.0);
        double maxError = 1000.0;
        int currRoot = 0;
        MatrixXcd mostRoots(actualPVField.rows(), actualN - 1);
        int maxIterations = 1000;
        int currIteration = 0;
        do
        {
            mostRoots.block(0, 0, mostRoots.rows(), currRoot) = roots.block(0, 0, roots.rows(), currRoot);
            mostRoots.block(0, currRoot, mostRoots.rows(), actualN - currRoot - 1) = roots.block(0, currRoot + 1, roots.rows(), actualN - currRoot - 1);
            VectorXcd numerator;
            polynomial_eval(actualPVField, roots.col(currRoot), numerator);
            VectorXcd denominator;
            polynomial_eval_from_roots(mostRoots, roots.col(currRoot), denominator);
            roots.col(currRoot).array() -= numerator.array() / denominator.array();
            rootError.col(currRoot) = numerator.cwiseAbs();
            maxError = rootError.lpNorm<Infinity>(); // optimize by only the relevant column?
            currRoot = (currRoot + 1) % actualN;
            currIteration++;
        } while ((maxError > rootTolerance) && (currIteration < maxIterations));

        if (currIteration >= maxIterations)
            return false;

        if (signSymmetry)
            roots = roots.cwiseSqrt();

        for (int f = 0; f < roots.rows(); f++)
        {
            RowVectorXcd rowRoots = roots.row(f);
            std::sort(rowRoots.data(), rowRoots.data() + rowRoots.size(), [](std::complex<double> a, std::complex<double> b)
                      { return arg(a) < arg(b); });
            roots.row(f) = rowRoots;
        }

        if (signSymmetry)
        {
            MatrixXcd actualRoots(roots.rows(), 2 * roots.cols());
            actualRoots << roots, -roots;
            roots = actualRoots;
        }

        for (int i = 0; i < N; i++)
            for (int f = 0; f < pvField.rows(); f++)
            {
                std::complex<double> root = roots(f, i);
                rawField.block<1, 3>(f, 3 * i) = B1.row(f) * root.real() + B2.row(f) * root.imag();
            }

        return true;
    }
};

#endif
