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

#ifndef CATMULL_CLARKABILITY
#define CATMULL_CLARKABILITY

#include "../solve/least_squares_lin_constraint_solver.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1, ScalarType L2, ScalarType L3) {
  ScalarType CC = L1 + L2 - L3;

  // CC+=std::min(L1,std::min(L2,L3));
  // if (CC<0)return 0;

  return (CC);
}

template <class ScalarType>
ScalarType OrderedCC(ScalarType L1, ScalarType L2, ScalarType L3, ScalarType L4,
                     ScalarType L5) {
  ScalarType CC = (L3 + L4 + L5) - (L1 + L2);

  // CC+=std::min(L1,std::min(L2,std::min(L3,std::min(L4,L5))));
  // if (CC<0)return 0;

  return CC;
}

template <class ScalarType>
ScalarType CatmullClarkability3(const ScalarType &L1, const ScalarType &L2,
                                const ScalarType &L3,
                                const ScalarType &SideThr) {
  ScalarType minL = std::min(L1, std::min(L2, L3));
  ScalarType CC0 = OrderedCC(L1, L2, L3);
  ScalarType CC1 = OrderedCC(L2, L3, L1);
  ScalarType CC2 = OrderedCC(L3, L1, L2);
  ScalarType MinCC = std::min(CC0, std::min(CC1, CC2));
  MinCC += minL * SideThr;
  // return (std::min(CC0,std::min(CC1,CC2)));
  return (MinCC);
}

template <class ScalarType>
ScalarType CatmullClarkability4(const ScalarType &L1, const ScalarType &L2,
                                const ScalarType &L3, const ScalarType &L4,
                                const ScalarType &SideThr) {
  ScalarType minL = std::min(L1, std::min(L2, std::min(L3, L4)));
  ScalarType CC0 = std::min(L1, L3) - std::max(L1, L3);
  ScalarType CC1 = std::min(L2, L4) - std::max(L2, L4);
  ScalarType CC = std::min(CC0, CC1);
  CC += minL * SideThr;
  // CC+=std::min(L1,std::min(L2,std::min(L3,L4)));
  return CC;
}

template <class ScalarType>
ScalarType CatmullClarkability5(const ScalarType &L1, const ScalarType &L2,
                                const ScalarType &L3, const ScalarType &L4,
                                const ScalarType &L5,
                                const ScalarType &SideThr) {
  ScalarType minL = std::min(L1, std::min(L2, std::min(L3, std::min(L4, L5))));
  ScalarType CC0 = OrderedCC(L1, L2, L3, L4, L5);
  ScalarType CC1 = OrderedCC(L2, L3, L4, L5, L1);
  ScalarType CC2 = OrderedCC(L3, L4, L5, L1, L2);
  ScalarType CC3 = OrderedCC(L4, L5, L1, L2, L3);
  ScalarType CC4 = OrderedCC(L5, L1, L2, L3, L4);
  ScalarType MinCC =
      std::min(CC0, std::min(CC1, std::min(CC2, std::min(CC3, CC4))));
  MinCC += minL * SideThr;
  return (MinCC);
}

template <class ScalarType>
ScalarType CatmullClarkability6(const ScalarType &L1, const ScalarType &L2,
                                const ScalarType &L3, const ScalarType &L4,
                                const ScalarType &L5, const ScalarType &L6,
                                const ScalarType &SideThr) {
  ScalarType minL =
      std::min(L1, std::min(L2, std::min(L3, std::min(L4, std::min(L5, L6)))));
  // OLD (one triangle inequality per triple, and SideThr ignored in the
  // returned value)
  // ScalarType CC0 = OrderedCC(L1, L3, L5);
  // ScalarType CC1 = OrderedCC(L2, L4, L6);
  // ScalarType MinCC = std::min(CC0, CC1);
  // MinCC += minL * SideThr;
  // return (std::min(CC0, CC1));
  // NEW (all three inequalities of each triple, SideThr applied): the even
  // sides (L1,L3,L5) and the odd ones (L2,L4,L6) must each be a valid
  // triangle, like CatmullClarkability3
  ScalarType CC0 = std::min(OrderedCC(L1, L3, L5),
                            std::min(OrderedCC(L3, L5, L1), OrderedCC(L5, L1, L3)));
  ScalarType CC1 = std::min(OrderedCC(L2, L4, L6),
                            std::min(OrderedCC(L4, L6, L2), OrderedCC(L6, L2, L4)));
  ScalarType MinCC = std::min(CC0, CC1);
  MinCC += minL * SideThr;
  return (MinCC);
}

template <class ScalarType>
bool IsCatmullClarkable(const std::vector<ScalarType> &EdgeL,
                        const ScalarType &SideThr, bool SkipValence4,
                        bool print_debug = false) {
  if (EdgeL.size() < 3)
    return false;
  if (EdgeL.size() > 6)
    return false;

  if ((SkipValence4) && (EdgeL.size() == 4))
    return true;

  ScalarType CC = 0;
  if (EdgeL.size() == 3)
    CC = (CatmullClarkability3(EdgeL[0], EdgeL[1], EdgeL[2], SideThr));
  if (EdgeL.size() == 4)
    CC =(CatmullClarkability4(EdgeL[0], EdgeL[1], EdgeL[2], EdgeL[3], SideThr));
  if (EdgeL.size() == 5)
    CC = (CatmullClarkability5(EdgeL[0], EdgeL[1], EdgeL[2], EdgeL[3], EdgeL[4],
                               SideThr));
  if (EdgeL.size() == 6)
    CC = (CatmullClarkability6(EdgeL[0], EdgeL[1], EdgeL[2], EdgeL[3], EdgeL[4],
                               EdgeL[5], SideThr));

  // CC+=SideThr;

  //  if (CC<0)
  // {
  //     std::cout<<"*** CC Computation ***"<<CC<<std::endl;
  //     std::cout<<"Side "<<EdgeL.size()<<std::endl;
  //     for (size_t i=0;i<EdgeL.size();i++)
  //         std::cout<<EdgeL[i]<<",";
  //     std::cout<<std::endl;

  //     std::cout<<"SideThr: "<<SideThr<<std::endl;
  //     std::cout<<"CC: "<<CC<<std::endl;
  // }
  // CC=std::max(CC,(ScalarType)0);
  //    if ((EdgeL.size()==4)&&(CC>Ratio4))
  //        return true;
  //    if ((EdgeL.size()==4)&&(CC<=Ratio4))
  //        return false;

  return (CC > 0);
}

// template <class ScalarType>
// size_t AddedSingularities(//const int NumF,
//                           const std::vector<ScalarType> &EdgeL,
//                           const ScalarType &SideThr,
//                           bool SkipValence4,
//                           bool print_debug=false)
//{
//     assert(SideThr>0);
//     //bool
//     CC=IsCatmullClarkable(NumF,EdgeL,SideThr,SkipValence4,print_debug); bool
//     CC=IsCatmullClarkable(EdgeL,SideThr,SkipValence4,print_debug);
//     //std::cout<<"test CC "<<CC<<std::endl;
//     if (CC)
//         return 0;
//     else
//         return 1;
// }

// template <class ScalarType>
// void CCSideSubds_6(const ScalarType &e0, const ScalarType &e1,
//                  const ScalarType &e2, const ScalarType &e3, 
//                  const ScalarType &e4,const ScalarType &e5,
//                  ScalarType &s0,ScalarType &s1,
//                  ScalarType &s2, ScalarType &s3, 
//                  ScalarType &s4,ScalarType &s5) {
//   Geo::LinearConstraint<ScalarType> Constr0;
//   Constr0.Weight = 1;
//   // s1+s5=e0
//   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
//   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(5, 0, 1));
//   Constr0.RHSide = e0;

//   Geo::LinearConstraint<ScalarType> Constr1;
//   Constr1.Weight = 1;
//   // s0+s2=e1
//   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
//   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
//   Constr1.RHSide = e1;

//   Geo::LinearConstraint<ScalarType> Constr2;
//   Constr2.Weight = 1;
//   // s1+s3=e2
//   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
//   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(3, 0, 1));
//   Constr2.RHSide = e2;

//   Geo::LinearConstraint<ScalarType> Constr3;
//   Constr3.Weight = 1;
//   // s2+s4=e3
//   Constr3.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
//   Constr3.LHSide.push_back(Geo::VarTerm<ScalarType>(4, 0, 1));
//   Constr3.RHSide = e3;

//   Geo::LinearConstraint<ScalarType> Constr4;
//   Constr4.Weight = 1;
//   // s3+s5=e4
//   Constr4.LHSide.push_back(Geo::VarTerm<ScalarType>(3, 0, 1));
//   Constr4.LHSide.push_back(Geo::VarTerm<ScalarType>(5, 0, 1));
//   Constr4.RHSide = e4;
 
//   Geo::LinearConstraint<ScalarType> Constr5;
//   Constr5.Weight = 1;
//   // s0+s4=e5
//   Constr5.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
//   Constr5.LHSide.push_back(Geo::VarTerm<ScalarType>(4, 0, 1));
//   Constr5.RHSide = e5;

//   std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;
//   LinConstr.push_back(Constr0);
//   LinConstr.push_back(Constr1);
//   LinConstr.push_back(Constr2);
//   LinConstr.push_back(Constr3);
//   LinConstr.push_back(Constr4);
//   LinConstr.push_back(Constr5);

//   Geo::LinearConstraintLeastSQSolver<ScalarType> LQSolver;
//   LQSolver.Init(LinConstr);
//   std::vector<ScalarType> Solution;
//   LQSolver.Solve(Solution);

//   assert(Solution.size() == 6);
//   s0 = Solution[0];
//   s1 = Solution[1];
//   s2 = Solution[2];
//   s3 = Solution[3];
//   s4 = Solution[4];
//   s5 = Solution[5];

// }

// template <class ScalarType>
// void CCSideSubds_5(const ScalarType &e0, const ScalarType &e1,
//                  const ScalarType &e2, const ScalarType &e3, 
//                  const ScalarType &e4,
//                  ScalarType &s0,ScalarType &s1,
//                  ScalarType &s2, ScalarType &s3, 
//                  ScalarType &s4) {
//   Geo::LinearConstraint<ScalarType> Constr0;
//   Constr0.Weight = 1;
//   // s1+s4=e0
//   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
//   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(4, 0, 1));
//   Constr0.RHSide = e0;

//   Geo::LinearConstraint<ScalarType> Constr1;
//   Constr1.Weight = 1;
//   // s0+s2=e1
//   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
//   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
//   Constr1.RHSide = e1;

//   Geo::LinearConstraint<ScalarType> Constr2;
//   Constr2.Weight = 1;
//   // s1+s3=e2
//   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
//   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(3, 0, 1));
//   Constr2.RHSide = e2;

//   Geo::LinearConstraint<ScalarType> Constr3;
//   Constr3.Weight = 1;
//   // s2+s4=e3
//   Constr3.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
//   Constr3.LHSide.push_back(Geo::VarTerm<ScalarType>(4, 0, 1));
//   Constr3.RHSide = e3;

//   Geo::LinearConstraint<ScalarType> Constr4;
//   Constr4.Weight = 1;
//   // s0+s3=e4
//   Constr4.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
//   Constr4.LHSide.push_back(Geo::VarTerm<ScalarType>(3, 0, 1));
//   Constr4.RHSide = e4;

//   std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;
//   LinConstr.push_back(Constr0);
//   LinConstr.push_back(Constr1);
//   LinConstr.push_back(Constr2);
//   LinConstr.push_back(Constr3);
//   LinConstr.push_back(Constr4);

//   Geo::LinearConstraintLeastSQSolver<ScalarType> LQSolver;
//   LQSolver.Init(LinConstr);
//   std::vector<ScalarType> Solution;
//   LQSolver.Solve(Solution);

//   assert(Solution.size() == 5);
//   s0 = Solution[0];
//   s1 = Solution[1];
//   s2 = Solution[2];
//   s3 = Solution[3];
//   s4 = Solution[4];

// }

// template <class ScalarType>
// void CCSideSubds_4(const ScalarType &e0, const ScalarType &e1,
//                  const ScalarType &e2, const ScalarType &e3, 
//                  ScalarType &s0,ScalarType &s1, 
//                  ScalarType &s2, ScalarType &s3) {
// s0=(e1+e3)/4;
// s1=(e0+e2)/4;
// s3=s1;
// s2=s0;
// //   Geo::LinearConstraint<ScalarType> Constr0;
// //   Constr0.Weight = 1;
// //   // s1+s3=e0
// //   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
// //   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(3, 0, 1));
// //   Constr0.RHSide = e0;

// //   Geo::LinearConstraint<ScalarType> Constr1;
// //   Constr1.Weight = 1;
// //   // s0+s2=e1
// //   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
// //   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
// //   Constr1.RHSide = e1;

// //   Geo::LinearConstraint<ScalarType> Constr2;
// //   Constr2.Weight = 1;
// //   // s1+s3=e2
// //   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
// //   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(3, 0, 1));
// //   Constr2.RHSide = e2;

// //   Geo::LinearConstraint<ScalarType> Constr3;
// //   Constr3.Weight = 1;
// //   // s0+s2=e3
// //   Constr3.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
// //   Constr3.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
// //   Constr3.RHSide = e3;
// //   //A
// //   std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;
// //   LinConstr.push_back(Constr0);
// //   LinConstr.push_back(Constr1);
// //   LinConstr.push_back(Constr2);
// //   LinConstr.push_back(Constr3);

// //   Geo::LinearConstraintLeastSQSolver<ScalarType> LQSolver;
// //   LQSolver.Init(LinConstr);
// //   std::vector<ScalarType> Solution;
// //   LQSolver.Solve(Solution);

// //   assert(Solution.size() == 4);
// //   s0 = Solution[0];
// //   s1 = Solution[1];
// //   s2 = Solution[2];
// //   s3 = Solution[3];
// }

// template <class ScalarType>
// void CCSideSubds_3(const ScalarType &e0, const ScalarType &e1,
//                  const ScalarType &e2, ScalarType &s0, ScalarType &s1,
//                  ScalarType &s2) {
//   Geo::LinearConstraint<ScalarType> Constr0;
//   Constr0.Weight = 1;
//   // s1+s2=e0
//   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
//   Constr0.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
//   Constr0.RHSide = e0;

//   Geo::LinearConstraint<ScalarType> Constr1;
//   Constr1.Weight = 1;
//   // s0+s2=e1
//   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
//   Constr1.LHSide.push_back(Geo::VarTerm<ScalarType>(2, 0, 1));
//   Constr1.RHSide = e1;

//   Geo::LinearConstraint<ScalarType> Constr2;
//   Constr2.Weight = 1;
//   // s0+s1=e2
//   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(0, 0, 1));
//   Constr2.LHSide.push_back(Geo::VarTerm<ScalarType>(1, 0, 1));
//   Constr2.RHSide = e2;

//   std::vector<Geo::LinearConstraint<ScalarType>> LinConstr;
//   LinConstr.push_back(Constr0);
//   LinConstr.push_back(Constr1);
//   LinConstr.push_back(Constr2);

//   Geo::LinearConstraintLeastSQSolver<ScalarType> LQSolver;
//   LQSolver.Init(LinConstr);
//   std::vector<ScalarType> Solution;
//   LQSolver.Solve(Solution);

//   assert(Solution.size() == 3);
//   s0 = Solution[0];
//   s1 = Solution[1];
//   s2 = Solution[2];
// }

// template <class ScalarType>
// void CCSideSubds_2(const ScalarType &e0, const ScalarType &e1, ScalarType &s0,
//                  ScalarType &s1) {
//   s1 = e0 / 2;
//   s0 = e1 / 2;
// }

// template <class ScalarType>
// bool CCSideSubd(const std::vector<ScalarType> &EdgeL,
//               std::vector<ScalarType> &SideSubd) {
//   if (EdgeL.size() == 2) {

//     SideSubd.clear();
//     SideSubd.resize(2, 0);
//     CCSideSubds_2(EdgeL[0], EdgeL[1], SideSubd[0], SideSubd[1]);
//     return true;
//   }

//   if (EdgeL.size() == 3) {

//     SideSubd.clear();
//     SideSubd.resize(3, 0);

//     CCSideSubds_3(EdgeL[0], EdgeL[1], EdgeL[2], SideSubd[0], SideSubd[1],
//                 SideSubd[2]);
//     return true;
//   }

//   if (EdgeL.size() == 4) {

//     SideSubd.clear();
//     SideSubd.resize(4, 0);
//     CCSideSubds_4(EdgeL[0], EdgeL[1], EdgeL[2], EdgeL[3], SideSubd[0],
//                 SideSubd[1], SideSubd[2], SideSubd[3]);
//     return true;
//   }

//   if (EdgeL.size() == 5) {
//     //std::cout<<"Size 5"<<std::endl;

//     SideSubd.clear();
//     SideSubd.resize(5, 0);
//     CCSideSubds_5(EdgeL[0], EdgeL[1], 
//                 EdgeL[2], EdgeL[3],  EdgeL[4],
//                 SideSubd[0],SideSubd[1], 
//                 SideSubd[2], SideSubd[3],SideSubd[4]);
//     return true;
//   }

//   if (EdgeL.size() == 6) {
//     //std::cout<<"Size 6"<<std::endl;

//     SideSubd.clear();
//     SideSubd.resize(6, 0);
//     CCSideSubds_6(EdgeL[0], EdgeL[1], 
//                 EdgeL[2], EdgeL[3], 
//                 EdgeL[4], EdgeL[5],
//                 SideSubd[0],SideSubd[1], 
//                 SideSubd[2], SideSubd[3],
//                 SideSubd[4],SideSubd[5]);
//     return true;
//   }

//   return false;
// }


// template <class ScalarType>
// bool CCSideSubdRatio(const std::vector<ScalarType> &EdgeL,
//                     std::vector<ScalarType> &SideSubRatio) {
//     std::vector<ScalarType> SideSubd;
//     bool found=CCSideSubd(EdgeL,SideSubd);
//     if (!found)return false;
//     SideSubRatio.clear();
//     if (EdgeL.size() == 2) {
//         //s1,s0
//         // SideSubRatio.push_back(SideSubd[1]/EdgeL[1]);
//         // SideSubRatio.push_back(SideSubd[0]/EdgeL[0]);
//         SideSubRatio.push_back(0.5);
//         SideSubRatio.push_back(0.5);
//         return true;
//     }
//     if (EdgeL.size() == 3) {
//         //s2,s0,s1
//         // SideSubRatio.push_back(SideSubd[2]/EdgeL[2]);
//         // SideSubRatio.push_back(SideSubd[0]/EdgeL[0]);
//         // SideSubRatio.push_back(SideSubd[1]/EdgeL[1]);
        
//         ScalarType s0= (-EdgeL[0]+EdgeL[1]+EdgeL[2])/2;
//         ScalarType s1= (EdgeL[0]-EdgeL[1]+EdgeL[2])/2;
//         ScalarType s2= (EdgeL[0]+EdgeL[1]-EdgeL[2])/2;
//         SideSubRatio.push_back(s2/EdgeL[0]);
//         SideSubRatio.push_back(s0/EdgeL[1]);
//         SideSubRatio.push_back(s1/EdgeL[2]);
//         return true;
//     }
//     if (EdgeL.size() == 4) {
//         //s3,s0,s1,s2
//         // SideSubRatio.push_back(SideSubd[3]/EdgeL[3]);
//         // SideSubRatio.push_back(SideSubd[0]/EdgeL[0]);
//         // SideSubRatio.push_back(SideSubd[1]/EdgeL[1]);
//         // SideSubRatio.push_back(SideSubd[2]/EdgeL[2]);
//         SideSubRatio.push_back(0.5);
//         SideSubRatio.push_back(0.5);
//         SideSubRatio.push_back(0.5);
//         SideSubRatio.push_back(0.5);
//         return false;
//     }
//     if (EdgeL.size() == 5) {
//         //s4,s0,s1,s2,s3
//         ScalarType s0= (EdgeL[0]+EdgeL[1]-EdgeL[2]-EdgeL[3]+EdgeL[4])/2;
//         ScalarType s1= (EdgeL[0]+EdgeL[1]+EdgeL[2]-EdgeL[3]-EdgeL[4])/2;
//         ScalarType s2= (-EdgeL[0]+EdgeL[1]+EdgeL[2]+EdgeL[3]-EdgeL[4])/2;
//         ScalarType s3= (-EdgeL[0]-EdgeL[1]+EdgeL[2]+EdgeL[3]+EdgeL[4])/2;
//         ScalarType s4= (EdgeL[0]-EdgeL[1]-EdgeL[2]+EdgeL[3]+EdgeL[4])/2;
//         SideSubRatio.push_back(s4/EdgeL[0]);
//         SideSubRatio.push_back(s0/EdgeL[1]);
//         SideSubRatio.push_back(s1/EdgeL[2]);
//         SideSubRatio.push_back(s2/EdgeL[3]);
//         SideSubRatio.push_back(s3/EdgeL[4]);
//         // SideSubRatio.push_back(SideSubd[4]/EdgeL[4]);
//         // SideSubRatio.push_back(SideSubd[0]/EdgeL[0]);
//         // SideSubRatio.push_back(SideSubd[1]/EdgeL[1]);
//         // SideSubRatio.push_back(SideSubd[2]/EdgeL[2]);
//         // SideSubRatio.push_back(SideSubd[3]/EdgeL[3]);
//         return true;
//     }
//     if (EdgeL.size() == 6) {
//         //s5,s0,s1,s2,s3,s4
//         ScalarType s0= (EdgeL[1]-EdgeL[3]+EdgeL[5])/2;
//         ScalarType s1= (EdgeL[0]+EdgeL[2]-EdgeL[4])/2;
//         ScalarType s2= (EdgeL[1]+EdgeL[3]-EdgeL[5])/2;
//         ScalarType s3= (-EdgeL[0]+EdgeL[2]+EdgeL[4])/2;
//         ScalarType s4= (-EdgeL[1]+EdgeL[3]+EdgeL[5])/2;
//         ScalarType s5= (EdgeL[0]-EdgeL[2]+EdgeL[4])/2;
//         SideSubRatio.push_back(s5/EdgeL[0]);
//         SideSubRatio.push_back(s0/EdgeL[1]);
//         SideSubRatio.push_back(s1/EdgeL[2]);
//         SideSubRatio.push_back(s2/EdgeL[3]);
//         SideSubRatio.push_back(s3/EdgeL[4]);
//         SideSubRatio.push_back(s4/EdgeL[5]);

//         // SideSubRatio.push_back(SideSubd[5]/EdgeL[5]);
//         // SideSubRatio.push_back(SideSubd[0]/EdgeL[0]);
//         // SideSubRatio.push_back(SideSubd[1]/EdgeL[1]);
//         // SideSubRatio.push_back(SideSubd[2]/EdgeL[2]);
//         // SideSubRatio.push_back(SideSubd[3]/EdgeL[3]);
//         // SideSubRatio.push_back(SideSubd[4]/EdgeL[4]);
//         return true;
//     }
//     return false;
// }

// template <class ScalarType>
// void MakeCCable(std::vector<ScalarType> &EdgeL)
// {

// }

// NEW (approximate CC split): true if all the ratios are in
// [MinRatio, 1-MinRatio]
template <class ScalarType>
bool InRatioRange(const std::vector<ScalarType> &Ratio,
                  const ScalarType MinRatio) {
  for (size_t k = 0; k < Ratio.size(); k++)
    if ((Ratio[k] < MinRatio) || (Ratio[k] > ((ScalarType)1 - MinRatio)))
      return false;
  return true;
}

// NEW (approximate CC split): when the patch is not exactly Catmull-Clarkable
// some exact ratios fall outside [0,1] and the split edge was snapped to a
// corner (degenerate quads). Here the split lengths x_k (from the start corner
// of side k) are fitted in least squares to the CC condition
// x_{k+1} = e_{k-1} - x_{k-1}, keeping every ratio in [MinRatio, 1-MinRatio]
// (projected Gauss-Seidel, the system has at most 6 unknowns). Exactly
// CC-able patches keep the exact ratios.
template <class ScalarType>
void CCSideSubdRatioFit(const std::vector<ScalarType> &EdgeL,
                        std::vector<ScalarType> &SideSubRatio,
                        const ScalarType MinRatio, const int NumIte = 200,
                        const ScalarType RegW = (ScalarType)0.1) {
  int n = (int)EdgeL.size();
  assert((int)SideSubRatio.size() == n);
  // start from the exact ratios clamped into the bounds
  std::vector<ScalarType> X0(n), X(n), Low(n), Up(n);
  for (int k = 0; k < n; k++) {
    Low[k] = MinRatio * EdgeL[k];
    Up[k] = ((ScalarType)1 - MinRatio) * EdgeL[k];
    X0[k] = std::min(std::max(SideSubRatio[k] * EdgeL[k], Low[k]), Up[k]);
  }
  X = X0;
  for (int it = 0; it < NumIte; it++)
    for (int k = 0; k < n; k++) {
      // x_k appears in x_k = e_{k-2} - x_{k-2} and x_k = e_k - x_{k+2}
      int km2 = (k - 2 + n) % n;
      int kp2 = (k + 2) % n;
      ScalarType A = EdgeL[km2] - X[km2];
      ScalarType B = EdgeL[k] - X[kp2];
      ScalarType V = (A + B + RegW * X0[k]) / ((ScalarType)2 + RegW);
      X[k] = std::min(std::max(V, Low[k]), Up[k]);
    }
  for (int k = 0; k < n; k++)
    SideSubRatio[k] = X[k] / EdgeL[k];
}

template <class ScalarType>
bool CCSideSubdRatio(const std::vector<ScalarType> &EdgeL,
                    std::vector<ScalarType> &SideSubRatio,
                    // NEW (approximate CC split): min distance of a split from
                    // a corner, as a fraction of the side (<= 0 disables)
                    const ScalarType MinRatio = (ScalarType)0.1) {
    std::vector<ScalarType> SideSubd;
   
    SideSubRatio.clear();
    // NEW (approximate CC split): zero length sides give no valid ratio
    for (size_t k = 0; k < EdgeL.size(); k++)
      if (EdgeL[k] <= 0) {
        SideSubRatio.assign(EdgeL.size(), (ScalarType)0.5);
        return false;
      }
    if (EdgeL.size() == 2) {
        SideSubRatio.push_back(0.5);
        SideSubRatio.push_back(0.5);
        return true;
    }
    if (EdgeL.size() == 3) {
        //s2,s0,s1
        ScalarType s0= (-EdgeL[0]+EdgeL[1]+EdgeL[2])/2;
        ScalarType s1= (EdgeL[0]-EdgeL[1]+EdgeL[2])/2;
        ScalarType s2= (EdgeL[0]+EdgeL[1]-EdgeL[2])/2;
        SideSubRatio.push_back(s2/EdgeL[0]);
        SideSubRatio.push_back(s0/EdgeL[1]);
        SideSubRatio.push_back(s1/EdgeL[2]);
        // NEW (approximate CC split): fit when some ratio is out of range
        if ((MinRatio > 0) && (!InRatioRange(SideSubRatio, MinRatio)))
          CCSideSubdRatioFit(EdgeL, SideSubRatio, MinRatio);
        return true;
    }
    if (EdgeL.size() == 4) {
        //s3,s0,s1,s2
        SideSubRatio.push_back(0.5);
        SideSubRatio.push_back(0.5);
        SideSubRatio.push_back(0.5);
        SideSubRatio.push_back(0.5);
        return false;
    }
    if (EdgeL.size() == 5) {
        //s4,s0,s1,s2,s3
        ScalarType s0= (EdgeL[0]+EdgeL[1]-EdgeL[2]-EdgeL[3]+EdgeL[4])/2;
        ScalarType s1= (EdgeL[0]+EdgeL[1]+EdgeL[2]-EdgeL[3]-EdgeL[4])/2;
        ScalarType s2= (-EdgeL[0]+EdgeL[1]+EdgeL[2]+EdgeL[3]-EdgeL[4])/2;
        ScalarType s3= (-EdgeL[0]-EdgeL[1]+EdgeL[2]+EdgeL[3]+EdgeL[4])/2;
        ScalarType s4= (EdgeL[0]-EdgeL[1]-EdgeL[2]+EdgeL[3]+EdgeL[4])/2;
        SideSubRatio.push_back(s4/EdgeL[0]);
        SideSubRatio.push_back(s0/EdgeL[1]);
        SideSubRatio.push_back(s1/EdgeL[2]);
        SideSubRatio.push_back(s2/EdgeL[3]);
        SideSubRatio.push_back(s3/EdgeL[4]);
        // NEW (approximate CC split): fit when some ratio is out of range
        if ((MinRatio > 0) && (!InRatioRange(SideSubRatio, MinRatio)))
          CCSideSubdRatioFit(EdgeL, SideSubRatio, MinRatio);
        return true;
    }
    if (EdgeL.size() == 6) {
        //s5,s0,s1,s2,s3,s4
        ScalarType s0= (EdgeL[1]-EdgeL[3]+EdgeL[5])/2;
        ScalarType s1= (EdgeL[0]+EdgeL[2]-EdgeL[4])/2;
        ScalarType s2= (EdgeL[1]+EdgeL[3]-EdgeL[5])/2;
        ScalarType s3= (-EdgeL[0]+EdgeL[2]+EdgeL[4])/2;
        ScalarType s4= (-EdgeL[1]+EdgeL[3]+EdgeL[5])/2;
        ScalarType s5= (EdgeL[0]-EdgeL[2]+EdgeL[4])/2;
        SideSubRatio.push_back(s5/EdgeL[0]);
        SideSubRatio.push_back(s0/EdgeL[1]);
        SideSubRatio.push_back(s1/EdgeL[2]);
        SideSubRatio.push_back(s2/EdgeL[3]);
        SideSubRatio.push_back(s3/EdgeL[4]);
        SideSubRatio.push_back(s4/EdgeL[5]);
        // NEW (approximate CC split): fit when some ratio is out of range
        if ((MinRatio > 0) && (!InRatioRange(SideSubRatio, MinRatio)))
          CCSideSubdRatioFit(EdgeL, SideSubRatio, MinRatio);
        return true;
    }
    return false;
}
#endif
