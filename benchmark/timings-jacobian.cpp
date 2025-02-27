//
// Copyright (c) 2015-2018 CNRS
//

#include "pinocchio/algorithm/frames.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/multibody/sample-models.hpp"

#include <iostream>

#include "pinocchio/utils/timer.hpp"

int main(int argc, const char ** argv)
{
  using namespace Eigen;
  using namespace pinocchio;

  PinocchioTicToc timer(PinocchioTicToc::US);
#ifdef NDEBUG
  const int NBT = 1000 * 2000;
#else
  const int NBT = 1;
  std::cout << "(the time score in debug mode is not relevant) " << std::endl;
#endif

  pinocchio::Model model;

  std::string filename = PINOCCHIO_MODEL_DIR + std::string("/simple_humanoid.urdf");
  if (argc > 1)
    filename = argv[1];

  bool with_ff = true;
  if (argc > 2)
  {
    const std::string ff_option = argv[2];
    if (ff_option == "-no-ff")
      with_ff = false;
  }

  if (filename == "HS")
    pinocchio::buildModels::humanoidRandom(model, true);
  else if (with_ff)
    pinocchio::urdf::buildModel(filename, JointModelFreeFlyer(), model);
  else
    pinocchio::urdf::buildModel(filename, model);
  std::cout << "nq = " << model.nq << std::endl;

  pinocchio::Model::JointIndex JOINT_ID = (Model::JointIndex)(model.njoints - 1);
  const SE3 & framePlacement = SE3::Random();
  auto FRAME_ID = model.addFrame(Frame("test_frame", JOINT_ID, 0, framePlacement, OP_FRAME));

  pinocchio::Data data(model);
  pinocchio::Data::Matrix6x J(6, model.nv);
  J.setZero();

  VectorXd q;
  VectorXd qmax = Eigen::VectorXd::Ones(model.nq);
  q = randomConfiguration(model, -qmax, qmax);
  forwardKinematics(model, data, q);
  computeJointJacobians(model, data, q);
  updateFramePlacements(model, data);

  // timer.tic();
  // SMOOTH(NBT)
  // {
  //   forwardKinematics(model, data, qs[_smooth]);
  // }
  // std::cout << "Zero Order Kinematics = \t";
  // timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    computeJointJacobian(model, data, q, JOINT_ID, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    computeJointJacobian(model, data, q, JOINT_ID, J);
  }
  std::cout << "computeJointJacobian = \t\t";
  timer.toc(std::cout, NBT);

  // timer.tic();
  // SMOOTH(NBT)
  // {
  //   computeJointJacobians(model, data, qs[_smooth]);
  // }
  // std::cout << "computeJointJacobians(q) = \t";
  // timer.toc(std::cout, NBT);
  //
  // timer.tic();
  // SMOOTH(NBT)
  // {
  //   computeJointJacobians(model, data);
  // }
  // std::cout << "computeJointJacobians() = \t";
  // timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    getJointJacobian(model, data, JOINT_ID, LOCAL, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    getJointJacobian(model, data, JOINT_ID, LOCAL, J);
  }
  std::cout << "getJointJacobian(LOCAL) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    getJointJacobian(model, data, JOINT_ID, WORLD, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    getJointJacobian(model, data, JOINT_ID, WORLD, J);
  }
  std::cout << "getJointJacobian(WORLD) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    getJointJacobian(model, data, JOINT_ID, LOCAL_WORLD_ALIGNED, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    getJointJacobian(model, data, JOINT_ID, LOCAL_WORLD_ALIGNED, J);
  }
  std::cout << "getJointJacobian(LOCAL_WORLD_ALIGNED) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    getFrameJacobian(model, data, FRAME_ID, LOCAL, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    getFrameJacobian(model, data, FRAME_ID, LOCAL, J);
  }
  std::cout << "getFrameJacobian(LOCAL) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    getFrameJacobian(model, data, FRAME_ID, WORLD, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    getFrameJacobian(model, data, FRAME_ID, WORLD, J);
  }
  std::cout << "getFrameJacobian(WORLD) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    getFrameJacobian(model, data, FRAME_ID, LOCAL_WORLD_ALIGNED, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    getFrameJacobian(model, data, FRAME_ID, LOCAL_WORLD_ALIGNED, J);
  }
  std::cout << "getFrameJacobian(LOCAL_WORLD_ALIGNED) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    computeFrameJacobian(model, data, q, FRAME_ID, LOCAL, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    computeFrameJacobian(model, data, q, FRAME_ID, LOCAL, J);
  }
  std::cout << "computeFrameJacobian(LOCAL) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    computeFrameJacobian(model, data, q, FRAME_ID, WORLD, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    computeFrameJacobian(model, data, q, FRAME_ID, WORLD, J);
  }
  std::cout << "computeFrameJacobian(WORLD) = \t";
  timer.toc(std::cout, NBT);

  SMOOTH(NBT / 100)
  {
    computeFrameJacobian(model, data, q, FRAME_ID, LOCAL_WORLD_ALIGNED, J);
  }
  timer.tic();
  SMOOTH(NBT)
  {
    computeFrameJacobian(model, data, q, FRAME_ID, LOCAL_WORLD_ALIGNED, J);
  }
  std::cout << "computeFrameJacobian(LOCAL_WORLD_ALIGNED) = \t";
  timer.toc(std::cout, NBT);

  std::cout << "--" << std::endl;
  return 0;
}
