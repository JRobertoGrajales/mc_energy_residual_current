/*
 * Copyright 2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_control/GlobalPlugin.h>
#include "LpfThreshold.h"

#include <RBDyn/Coriolis.h>
#include <RBDyn/FA.h>
#include <RBDyn/FK.h>
#include <RBDyn/FV.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

namespace mc_plugin
{

struct EnergyResidualCurrent : public mc_control::GlobalPlugin
{
  void init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config) override;

  void reset(mc_control::MCGlobalController & controller) override;

  void before(mc_control::MCGlobalController & controller) override;

  void after(mc_control::MCGlobalController & controller) override;

  void energy_residual_computation(mc_control::MCGlobalController & controller);

  void addGui(mc_control::MCGlobalController & controller);
  void addPlot(mc_control::MCGlobalController & controller);
  void addLog(mc_control::MCGlobalController & controller);

  mc_control::GlobalPlugin::GlobalPluginConfiguration configuration() override;

  ~EnergyResidualCurrent() override;

private:
  double dt_;
  double counter_;

  double ko;
  int jointNumber;
  rbd::Coriolis * coriolis;
  rbd::ForwardDynamics forwardDynamics;
  double gear_ratio = 100.0;
  std::map<std::string, double> kt;
  Eigen::VectorXd tau_fric;
  double integralTerm;
  double energy_input;
  double t_kinetic;
  double tzero; // kinetic energy init
  double residual;

  LpfThreshold lpf_threshold_;
  double threshold_offset_;
  double threshold_filtering_;
  double residual_high_;
  double residual_low_;

  bool activate_plot_ = false;
  bool plot_added_ = false;
  bool collision_stop_activated_ = false;
  bool collision_stop_activated_zurlo_ = false;
  bool obstacle_detected_ = false;
  bool activate_verbose = false;
};

} // namespace mc_plugin
