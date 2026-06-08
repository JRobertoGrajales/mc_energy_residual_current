#include "EnergyResidualCurrent.h"

#include <mc_control/GlobalPluginMacros.h>

namespace mc_plugin
{

EnergyResidualCurrent::~EnergyResidualCurrent() = default;

void EnergyResidualCurrent::init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);

  auto & robot = ctl.robot(ctl.robots()[0].name());
  auto & realRobot = ctl.realRobot(ctl.robots()[0].name());
  auto & rjo = robot.refJointOrder();

  // Make sure to have obstacle detection
  if(!ctl.controller().datastore().has("Obstacle detected"))
  {
    ctl.controller().datastore().make<bool>("Obstacle detected", false);
  }

  ctl.controller().datastore().make<bool>("Energy Residual Current Obstacle detected", false);

  dt_ = ctl.timestep();
  counter_ = 0.0;
  jointNumber = ctl.robot(ctl.robots()[0].name()).refJointOrder().size();

  auto plugin_config = config("energy_residual");

  gear_ratio = plugin_config("gear_ratio", 100.0);
  ko = plugin_config("ko", 1.0);

  kt = plugin_config("kt");
  if(kt.empty())
  {
    kt = {
      {"joint_1", 0.11}, 
      {"joint_2", 0.11}, 
      {"joint_3", 0.11},
      {"joint_4", 0.11}, 
      {"joint_5", 0.076}, 
      {"joint_6", 0.076},
      {"joint_7", 0.076}
    };
  }

  threshold_filtering_ = plugin_config("threshold_filtering", 0.05);
  threshold_offset_ = plugin_config("threshold_offset", 2.0);
  lpf_threshold_.setValues(threshold_offset_, threshold_filtering_, jointNumber);

  if(ctl.controller().datastore().has("torque_fric"))
  {
    tau_fric = ctl.controller().datastore().get<Eigen::VectorXd>("torque_fric");
    for(int i = 0; i < jointNumber; i++)
    {
      mc_rtc::log::info("[EnergyResidualCurrent] Joint {} has friction torque {} N.m", rjo[i], tau_fric[i]);
    }
  }
  else
  {
    tau_fric.setZero(jointNumber);
    mc_rtc::log::error("[EnergyResidualCurrent] No torque_fric in datastore");
  }

  for(auto const& [key, val] : kt)
  {
    double tau_mot = val * gear_ratio * realRobot.jointJointSensor(key).motorCurrent();
    mc_rtc::log::info("[EnergyResidualCurrent] Joint {} has kt value {}, motor torque current based {} N.m", key, val, tau_mot);
  }

  Eigen::VectorXd qdot(jointNumber);
  for(size_t i = 0; i < jointNumber; i++)
  {
    qdot[i] = robot.alpha()[robot.jointIndexByName(rjo[i])][0];
  }
  residual = 0.0;
  residual_high_ = 0.0;
  residual_low_ = 0.0;
  integralTerm = 0.0;
  coriolis = new rbd::Coriolis(robot.mb());
  forwardDynamics = rbd::ForwardDynamics(robot.mb());

  forwardDynamics.computeH(robot.mb(), robot.mbc());
  auto inertiaMatrix = forwardDynamics.H();
  tzero = 0.5 * qdot.transpose() * inertiaMatrix * qdot;

  addGui(ctl);
  addLog(ctl);

  mc_rtc::log::info("EnergyResidualCurrent::init called with configuration:\n{}", config.dump(true, true));
}

void EnergyResidualCurrent::reset(mc_control::MCGlobalController & controller)
{
  // mc_rtc::log::info("EnergyResidualCurrent::reset called");
}

void EnergyResidualCurrent::before(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  counter_ += dt_;

  if(activate_plot_ && !plot_added_)
  {
    addPlot(controller);
    plot_added_ = true;
  }

  if(ctl.controller().datastore().has("Zurlo Collision Detection"))
  {
    collision_stop_activated_zurlo_ = ctl.controller().datastore().get<bool>("Zurlo Collision Detection");
  }

  energy_residual_computation(controller);
  residual_high_ = lpf_threshold_.adaptiveThreshold(residual, true);
  residual_low_ = lpf_threshold_.adaptiveThreshold(residual, false);
  obstacle_detected_ = false;
  if(residual > residual_high_ || residual < residual_low_)
  {
    obstacle_detected_ = true;
    if(activate_verbose) mc_rtc::log::info("[Energy Residual Current] Obstacle detected");
    if(collision_stop_activated_)
    {
      ctl.controller().datastore().get<bool>("Obstacle detected") = obstacle_detected_;
    }
  }
  ctl.controller().datastore().get<bool>("Energy Residual Current Obstacle detected") = obstacle_detected_;
}

void EnergyResidualCurrent::after(mc_control::MCGlobalController & controller)
{
  // mc_rtc::log::info("EnergyResidualCurrent::after");
}

mc_control::GlobalPlugin::GlobalPluginConfiguration EnergyResidualCurrent::configuration()
{
  mc_control::GlobalPlugin::GlobalPluginConfiguration out;
  out.should_run_before = true;
  out.should_run_after = false;
  out.should_always_run = false;
  return out;
}

void EnergyResidualCurrent::energy_residual_computation(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);

  auto & robot = ctl.robot(ctl.robots()[0].name());
  auto & realRobot = ctl.realRobot(ctl.robots()[0].name());
  auto & rjo = robot.refJointOrder();

  Eigen::VectorXd tau_m = Eigen::VectorXd::Zero(jointNumber);

  int jointIndex = 0;
  for(auto const& [key, val] : kt)
  {
      if (jointIndex >= 0 && jointIndex < jointNumber) {
          double tau_mot = val * gear_ratio * realRobot.jointJointSensor(key).motorCurrent();
          tau_m[jointIndex] = tau_mot;
      } else {
          mc_rtc::log::error("[EnergyResidualCurrent] Invalid joint name: {} or index out of bounds", key);
      }
      jointIndex++;
  }

  Eigen::VectorXd qdot(jointNumber);
  rbd::paramToVector(realRobot.alpha(), qdot);
 
  forwardDynamics.computeC(realRobot.mb(), realRobot.mbc());
  forwardDynamics.computeH(realRobot.mb(), realRobot.mbc());
  auto coriolisMatrix = coriolis->coriolis(realRobot.mb(), realRobot.mbc());
  auto coriolisGravityTerm = forwardDynamics.C();
  auto negative_gravity = coriolisMatrix*qdot - coriolisGravityTerm;
  auto inertiaMatrix = forwardDynamics.H();
  double t_kinetic = 0.5 * qdot.transpose() * inertiaMatrix * qdot;

  integralTerm += (qdot.transpose()*(tau_m + negative_gravity - tau_fric) + residual) * ctl.timestep();
  residual = ko * (t_kinetic - integralTerm - tzero);

  if(!ctl.controller().datastore().has("energy_residual"))
  {
    ctl.controller().datastore().make<double>("energy_residual", residual);
  }
  else
  {
    ctl.controller().datastore().assign("energy_residual", residual);
  }
}

void EnergyResidualCurrent::addGui(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);

  ctl.controller().gui()->addElement({"Plugins", "EnergyResidualCurrent"},
  mc_rtc::gui::NumberInput(
      "Gain", [this]() { return this->ko; },
      [this](double gain)
      {
      this->integralTerm = 0.0;
      this->residual = 0.0;
      this->ko = gain;
      }),
      mc_rtc::gui::Button("Add plot", [this]() { return activate_plot_ = true; }),
    // Add checkbox to activate the collision stop
    mc_rtc::gui::Checkbox("Collision stop", collision_stop_activated_),
    mc_rtc::gui::Checkbox("Verbose", activate_verbose), 
    // Add Threshold offset input
    mc_rtc::gui::NumberInput("Threshold offset", [this](){return this->threshold_offset_;},
        [this](double offset)
      { 
        threshold_offset_ = offset;
        lpf_threshold_.setOffset(threshold_offset_); 
      }),
    // Add Threshold filtering input
    mc_rtc::gui::NumberInput("Threshold filtering", [this](){return this->threshold_filtering_;},
        [this](double filtering)
      { 
        threshold_filtering_ = filtering;
        lpf_threshold_.setFiltering(threshold_filtering_); 
      })       
    );

}

void EnergyResidualCurrent::addLog(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);

  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_residual", [this]() -> const double & { return residual; });
  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_threshold_high", [this]() -> const double & { return residual_high_; });
  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_threshold_low", [this]() -> const double & { return residual_low_; });
  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_threshold_offset", [this]() -> const double & { return threshold_offset_; });
  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_threshold_filtering", [this]() -> const double & { return threshold_filtering_; });
  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_ko", [this]() -> const double & { return ko; });
  ctl.controller().logger().addLogEntry("EnergyResidualCurrent_obstacleDetected", [this]() -> const bool & { return obstacle_detected_; });
}

void EnergyResidualCurrent::addPlot(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  auto & gui = *ctl.controller().gui();

  gui.addPlot(
    "EnergyResidualCurrent",
    mc_rtc::gui::plot::X("t", [this]() { return counter_; }),
    mc_rtc::gui::plot::Y("Residual", [this]() { return residual; }, mc_rtc::gui::Color::Red),
    mc_rtc::gui::plot::Y("Threshold high", [this]() { return residual_high_; }, mc_rtc::gui::Color::Green),
    mc_rtc::gui::plot::Y("Threshold low", [this]() { return residual_low_; }, mc_rtc::gui::Color::Blue)
  );
}

} // namespace mc_plugin

EXPORT_MC_RTC_PLUGIN("EnergyResidualCurrent", mc_plugin::EnergyResidualCurrent)