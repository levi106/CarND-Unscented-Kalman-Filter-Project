#include "logger.h"
#include <fstream>

using namespace std;

class Logger::impl {
public:
  impl();
  ~impl();

  void Init();
  void Log(const MeasurementPackage& meas_package, const VectorXd& estimate, const VectorXd& gt_values, const VectorXd& rmse);

private:
  static const string filename_;
  ofstream os_;
};

Logger::Logger() : pImpl{make_unique<impl>()} {
}

Logger::~Logger() = default;

void Logger::Init() {
  pImpl->Init();
}

void Logger::Log(const MeasurementPackage& meas_package, const VectorXd& estimate, const VectorXd& gt_values, const VectorXd& rmse) {
  pImpl->Log(meas_package, estimate, gt_values, rmse);
}

const std::string Logger::impl::filename_{"ukf.csv"};

Logger::impl::impl() {
  os_.open(filename_, ios::trunc);
}

Logger::impl::~impl() {
  os_.close();
}

void Logger::impl::Init() {
  os_ << "sensor_type,timestamp,sensor_0,sensor_1,sensor_2,";
  os_ << "estimate_x,estimate_y,estimate_vx,estimate_vy,";
  os_ << "ground_truth_x,ground_truth_y,ground_truth_vx,ground_truth_vy,";
  os_ << "rmse_x,rmse_y,rmse_vx,rmse_vy" << endl;
}

void Logger::impl::Log(const MeasurementPackage& meas_package, const VectorXd& estimate, const VectorXd& gt_values, const VectorXd& rmse) {
  os_ << meas_package.sensor_type_ << ",";
  os_ << meas_package.timestamp_ << ",";
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    os_ << meas_package.raw_measurements_(0) << ",";
    os_ << meas_package.raw_measurements_(1) << ",,";
  } else {
    os_ << meas_package.raw_measurements_(0) << ",";
    os_ << meas_package.raw_measurements_(1) << ",";
    os_ << meas_package.raw_measurements_(2) << ",";    
  }
  os_ << estimate(0) << ",";
  os_ << estimate(1) << ",";
  os_ << estimate(2) << ",";
  os_ << estimate(3) << ",";
  os_ << gt_values(0) << ",";
  os_ << gt_values(1) << ",";
  os_ << gt_values(2) << ",";
  os_ << gt_values(3) << ",";
  os_ << rmse(0) << ",";
  os_ << rmse(1) << ",";
  os_ << rmse(2) << ",";
  os_ << rmse(3) << endl;
}
