#ifndef LOGGER_H_
#define LOGGER_H_

#include <memory>
#include <experimental/propagate_const>
#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::VectorXd;

class Logger {
public:
  /**
  * Constructor.
  */
  Logger();

  /**
  * Destructor
  */
  virtual ~Logger();

  /**
  * Initialize Logger
  */
  void Init();

  /**
  * Output Log
  */
  void Log(const MeasurementPackage& meas_package, const VectorXd& estimate, const VectorXd& gt_values, const VectorXd& rmse);

private:
  class impl;
  std::experimental::propagate_const<std::unique_ptr<impl>> pImpl;
};

#endif /* LOGGER_H_ */
