#include <iostream>
#include <cmath>
#include <pybind11/pybind11.h>
#include <string>
 
class EuropeanOption {
    public:
        EuropeanOption(double S, double K, double r, double T, double sigma)
            : S_(S), K_(K), r_(r), T_(T), sigma_(sigma) {
                if (S <= 0 || K <= 0 || r < 0 || T <= 0 || sigma <= 0) {
                    throw std::invalid_argument("Invalid parameters for European Option");
                }
                d1_ = (std::log(S_ / K_) + (r_ + 0.5 * sigma_ * sigma_) * T_) / (sigma_ * std::sqrt(T_));
                d2_ = d1_ - sigma_ * std::sqrt(T_);
            }
            
            double priceCall() const {
                return S_ * normalCDF(d1_) - K_ * std::exp(-r_ * T_) * normalCDF(d2_);
            }
            double pricePut() const {
                return K_ * std::exp(-r_ * T_) * normalCDF(-d2_) - S_ * normalCDF(-d1_);
            }
            //---Greeks---
            double deltaCall() const {
                return normalCDF(d1_);
            }
            double deltaPut() const {
                return normalCDF(d1_) - 1;}
            double gamma() const {
                return normalPDF(d1_) / (S_ * sigma_ * std::sqrt(T_));
            }
            double vega() const {
                return S_ * normalPDF(d1_) * std::sqrt(T_);
            }
            double thetaCall() const {
                return -(S_ * normalPDF(d1_) * sigma_) / (2 * std::sqrt(T_)) - r_ * K_ * std::exp(-r_ * T_) * normalCDF(d2_);
            }
            double thetaPut() const {
                return -(S_ * normalPDF(d1_) * sigma_) / (2 * std::sqrt(T_)) + r_ * K_ * std::exp(-r_ * T_) * normalCDF(-d2_);
            }   
            double rhoCall() const {
                return K_ * T_ * std::exp(-r_ * T_) * normalCDF(d2_);
            }
            double rhoPut() const {
                return -K_ * T_ * std::exp(-r_ * T_) * normalCDF(-d2_);
            }   
                    
    private:    
        double S_;  // Underlying asset price
        double K_;  // Strike price
        double r_;  // Risk-free interest rate
        double T_;  // Time to maturity
        double sigma_;  // Volatility
        double d1_;
        double d2_;

        double normalCDF(double x) const {
            return 0.5 * std::erfc(-x / std::sqrt(2));
        }
        double normalPDF(double x) const {
            return std::exp(-0.5 * x * x) / std::sqrt(2 * M_PI);
        }
};

namespace py = pybind11;

PYBIND11_MODULE(bs, m) {
    py::class_<EuropeanOption>(m, "EuropeanOption")
        .def(py::init<double, double, double, double, double>())
        .def("priceCall", &EuropeanOption::priceCall)
        .def("pricePut", &EuropeanOption::pricePut)
        .def("deltaCall", &EuropeanOption::deltaCall)
        .def("deltaPut", &EuropeanOption::deltaPut)
        .def("gamma", &EuropeanOption::gamma)
        .def("vega", &EuropeanOption::vega)
        .def("thetaCall", &EuropeanOption::thetaCall)
        .def("thetaPut", &EuropeanOption::thetaPut)
        .def("rhoCall", &EuropeanOption::rhoCall)
        .def("rhoPut", &EuropeanOption::rhoPut);
}
