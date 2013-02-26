/*
 * Copyright (c) 2013 Roland Philippsen. All rights reserved.
 *
 * BSD license:
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of
 *    contributors to this software may be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR THE CONTRIBUTORS TO THIS SOFTWARE BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <iostream>
#include <cmath>

using namespace std;


class Pose
{
public:
  double x_, y_, theta_;
  
  void rtr1(double B, double mu, double t) {
    double const kappa((1.0 - 2.0 * mu) / B);
    if (fabs(kappa) < 1e-6) {
      x_ = 0.5 * t;
      y_ = 0.0;
      theta_ = 0.0;
      return;
    }
    
    double const R(1.0 / kappa);
    double const bdot((1.0 - 2.0 * mu) / 2.0 / B);
    theta_ = bdot * t;
    x_ = R * sin(theta_);
    y_ = R * (1.0 - cos(theta_));
  }
  
  void rtr2(double B, double mu, double t) {
    double const R(B * (1.0 - 2.0 * mu));
    double const bdot(1.0 / 2.0 / B);
    theta_ = bdot * t;
    x_ = R * sin(theta_);
    y_ = R * (1.0 - cos(theta_));
  }
  
  void rtr3(double B, double mu, double t) {
    double const R(B * (2.0 * mu - 1.0));
    double const bdot(-1.0 / 2.0 / B);
    theta_ = bdot * t;
    x_ = R * sin(theta_);
    y_ = R * (1.0 - cos(theta_));
  }
  
  void rtr4(double B, double mu, double t) {
    double const kappa((1.0 - 2.0 * mu) / B);
    if (fabs(kappa) < 1e-6) {
      x_ = -0.5 * t;
      y_ = 0.0;
      theta_ = 0.0;
      return;
    }
    
    double const R(1.0 / kappa);
    double const bdot((2.0 * mu - 1.0) / 2.0 / B);
    theta_ = bdot * t;
    x_ = R * sin(theta_);
    y_ = R * (1.0 - cos(theta_));
  }
};


int main(int argc, char ** argv)
{
  double const B(1.0);
  size_t const nmu(20);
  double const tmax(3.0);
  size_t const nt(20);
  Pose pp;
  
  cout << "# case 1: t mu x y theta\n";
  for (size_t imu(0); imu <= nmu; ++imu) {
    double const mu((double) imu / nmu);
    for (size_t it(0); it <= nt; ++it) {
      double const t(tmax * it / nt);
      pp.rtr1(B, mu, t);
      cout << t << "  " << mu << "  " << pp.x_ << "  " << pp.y_ << "  " << pp.theta_ << "\n";
    }
    cout << "\n";
  }
  
  cout << "\n\n# case 2: t mu x y theta\n";
  for (size_t imu(0); imu <= nmu; ++imu) {
    double const mu((double) imu / nmu);
    for (size_t it(0); it <= nt; ++it) {
      double const t(tmax * it / nt);
      pp.rtr2(B, mu, t);
      cout << t << "  " << mu << "  " << pp.x_ << "  " << pp.y_ << "  " << pp.theta_ << "\n";
    }
    cout << "\n";
  }
  
  cout << "\n# case 3: t mu x y theta\n";
  for (size_t imu(0); imu <= nmu; ++imu) {
    double const mu((double) imu / nmu);
    for (size_t it(0); it <= nt; ++it) {
      double const t(tmax * it / nt);
      pp.rtr3(B, mu, t);
      cout << t << "  " << mu << "  " << pp.x_ << "  " << pp.y_ << "  " << pp.theta_ << "\n";
    }
    cout << "\n";
  }
  
  cout << "\n# case 4: t mu x y theta\n";
  for (size_t imu(0); imu <= nmu; ++imu) {
    double const mu((double) imu / nmu);
    for (size_t it(0); it <= nt; ++it) {
      double const t(tmax * it / nt);
      pp.rtr4(B, mu, t);
      cout << t << "  " << mu << "  " << pp.x_ << "  " << pp.y_ << "  " << pp.theta_ << "\n";
    }
    cout << "\n";
  }
  
  cout << "#  set view equal xy\n"
       << "#  set hidden3d\n"
       << "#  splot 'data' u 3:4:5 w l\n";
}
