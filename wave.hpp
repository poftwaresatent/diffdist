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

#ifndef DIFFDIST_WAVE_HPP
#define DIFFDIST_WAVE_HPP

#include "pose.hpp"


namespace diffdist {
  
  
  inline Pose computeWaveSolution1(double base, double mu, double t) {
    double const kappa((1.0 - 2.0 * mu) / base);
    if (fabs(kappa) < 1e-6) {
      Pose pp(0.5 * t, 0.0, 0.0);
      return pp;
    }
    Pose pp(1.0 / kappa, t * (1.0 - 2.0 * mu) / 2.0 / base);
    return pp;
  }
  
  
  inline Pose computeWaveSolution2(double base, double mu, double t) {
    Pose pp(base * (1.0 - 2.0 * mu), t / 2.0 / base);
    return pp;
  }
  
  
  inline Pose computeWaveSolution3(double base, double mu, double t) {
    Pose pp(base * (2.0 * mu - 1.0), -t / 2.0 / base);
    return pp;
  }
  
  
  inline Pose computeWaveSolution4(double base, double mu, double t) {
    double const kappa((1.0 - 2.0 * mu) / base);
    if (fabs(kappa) < 1e-6) {
      Pose pp(-0.5 * t, 0.0, 0.0);
      return pp;
    }
    Pose pp(1.0 / kappa, t * (2.0 * mu - 1.0) / 2.0 / base);
    return pp;
  }
  
  
  inline Pose computeWave(double base, int segment, double mu, double t) {
    switch (segment % 4) {
    case 0:
      return computeWaveSolution1(base, mu, t);
    case 1:
      return computeWaveSolution3(base, 1.0 - mu, t);
    case 2:
      return computeWaveSolution4(base, mu, t);
      /* default:
	 fall through */
    }
    return computeWaveSolution2(base, 1.0 - mu, t);
  }
  
}

#endif // DIFFDIST_WAVE_HPP
