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

#ifndef DIFFDIST_POSEGRID_HPP
#define DIFFDIST_POSEGRID_HPP

#include "pose.hpp"
#include <vector>

namespace diffdist {

  
  class Posegrid
  {
  public:
    struct index {
      index(size_t ix, size_t iy, size_t itheta);
      size_t ix, iy, itheta;
    };
    
    Posegrid(double x0, double x1, size_t nx,
	     double y0, double y1, size_t ny,
	     size_t ntheta);
    
    ~Posegrid();
    
    static inline double normangle(double theta) {
      theta = fmod(theta, 2.0 * M_PI);
      if (theta <= -M_PI) {
	return theta + 2.0 * M_PI;
      }
      if (theta > M_PI) {
	return theta - 2.0 * M_PI;
      }
      return theta;
    }
    
    Pose const & get(size_t ix, size_t iy, size_t itheta) const
    { return *pose_[ix][iy][itheta]; }
    
    Pose const & get(index const & idx) const
    { return *pose_[idx.ix][idx.iy][idx.itheta]; }
    
    index snap(double xx, double yy, double theta) const;
    
    inline index snap(Pose const & pose) const
    { return snap(pose.x(), pose.y(), pose.theta()); }
    
    inline Pose const & closest(Pose const & pose) const
    { return get(snap(pose.x(), pose.y(), pose.theta())); }
    
    size_t const nx_;
    size_t const ny_;
    size_t const ntheta_;
    
  private:
    double x0_, x1_, dx_;
    double y0_, y1_, dy_;
    double dtheta_;
    vector<vector<vector<Pose*> > > pose_;
  };
  
}

#endif // DIFFDIST_POSEGRID_HPP
