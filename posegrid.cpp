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

#include "posegrid.hpp"

namespace diffdist {
  
  Posegrid::index::
  index(size_t ix_, size_t iy_, size_t itheta_)
    : ix(ix_), iy(iy_), itheta(itheta_)
  {
  }
  
  
  Posegrid::
  Posegrid(double x0, double x1, size_t nx,
	   double y0, double y1, size_t ny,
	   size_t ntheta)
    : nx_(nx), ny_(ny), ntheta_(ntheta),
      x0_(x0), x1_(x1), dx_((x1_ - x0) / (nx - 1)),
      y0_(y0), y1_(y1), dy_((y1_ - y0) / (ny - 1)),
      dtheta_(2.0 * M_PI / ntheta) // theta wraps around, thus no "-1" here
  {
    pose_.resize(nx);
    for (size_t ix(0); ix < nx; ++ix) {
      pose_[ix].resize(ny);
      for (size_t iy(0); iy < ny; ++iy) {
	pose_[ix][iy].resize(ntheta);
	for (size_t itheta(0); itheta < ntheta; ++itheta) {
	  pose_[ix][iy][itheta] = new Pose(x0 + ix * dx_, y0 + iy * dy_, itheta * dtheta_);
	}
      }
    }
  }
  
  
  Posegrid::
  ~Posegrid()
  {
    for (size_t ix(0); ix < nx_; ++ix) {
      for (size_t iy(0); iy < ny_; ++iy) {
	for (size_t itheta(0); itheta < ntheta_; ++itheta) {
	  delete pose_[ix][iy][itheta];
	}
      }
    }
  }
  
  
  Posegrid::index Posegrid::
  snap(double xx, double yy, double theta) const
  {
    if (xx < x0_) {
      xx = x0_;
    }
    else if (xx > x1_) {
      xx = x1_;
    }

    if (yy < y0_) {
      yy = y0_;
    }
    else if (yy > y1_) {
      yy = y1_;
    }
    
    index idx(static_cast<size_t>(rint((xx - x0_) / dx_)),
	      static_cast<size_t>(rint((yy - y0_) / dy_)),
	      static_cast<size_t>(rint(normangle(theta) / dtheta_)));
    return idx;
  }
  
}
