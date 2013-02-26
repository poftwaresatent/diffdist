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
#include "wave.hpp"
#include <stdio.h>
#include <cmath>

using namespace diffdist;


int main(int argc, char ** argv)
{
  double b0(0.25);
  double b1(1.5);
  int nb(5);
  int nmu(20);
  double tmax(3.0);
  int nt(20);
  
  Posegrid grid(-1.0, 1.0, 20, -0.4, 0.4, 20, 20);
  
  printf("# base t segment mu x y theta\n");
  for (int ib(0); ib <= nb; ++ib) {
    double const base(b0 + ib * (b1 - b0) / nb);
    for (int iseg(0); iseg < 4; ++iseg) {
      for (int imu(0); imu <= nmu; ++imu) {
	double const mu((double) imu / nmu);
	for (int it(0); it <= nt; ++it) {
	  double const t(tmax * it / nt);
	  Pose const real(computeWave(base, iseg, mu, t));
	  Pose const & pp(grid.closest(real));
	  printf("%8g  %8g  %d  %8g    %8g  %8g  %8g\n", base, t, iseg, mu, pp.x(), pp.y(), pp.theta());
	}
	printf("\n");
      }
    }
    printf("\n");
  }
  
  printf("#  set view equal xy\n"
	 "#  set hidden3d\n"
	 "#  splot 'data' u 5:6:7 w l lc 0\n");
}
