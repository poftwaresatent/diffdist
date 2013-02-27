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

#include "pose.hpp"
#include "wave.hpp"
#include <stdlib.h>
#include <err.h>
#include <stdio.h>
#include <signal.h>

using namespace diffdist;


static vector<vector<posevector_t> > poses;
static bool interrupt(false);


static void handle(int signal)
{
  if (interrupt) {
    errx(EXIT_FAILURE, "handle(%d)", signal);
  }
  fprintf(stderr, "# interrupting...\n");
  interrupt = true;
}


static void init_poses(double base, int nseg, int nmu, double tmax, int nt)
{
  poses.resize(nseg);
  for (int iseg(0); iseg < nseg; ++iseg) {
    poses[iseg].resize(nmu+1);
    for (int imu(0); imu <= nmu; ++imu) {
      double const mu((double) imu / nmu);
      for (int it(0); it <= nt; ++it) {
	double const tt(tmax * it / nt);
	poses[iseg][imu].push_back(computeWave(base, iseg, mu, tt));
      }
    }
  }
}


int main(int argc, char ** argv)
{
  if (SIG_ERR == signal(SIGINT, handle)) {
    err (EXIT_FAILURE, "signal");
  }
  
  printf("initializing poses\n");
  fflush(stdout);
  double base(1.25);
  int nseg(4);
  int nmu(2);
  double tmax(0.5);
  int nt(3);
  init_poses(base, nseg, nmu, tmax, nt);
  
  int nsteps(5);
  posevector_t sources;
  sources.push_back(Pose(0.0, 0.0, 0.0));
  
  for (int step(0); step < nsteps; ++step) {
    if (interrupt) {
      break;
    }
    
    char dfname[64];
    snprintf(dfname, 64, "bfl%03d.data", step);
    FILE * data(fopen(dfname, "w"));
    if (0 == data) {
      err(EXIT_FAILURE, "fopen bfl.data");
    }
    
    posevector_t next_sources;
    for (size_t isrc(0); isrc < sources.size(); ++isrc) {
      Pose const & src(sources[isrc]);
      printf("step %d  source %zu / %zu   %g %g %g\n",
	     step, isrc, sources.size(),
	     src.x(), src.y(), src.theta());
      
      for (int iseg(0); iseg < nseg; ++iseg) {
	for (int imu(0); imu <= nmu; ++imu) {
	  for (int it(0); it <= nt; ++it) {
	    if (interrupt) {
	      break;
	    }
	    
	    Pose pp(src + poses[iseg][imu][it]);
	    fprintf(data, "%8g  %8g  %8g\n", pp.x(), pp.y(), pp.theta());
	    if (nt == it) {
	      next_sources.push_back(pp);
	    }
	  }
	  fprintf(data, "\n");
	}
      }
      fprintf(data, "\n");
    }
    fclose(data);
    
    sources.swap(next_sources);
  }
  
  if (interrupt) {
    fprintf(stderr, "...interrupted\n");
  }
  
  FILE * plot(fopen("bfl.plot", "w"));
  if (0 == plot) {
    err(EXIT_FAILURE, "fopen bfl.plot");
  }
  fprintf(plot,
	  "set view equal xy\n"
	  "set view 56,210\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'y'\n"
	  "set zlabel 'theta'\n"
	  "set hidden3d\n"
	  "set term png\n");
  for (int step(0); step < nsteps; ++step) {
    fprintf(plot,
	    "set output 'bfl%03d.png'\n"
	    "splot 'bfl%03d.data' u 1:2:3 w l lc 0 t 'step %d'\n",
	    step, step, step);
  }
  fclose(plot);
  
  plot = fopen("bfl-anim.plot", "w");
  if (0 == plot) {
    err(EXIT_FAILURE, "fopen bfl-anim.plot");
  }
  fprintf(plot,
	  "set view equal xy\n"
	  "set xlabel 'x'\n"
	  "set ylabel 'y'\n"
	  "set zlabel 'theta'\n"
	  "set hidden3d\n"
	  "set term png\n");
  for (int deg(0); deg < 360; ++deg) {
    fprintf(plot,
	    "set view %d,210\n"
	    "set output 'bfl-anim-%03d.png'\n",
	    deg, deg);
    for (int step(0); step < nsteps; ++step) {
      if (0 == step) {
	fprintf(plot, "splot ");
      }
      else {
	fprintf(plot, ", ");
      }
      fprintf(plot, "'bfl%03d.data' u 1:2:3 w l lc %d notitle", step, step);
    }
    fprintf(plot, "\n");
  }
  fclose(plot);
}
