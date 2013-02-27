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
#include "heap.h"
#include <stdlib.h>
#include <err.h>
#include <stdio.h>
#include <signal.h>

using namespace diffdist;


struct Sprite : public Pose {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  
  Sprite(double tt, Pose const & pose)
    : Pose(pose), tt_(tt)
  {}
  
  double tt_;
};


struct Entry {
  Entry()
    : pose_(0),
      tt_(numeric_limits<double>::max()),
      key_(-1.0)
  {}
  
  Pose const * pose_;
  double tt_, key_;
};


static heap_t * heap;
static vector<Sprite, Eigen::aligned_allocator<Sprite> > sprite;
static vector<vector<vector<Entry> > > lookup;
static size_t const nn(25);
static bool interrupt(false);


static void handle(int signal)
{
  if (interrupt) {
    errx(EXIT_FAILURE, "handle(%d)", signal);
  }
  fprintf(stderr, "# interrupting...\n");
  interrupt = true;
}


static void cleanup()
{
  if (heap) {
    heap_destroy(heap);
  }
}


static void init_sprite(double base, int nmu, double tmax, int nt)
{
  sprite.push_back(Sprite(0.0, computeWave(base, 0, 0.0, 0.0)));
  
  for (int iseg(0); iseg < 4; ++iseg) {
    for (int imu(0); imu <= nmu; ++imu) {
      double const mu((double) imu / nmu);
      for (int it(1); it <= nt; ++it) {
	double const tt(tmax * it / nt);
	sprite.push_back(Sprite(tt, computeWave(base, iseg, mu, tt)));
      }
    }
  }
}


int main(int argc, char ** argv)
{
  if (SIG_ERR == signal(SIGINT, handle)) {
    err (EXIT_FAILURE, "signal");
  }
  
  if (0 != atexit(cleanup)) {
    err(EXIT_FAILURE, "atexit");
  }
  
  heap = minheap_create(1000);
  if ( ! heap) {
    err(EXIT_FAILURE, "maxheap_create");
  }
  
  printf("# initializing sprite\n");
  fflush(stdout);
  init_sprite(1.25, 50, 5.0, 50);
  
  printf("# allocating grid\n");
  fflush(stdout);
  Posegrid grid(0.0, 1.0, nn,
		0.0, 1.0, nn,
		0.0, M_PI, nn);

  printf("\n");
  for (size_t itheta(0); itheta < nn; ++itheta) {
    printf("# splot 'file' matrix i %zu u (%g+$1*%g):(%g+$2*%g):3 w l lc 0 t '%g deg'\n",
	   itheta, grid.x0_, grid.dx_, grid.y0_, grid.dy_, grid.theta0_ + itheta * grid.dtheta_ * 180.0 / M_PI);
  }

  printf("# initializing lookup\n");
  fflush(stdout);
  lookup.resize(nn);
  for (size_t ix(0); ix < nn; ++ix) {
    lookup[ix].resize(nn);
    for (size_t iy(0); iy < nn; ++iy) {
      lookup[ix][iy].resize(nn);
      for (size_t itheta(0); itheta < nn; ++itheta) {
	lookup[ix][iy][itheta].pose_ = &grid.get(ix, iy, itheta);
      }
    }
  }
  Posegrid::index idx(grid.snap(sprite[0]));
  lookup[idx.ix][idx.iy][idx.itheta].tt_ = 0.0;
  lookup[idx.ix][idx.iy][idx.itheta].key_ = 0.0;
  
  printf("# propagating...\n");
  fflush(stdout);
  if (0 != heap_insert(heap, 0.0, &lookup[idx.ix][idx.iy][idx.itheta])) {
    err(EXIT_FAILURE, "heap_insert");
  }
  
  for (Entry * src((Entry*) heap_pop(heap)); (src != 0) && ( ! interrupt); src = (Entry*) heap_pop(heap)) {
    
    printf("# %zu: %g\n", HEAP_LENGTH(heap), src->tt_);
    fflush(stdout);
    
    src->key_ = -1.0;
    for (size_t is(1); (is < sprite.size()) && ( ! interrupt); ++is) {
      
      idx = grid.snap(*src->pose_ + sprite[is]);
      Pose const & candidate(grid.get(idx));
      if (&candidate == src->pose_) {
	continue;
      }
      
      double const tt(src->tt_ + sprite[is].tt_);
      Entry * dst(&lookup[idx.ix][idx.iy][idx.itheta]);
      if (dst->tt_ <= tt) {
	continue;
      }
      dst->tt_ = tt;
      
      if (dst->key_ >= 0.0) {
	if (0 != heap_change_key(heap, dst->key_, tt, dst)) {
          err(EXIT_FAILURE, "heap_change_key");
        }
      }
      else {
	if (0 != heap_insert(heap, tt, dst)) {
	  err(EXIT_FAILURE, "heap_insert");
	}
      }
      dst->key_ = tt;
    }
  }

  if (interrupt) {
    fprintf(stderr, "# ...interrupted\n");
  }
  
  for (size_t itheta(0); itheta < nn; ++itheta) {
    printf("\n\n# itheta: %zu\n", itheta);
    for (size_t ix(0); ix < nn; ++ix) {
      for (size_t iy(0); iy < nn; ++iy) {
	printf("%8g  ", lookup[ix][iy][itheta].tt_ == numeric_limits<double>::max() ? -1.0 : lookup[ix][iy][itheta].tt_);
      }
      printf("\n");
      fflush(stdout);
    }
  }
  
  printf("\n");
  for (size_t itheta(0); itheta < nn; ++itheta) {
    printf("# splot 'file' matrix i %zu u (%g+$1*%g):(%g+$2*%g):3 w l lc 0 t '%g deg'\n",
	   itheta, grid.x0_, grid.dx_, grid.y0_, grid.dy_, grid.theta0_ + itheta * grid.dtheta_ * 180.0 / M_PI);
  }
}
