/**
 * @file  mri_segment.c
 * @brief segments white matter and gray matter` from a brain volume
 * prior based segmentation
 * Arman Eshaghi 
 * contact: arman.eshaghi@me.com
 * For evaluation with Doug Greve, started in Sep2014
 * */
const char *MRI_SEGMENT_VERSION = "$Revision: Arman_evaluation";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "diag.h"
#include "error.h"
#include "macros.h"
#include "utils.h"
#include "proto.h"
#include "classify.h"
#include "mrisegment.h"
#include "mri.h"
#include "tags.h"
#include "mrinorm.h"
#include "timer.h"
#include "version.h"

const char *Progname ;
int main(int argc, char *argv[]) ;
int totalNumberOfClasses;

static int get_option(int argc, char *argv[]) ;



MRI *MRIsumPriorProbability(MRI *mri_prior_wm, *mri_prior_gm) ; 

static MRI *MRIrecoverBrightWhite(MRI *mri_T1, MRI *mri_src, MRI *mri_dst,
                                  float wm_low, float wm_hi, float slack,
                                  float pct_thresh) ;
static int is_diagonal(MRI *mri, int x, int y, int z) ;
/* Labels for prior labels (GM and WM)
 *
 */
int
main(int argc, char *argv[])
{ 
  MRI     *mri_src, *mri_dst, *mri_tmp, *mri_sum ;
  char    *input_file_name, *output_file_name, *gm_prior_probability_file_name, *wm_prior_probability_file_name ;
  int     nargs, i, msec ;
  struct timeb  then ;
  float   white_mean, white_sigma, gray_mean, gray_sigma ;
 
  char cmdline[CMD_LINE_LEN] ;
 
  TAGmakeCommandLineString(argc, argv, cmdline) ;

  /* initializing volumes from command line
   */
  TimerStart(&then) ;
  input_file_name = argv[1] ;
  output_file_name = argv[2] ;
  gm_prior_probability_file_name = argv[3] ;
  wm_prior_probability_file_name = argv[4] ;
  /* reading scans from files */
  mri_src = MRIread(input_file_name) ;
  prior_gm =  MRIread(gm_prior_probability_file_name) ;
  prior_wm =  MRIread(wm_prior_probability_file_name) ;

  if (!mri_src )
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
               Progname, input_file_name) ;
    MRIaddCommandLine(mri_src, cmdline) ;
  }
  else if (!prior_gm )
  {

    ErrorExit(ERROR_NOFILE, "%s: could not read prior GM volume from %s",
                Progname, input_file_name) ;
    MRIaddCommandLine(prior_gm, cmdline) ;
  }
  else if (!prior_wm )
  {
  
    ErrorExit(ERROR_NOFILE, "%s: could not read prior WM volume from %s",
                Progname, input_file_name) ;
    MRIaddCommandLine(prior_wm, cmdline) ;
  }

  if (mri_src->type != MRI_UCHAR)
  {
    MRI *mri_tmp ;
    printf("changing input type from %d to UCHAR\n", mri_src->type) ;
    mri_tmp = MRIchangeType(mri_src, MRI_UCHAR, 0, 1000, 1) ;
    MRIfree(&mri_src) ;
    mri_src = mri_tmp ;
  }

  if (prior_gm->type != MRI_UCHAR)
  {
    MRI *mri_tmp ;
    printf("changing input type from %d to UCHAR\n", prior_gm->type) ;
    mri_tmp = MRIchangeType(prior_gm, MRI_UCHAR, 0, 1000, 1) ;
    MRIfree(&prior_gm) ;
    prior_gm = mri_tmp ;
  }

  if (prior_wm->type != MRI_UCHAR)
  {
     MRI *mri_tmp ;
     printf("changing input type from %d to UCHAR\n", prior_wm->type) ;
     mri_tmp = MRIchangeType(prior_wm, MRI_UCHAR, 0, 1000, 1) ;
     MRIfree(&prior_wm) ;
     prior_wm= mri_tmp ;
  }
  /*Finite mixture model to calculate likelihood
   * we assume 2 clusters (K) for each probability map
   * 1. Caculating sum of prior probability at each voxel */
  MRIsumPriorProbability(prior_gm, prior_wm, mri_sum)
  MRIwrite(mri_dst, output_file_name) ;


  MRIfree(&mri_src) ;
  MRIfree(&prior_gm ) ;
  MRIfree(&prior_wm ) ;
  msec = TimerStop(&then) ;
  exit(0) ;
}

MRI *
MRIsumPriorProbability(MRI *mri_prior_gm, MRI *mri_prior_wm,  MRI *mri_sum)
{
  int  x, y, z, width, height, depth, nchanged, ntested ;
  float wmVal, gmVal,  sumVal ;
 
  fprintf(stderr, "Summing prior probability images\n") ;
 
  nchanged = ntested = 0 ;
  width = mri_prior_gm->width ;
  height = mri_prior_gm->height ;
  depth = mri_prior_gm->depth ;
 
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
       wmVal = MRIgetVoxVal(mri_prior_wm, x, y, z, 0) ; 
       gmVal = MRIgetVoxVal(mri_prior_gm, x, y, z, 0) ; 
       sumVal = gmVal + wmVal ;
       MRIsetVoxVal(mri_sum, x, y, z, 0, sumVal) ;
      }
    }
  }




























