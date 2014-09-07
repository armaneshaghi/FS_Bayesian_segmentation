/**
 * @file  mri_segment.c
 * @brief segments white matter from a brain volume
 * evaluation for prior based segmentation
 * Arman Eshaghi 
 * contact: arman.eshaghi@me.com
 * Evaluation with Doug, Sep2014
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
MRI *MRIremoveWrongDirection(MRI *mri_src, MRI *mri_dst, int wsize,
                             float low_thresh, float hi_thresh,
                             MRI *mri_labels) ;
MRI *MRIfindBrightNonWM(MRI *mri_T1, MRI *mri_wm) ;
MRI *MRIfilterMorphology(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfillBasalGanglia(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIfillVentricles(MRI *mri_src, MRI *mri_dst) ;
MRI *MRIremove1dStructures(MRI *mri_src, MRI *mri_dst, int max_iter,
                           int thresh, MRI *mri_labels) ;

static MRI *MRIrecoverBrightWhite(MRI *mri_T1, MRI *mri_src, MRI *mri_dst,
                                  float wm_low, float wm_hi, float slack,
                                  float pct_thresh) ;
static int is_diagonal(MRI *mri, int x, int y, int z) ;
/* For now I will define only two labels for prior labels (GM and WM)
 *
 */
int
main(int argc, char *argv[])
{  MRI     *mri_src, *mri_dst, *mri_tmp, *mri_labeled, *mri_labels;
   char    *input_file_name, *output_file_name ;
   int     nargs, i, msec ;
   struct timeb  then ;
   float   white_mean, white_sigma, gray_mean, gray_sigma ;
 
   char cmdline[CMD_LINE_LEN] ;
 
   TAGmakeCommandLineString(argc, argv, cmdline) ;

   /* initializing prior probability images
    */
   TimerStart(&then) ;
   input_file_name = argv[1] ;
   output_file_name = argv[2] ;
   gm_prior_probability_file_name = argv[3] ;
   wm_prior_probability_file_name = argv[4] ;

   mri_src = MRIread(input_file_name) ;
   mri_gm = 
   mri_wm = 
   if (!mri_src)
     ErrorExit(ERROR_NOFILE, "%s: could not read source volume from %s",
               Progname, input_file_name) ;
   MRIaddCommandLine(mri_src, cmdline) ;
   if (mri_src->type != MRI_UCHAR)
   {
     MRI *mri_tmp ;
     printf("changing input type from %d to UCHAR\n", mri_src->type) ;
     mri_tmp = MRIchangeType(mri_src, MRI_UCHAR, 0, 1000, 1) ;
     MRIfree(&mri_src) ;
     mri_src = mri_tmp ;
   }
