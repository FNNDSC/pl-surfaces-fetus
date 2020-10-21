#!/usr/bin/perl
#
# Author: Jennings Zhang <jenni_zh@protonmail.com>
#         Claude Lepage <claude@bic.mni.mcgill.ca>
#
# May 2011
# Updated Aug 2019
#
# Copyright Alan C. Evans
# Professor of Neurology
# McGill University
#

# Extract a white matter surface using a marching-cubes
# algorithm, then run ASP to complete the convergence.

# removed from marching_cubes by jennings
# 16 May 2019
# removed a bunch of unused code
# - -calibrate
# - t1.mnc
# - skull_mask.mnc
# - -refine (if not set, asp stops when resolution is 326680)
# - cls.mnc
# surf_reg (surface registration)
# 20 May 2019
# removed $out_dist from schedule, -boundary not used
# 22 May 2019
# removed mask-label, not used

use strict;
use warnings "all";
use File::Basename;
use File::Spec;
use File::Temp qw/ tempdir /;
use File::Copy;

use Getopt::Tabular;
use MNI::Startup;
use MNI::FileUtilities;
use MNI::DataDir;

# --- set the help & usage strings ---
my $help = <<HELP;
Required parameters:
  wm_mask.mnc  : white matter mask
  white.obj    : white matter surface (output)
HELP

my $license = <<LICENSE;
Copyright Alan C. Evans
Professor of Neurology
McGill University

LICENSE

my $usage = <<USAGE;
Usage: $ProgramName [-left|-right] [-label N] wm_mask.mnc white.obj
       $ProgramName -help to list options

$license
USAGE

Getopt::Tabular::SetHelp( $help, $usage );

my $side = undef;
my $subsample = "";
my $xfm = undef;
my $label = 0;
my $save_chamfer = undef;
my $age = 20.0;
my @options = (
  ['-left', 'const', "Left", \$side, "Extract left surface"],
  ['-right', 'const', "Right", \$side, "Extract right surface"],
  ['-subsample', 'const', "-subsample", \$subsample,
   "Subsample white matter mask at half voxel size"],
  ['-xfm', 'string', 1, \$xfm,
   "Transformation from t1.mnc space to mni stereotaxic space"],
  ['-label', 'integer', 1, \$label,
   "Indicates that the given mask is actually a painted volume segmentation.\n"
   . "A binary mask is created for the supplied label value (integer).\n"
   . "\n1=CSF,"
   . "\n2=gray matter (corticol plate),"
   . "\n3=white matter (subplate zone),"
   . "\n4=intermediate zone,"
   . "\n5=subventricular zone,"
   . "\n6=ventricle."],
   ['-age', 'float', 1, \$age,
   "Prevent overfitting by increasing voxel size to match edge lengths."],
   # ['-sw', 'float', 1, \$sw,
   # "ASP stretch weight regulates edge length and causes mesh shrinkage."],
   # ['-lw', 'float', 1, \$lw,
   # "ASP laplacian weight sets how tightly to fit the surface."],
  );

GetOptions( \@options, \@ARGV ) or exit 1;
die "$usage\n" unless @ARGV == 2;
my $original_white_matter_mask = shift;
my $white_surface = shift;

if( !( -e $original_white_matter_mask ) ) {
  print "$usage\n";
  die "White matter mask must exist.\n";
}

if( !( defined $side ) ) {
  die "You must specify -left or -right hemisphere.\n";
}


my $tmpdir = &tempdir( "mcubes-XXXXXX", TMPDIR => 1, CLEANUP => 1 );

if ( $label > 0 ) {
  print "Creating a binary mask from the input segmentation.\n";
  $label -= 0.5;
  &run( 'minccalc', '-quiet', '-byte', '-clob', '-expr', "A[0]>$label",
        $original_white_matter_mask, "${tmpdir}/original_wm_mask.mnc");
  $original_white_matter_mask = "${tmpdir}/original_wm_mask.mnc";
}
undef $label;


my $ICBM_white_model = MNI::DataDir::dir("surface-extraction") .
                       "/white_model_320.obj";
my $initial_model = "${tmpdir}/initial_white_model.obj";
if( $side eq "Left" ) {
  &run( "param2xfm", "-translation", -25, 0, 0,
        "${tmpdir}/slide_left.xfm" );
  &run( "transform_objects", $ICBM_white_model,
        "${tmpdir}/slide_left.xfm", $initial_model );
  unlink( "${tmpdir}/slide_left.xfm" );
} elsif ( $side eq "Right" ) {
  &run( "param2xfm", "-scales", -1, 1, 1,
        "${tmpdir}/flip.xfm" );
  &run( "transform_objects", $ICBM_white_model,
        "${tmpdir}/flip.xfm", $initial_model );
  unlink( "${tmpdir}/flip.xfm" );
  &run( "param2xfm", "-translation", 25, 0, 0,
        "${tmpdir}/slide_right.xfm" );
  &run( "transform_objects", $initial_model,
        "${tmpdir}/slide_right.xfm", $initial_model );
  unlink( "${tmpdir}/slide_right.xfm" );
}


if( defined( $xfm ) && ( -e $xfm ) ) {
  &run( 'xfminvert', '-clobber', $xfm, "${tmpdir}/mcubes_inv.xfm" );
  &run( "transform_objects", $initial_model,
        "${tmpdir}/mcubes_inv.xfm", $initial_model );
  unlink( "${tmpdir}/mcubes_inv.xfm" );
}

# Make sure that the mask is not empty.
my $sum = `mincstats -quiet -sum $original_white_matter_mask`;
chomp( $sum );
die "Empty white matter mask.\n" if( $sum == 0 );


# Remove loosely connected voxels and fill-in tightly connected voxels
# to have a smoother surface.

open KRNL, "> ${tmpdir}/ngh_count.kernel";
print KRNL "MNI Morphology Kernel File\n";
print KRNL "Kernel_Type = Normal_Kernel;\n";
print KRNL "Kernel =\n";
print KRNL " 1.0  0.0  0.0  0.0  0.0  1.0\n";
print KRNL "-1.0  0.0  0.0  0.0  0.0  1.0\n";
print KRNL " 0.0  1.0  0.0  0.0  0.0  1.0\n";
print KRNL " 0.0 -1.0  0.0  0.0  0.0  1.0\n";
print KRNL " 0.0  0.0  1.0  0.0  0.0  1.0\n";
print KRNL " 0.0  0.0 -1.0  0.0  0.0  1.0;\n";
close KRNL;

my $ngh_count = "${tmpdir}/ngh_count.mnc";
my $white_matter_mask = "${tmpdir}/wm_mask_clean.mnc";
my $tmp_wm_mask = "${tmpdir}/tmp_wm_mask.mnc";
my $wm_mask_defragged = "${tmpdir}/wm_mask_defragged.mnc";

&run( 'minccalc', '-quiet', '-clobber', '-unsigned', '-byte',
      '-expression', 'abs(A[0]-1)<0.5',
      $original_white_matter_mask, $wm_mask_defragged );
&run( 'mincdefrag', $wm_mask_defragged, $wm_mask_defragged, 1, 6 );
# might need to defrag 0 too, to remove holes inside mask

# We will use a white matter mask with simplified connectivity
# for the raw marching-cubes surface, but we will later fit the
# surface to the original white matter mask.
copy($wm_mask_defragged, $white_matter_mask);

my $lower = 2.5;
# May 29, 2019: changed 5 iterations down to 4
# Aug 23, 2019: this made no sense before
for( my $iter = 1; $iter <= 5; $iter++ ) {
  my $ret = `mincstats -sum $white_matter_mask`; chomp( $ret );
  print "Connectivity iteration $iter: $ret\n";
  &run( 'mincmorph', '-clobber', '-unsigned', '-byte', '-convolve', '-kernel',
        "${tmpdir}/ngh_count.kernel", $white_matter_mask, $ngh_count );
  &run( 'minccalc', '-clob', '-quiet', '-expr',
        "if(A[1]<$lower){0}else{if(A[1]>4.5){1}else{A[0]}}",
        $white_matter_mask, $ngh_count, $tmp_wm_mask );
  &run( 'mv', '-f', $tmp_wm_mask, $white_matter_mask );
}

unlink( $ngh_count );
unlink( "${tmpdir}/ngh_count.kernel" );

# Preparation of the white matter mask for extraction of the surface
# using marching cubes.

&run( 'minccalc', '-clobber', '-quiet', '-expression', 'out=1',
      '-unsigned', '-byte', $white_matter_mask, "${tmpdir}/filled.mnc" );
&run( 'surface_mask2', '-binary', "${tmpdir}/filled.mnc",
      $initial_model, "${tmpdir}/s0.mnc" );
&run( 'mincresample', '-clobber', '-quiet', '-like',
      "${tmpdir}/filled.mnc", "${tmpdir}/s0.mnc",
      "${tmpdir}/s1.mnc" );
&run( 'minccalc', '-clobber', '-quiet', '-expression',
      'if(A[0]>0.5||A[1]>0.5){1}else{0}',
      $white_matter_mask, "${tmpdir}/s1.mnc", "${tmpdir}/s0.mnc" );
&run( 'dilate_volume', "${tmpdir}/s0.mnc", "${tmpdir}/s1.mnc", 1, 26, 1 );
&run( 'minccalc', '-clobber', '-quiet', '-expression', 'A[0]+A[1]',
      $white_matter_mask, "${tmpdir}/s1.mnc", "${tmpdir}/s0.mnc" );
&run( 'mincreshape', '-quiet', '-clobber', '-unsigned', '-byte',
      '-image_range', 0, 255, '-valid_range', 0, 255, "${tmpdir}/s0.mnc",
      "${tmpdir}/s1.mnc" );
&run( 'mincdefrag', "${tmpdir}/s1.mnc", "${tmpdir}/s0.mnc", 2, 19 );  ## was 6 before
unlink( "${tmpdir}/filled.mnc" );
unlink( "${tmpdir}/s1.mnc" );

# Do not use sub-sampling if voxel resolution is already below 1mm.
# Too slow.

my $dx = `mincinfo -attvalue xspace:step $white_matter_mask`; chomp($dx);
$subsample = "" if( $dx < 0.90 );

# Crop volume to smallest extents to speed up marching-cubes.
my @crop = split( ' ', `mincbbox ${tmpdir}/s0.mnc -mincreshape` );
&run( 'mincreshape', '-quiet', '-clobber', @crop, "${tmpdir}/s0.mnc",
      "${tmpdir}/s1.mnc" );
unlink( "${tmpdir}/s0.mnc" );

# Run the marching-cubes algorithm on the mask.
#copy( "${tmpdir}/s1.mnc", "./s1.mnc"); #debug
& run ('sphere_mesh', "${tmpdir}/s1.mnc", $white_surface, $subsample);
unlink( "${tmpdir}/s1.mnc" );

# Transform the current raw surface to stereotaxic space for
# resampling and registration.

if( defined( $xfm ) && ( -e $xfm ) ) {
  &run( "transform_objects", $white_surface, $xfm, $white_surface );
}

# Coarsen and smooth the original marching-cubes surface.
# Wed, 22 May 2019 I don't think we need smoothing here.
# sphere_mesh already produces a smooth surface.
my $white_surface_sm = "${tmpdir}/white_surf_sm.obj";
# &run( 'adapt_object_mesh', $white_surface, $white_surface_sm,
#       120000, 10, 50, 1 );
copy( $white_surface, $white_surface_sm );

# Inflate the white surface onto a unit sphere.
my $white_sphere_sm = "${tmpdir}/white_sphere_sm.obj";
&run( 'inflate_to_sphere', $white_surface_sm, $white_sphere_sm );

# Interpolate from sphere-to-sphere to resample the white surface
# using the 40962 vertices on the standard ICBM surface average
# template. This unit sphere is the one used for surface registration.

my $unit_sphere = "${tmpdir}/unit_sphere.obj";
&run( 'create_tetra', $unit_sphere, 0, 0, 0, 1, 1, 1, 81920 );
if( $side eq "Right" ) {
  &run( "param2xfm", "-scales", -1, 1, 1,
        "${tmpdir}/flip.xfm" );
  &run( "transform_objects", $unit_sphere,
        "${tmpdir}/flip.xfm", $unit_sphere );
  unlink( "${tmpdir}/flip.xfm" );
}

# Evaluate the white surface from the marching-cubes surface.

&resample_white_surface( $white_surface_sm, $white_sphere_sm,
                         $unit_sphere, $white_surface, $xfm );

# Do surface registration on to the average white population model.
# This will be a first alignment that will ensure better isotropic mesh.
# SURFACE REGISTRAION REMOVED

unlink( $unit_sphere );

# Check for self-intersections in the marching-cubes surface.
my $prev_num_inter = 999999;
my $increasing = 0;
for( my $i = 1; $i <= 500; $i++ ) {
  my @ret = `check_self_intersect $white_surface -fix $white_surface`;
  $ret[1] =~ /Number of self-intersecting triangles = (\d+)/;
  my $num_inter = $1;
  printf "Iter = $i  Self-Intersections = $num_inter\n";
  if( $num_inter == 0 ) {
    $prev_num_inter = 0;
    last;
  }
  if( $num_inter >= $prev_num_inter ) {
    if( $increasing >= 0 ) {
      $increasing++;
    } else {
      $increasing = 1;
    }
  } else {
    $increasing = 0;
  }
  $prev_num_inter = $num_inter;
  last if( $increasing >= 5 && $i > 50 );
}
if( $prev_num_inter > 0 ) {
  my $failed_surface = $white_surface;
  $failed_surface =~ s/\.obj$/-failed\.obj/;
  `mv -f $white_surface $failed_surface`;
  # unlink( $white_surface );   ## is this desirable???
  die "Failed interpolation of marching-cubes surface with $prev_num_inter self-intersections.\n";
}

# Now we can transform back the resampled surface to native space
# to do the fitting against the volume.

if( defined( $xfm ) && ( -e $xfm ) ) {
  &run( 'xfminvert', '-clobber', $xfm, "${tmpdir}/mcubes_inv.xfm" );
  &run( "transform_objects", $white_surface,
        "${tmpdir}/mcubes_inv.xfm", $white_surface );
  unlink( "${tmpdir}/mcubes_inv.xfm" );
}


# Run ASP on the resampled white surface to converge it fully
# to the white matter mask. The ICBM model is used by ASP to
# define the distribution of the edge lengths for the stretch
# constraint. Here, we use the original white matter mask
# without the changes for dangling white voxels.

&run_asp( $white_surface, $wm_mask_defragged, $initial_model );


# Resample the white surface from a standard sphere from the
# hi-res marching-cubes white surface. The distribution of
# vertices on the sphere is adapted such as to produce an
# interpolated surface with triangles of nearly the same size.

sub resample_white_surface {

  my $white_mc = shift;     # hi-res raw marching-cubes surface
  my $sphere_mc = shift;    # sphere corresponding to white_mc
  my $unit_sphere = shift;  # standard sphere
  my $output = shift;       # output white surface with uniform triangles
  my $xfm = shift;          # transform from native to stx space

  my @conf = ( { 'size' => 320,     # most of the motion occurs early
                 'fwhm' => 20.0,
                 'niter' => 500 },
               { 'size' => 1280,
                 'fwhm' => 10.0,
                 'niter' => 500 },
               { 'size' => 5120,
                 'fwhm' => 5.0,
                 'niter' => 300 },
               { 'size' => 20480,
                 'fwhm' => 2.0,
                 'niter' => 150 } );

  my $start = 320;
  my $end = 20480;

  my $npolys = `print_n_polygons $unit_sphere`;
  chomp( $npolys );

  my $current_sphere = "${tmpdir}/current_sphere.obj";
  &run( 'cp', '-f', $unit_sphere, $current_sphere );

  # obtain initial white surface

  &run( 'interpolate_sphere', $white_mc, $sphere_mc,
        $current_sphere, $output );

  # Multi-resolution approach from coarse to fine mesh.

  for( my $idx = 0; $idx <= $#conf; $idx++ ) {

    my $size = $conf[$idx]{size};

    next if( $size < $start );
    last if( $size > $end );

    print "Sphere adaptation at $size vertices...\n";

    # Obtain the triangle areas from current white surface to
    # the current sphere at size npolys.

    my $white_area = "${tmpdir}/white_area.txt";
    my $sphere_area = "${tmpdir}/sphere_area.txt";
    &run( 'depth_potential', '-area_simple', $output, $white_area );
    &run( 'depth_potential', '-area_simple', $current_sphere, $sphere_area );
    &run( 'vertstats_math', '-old_style_file', '-div', $white_area,
          $sphere_area, $white_area );
    unlink( $sphere_area );
    if( $conf[$idx]{fwhm} > 0 ) {
      &run( 'depth_potential', '-smooth', $conf[$idx]{fwhm},
            $white_area, $output, $white_area );
    }

    # adapt the current_sphere at this size based on the areas.

    &subdivide_mesh( $current_sphere, $size, $current_sphere );
    &run( 'adapt_metric', $current_sphere, $white_area,
          $current_sphere, $conf[$idx]{niter} );
    unlink( $white_area );

    # interpolate relative to the original white surface at npolys.

    &subdivide_mesh( $current_sphere, $npolys, $current_sphere );
    &run( 'interpolate_sphere', $white_mc, $sphere_mc,
          $current_sphere, $output );

  }

  # Create a new hi-res background mesh with uniform triangles.
  # This new background mesh will be used for interpolating the
  # resampled white surface after surface registration. This way,
  # we can associate the standard sphere to this new bg mesh
  # since the standard sphere is used during surface registration.

  $npolys *= 4;
  &subdivide_mesh( $current_sphere, $npolys, $current_sphere );
  &run( 'interpolate_sphere', $white_mc, $sphere_mc,
        $current_sphere, $white_mc );
  &subdivide_mesh( $unit_sphere, $npolys, $sphere_mc );
  unlink( $current_sphere );

}

# Run ASP on the resampled white surface to converge it fully
# to the white matter mask. This is a simplified version of
# extract_white_surface without the coarse steps.
sub run_asp {

  my $surface = shift;
  my $wm_mask = shift;
  my $white_model = shift;

  my $self_dist2 = 0.001;
  my $self_weight2 = 1e08;
  my $n_selfs = 9;

  my $stop_threshold = 1e-3;
  my $stop_iters = 1000;
  my $n_per = 5;
  my $tolerance = 1.0e-03;
  my $f_tolerance = 1.0e-06;
  my $oo_scale = 0.5;

  # size  number of triangles
  # sw    weight for mesh stretching (bigger means less stretch)
  #       small value will create a highly voxelated/bumpy surface,
  #       but with good detail.
  #       large value for sw creates a very smooth surface that
  #       destroys topological detail.
  #       sw>40 will not fit into concavities that are 2 voxels wide.
  # n_itr number of iterations
  # inc   save every 20 iters
  # l_w   weight for Laplacian field: large value means tighter fit,
  #       but it seems that a small value is better convergence wise
  #       to ensure mesh quality.
  # iso   target value of the LaPlacian field to converge to
  # si    max step increment (large is faster, but risk of intersection)
  # l_s   do subsampling (0=none, n=#points extra along edge);
  #       (helps to get through isolated csf voxels in insula region,
  #       but must finish off the job without subsampling)
  # iw    weight for surface self-intersection
  # self  min distance to check for surface self-intersection
  #       (0.0625 found to be too high)

  my @schedule = (
  #  size    sw  n_itr  inc    l_w  iso    si   l_s    iw    self
  # -----   ---  -----  ---   ----  ---  ----  ----  ----   -----
    81920,   20,  400,   50,  8e-6,  10, 0.10,  0.0,  1e0,   0.01,
  );

  # parameters to fit the surface from sphere interpolation is simpler
  # because the surface is already close to its target.
  # sw maintains good vertex distribution, but shrinkage causes
  # surface to escape from narrow sulci so don't set sw too high

  my $chamfer_range = 5;
  my $slope = 1;
  my $scale_xfm = 0;

  if ( $age < 30.0 ) {
    print "Increasing mask volume.\n";
    # $scale should depend on the ratio of average edge length to voxel side
    # `surface-stats -edge_length $surface`
    my $scale = 3;
    $scale_xfm = "${tmpdir}/make_bigger.xfm";
    &run( "param2xfm", "-scale", $scale, $scale, $scale, $scale_xfm );
    &run( "transform_volume", $wm_mask, $scale_xfm, $wm_mask);
    &run( "transform_objects", $surface, $scale_xfm, $surface);
    $chamfer_range = $scale * 5;
    $slope = 1 / $scale;
  }

  my $chamfer_map = "${tmpdir}/simple_chamfer.mnc";
  simple_chamfer( $wm_mask, $chamfer_map, $tmpdir, $chamfer_range, $slope );
  copy( $chamfer_map, $save_chamfer ) if ( defined( $save_chamfer) );

  copy( $white_model, "${tmpdir}/white_model_tmp.obj" );
  $white_model = "${tmpdir}/white_model_tmp.obj";
  subdivide_mesh( $white_model, $schedule[0], $white_model );

  # Do the fitting stages like gray surface expansion.
  my $sched_size = 10;
  my $num_steps = @schedule / $sched_size;
  my $num_rows = @schedule / $sched_size;

  for( my $i = 0;  $i < @schedule;  $i += $sched_size ) {
    my $row = $i / $sched_size + 1;

    my ( $size, $sw, $n_iters, $iter_inc, $laplacian_weight, $iso,
         $step_increment, $oversample, $self_weight, $self_dist,
         ) = @schedule[$i..$i+$sched_size-1];

    $oversample *= $oo_scale;
    my $self2 = get_self_intersect( $self_weight, $self_weight2, $n_selfs,
                                    $self_dist, $self_dist2 );

    for( my $iter = 0;  $iter < $n_iters;  $iter += $iter_inc ) {
      print "echo Step ${size}: $iter / $n_iters    sw=$sw,  "
            . "Schedule row ${row} / ${num_rows}\n";

      my $ni = $n_iters - $iter;
      $ni = $iter_inc if( $ni > $iter_inc );

      # don't use taubin smoothing

      my $command = "surface_fit " .
                    "-mode two -surface  $surface $surface" .
                    # "-taubin_smoothing" .
                    # white_model is rightfully ignored in surface_fit
                    " -stretch $sw $white_model -.9 0 0 0" .
                    " -laplacian $chamfer_map $laplacian_weight 0 " .
                    " $iso $oversample " .
                    " $self2 -step $step_increment " .
                    " -fitting $ni $n_per $tolerance " .
                    " -ftol $f_tolerance " .
                    " -stop $stop_threshold $stop_iters ";
      print $command . "\n";
      system( $command ) == 0 or die "Command $command failed with status: $?";
    }
  }
  unlink( $white_model );
  if ( $scale_xfm ) {
    &run( "xfminvert", "-clobber", $scale_xfm, $scale_xfm);
    &run( "transform_objects", $surface, $scale_xfm, $surface);
  }
}

# Check if the input surface has the same side orientation (left)
# as the default template model.

sub CheckFlipOrientation {

  my $obj = shift;

  my $npoly = `print_n_polygons $obj`;
  chomp( $npoly );

  my $ret = `tail -5 $obj`;
  my @verts = split( ' ', $ret );
  my @last3 = ( $verts[$#verts-2], $verts[$#verts-1], $verts[$#verts] );

  my $dummy_sphere = "${tmpdir}/dummy_sphere.obj";
  &run('create_tetra',$dummy_sphere,0,0,0,1,1,1,$npoly);
  $ret = `tail -5 $dummy_sphere`;
  unlink( $dummy_sphere );
  @verts = split( ' ', $ret );
  my @sphere3 = ( $verts[$#verts-2], $verts[$#verts-1], $verts[$#verts] );
  if( $last3[0] == $verts[$#verts-2] &&
      $last3[1] == $verts[$#verts-1] &&
      $last3[2] == $verts[$#verts-0] ) {
    return 0;
  } else {
    return 1;
  }
}

# subdivide a surface taking into account if it's a left or right hemisphere.

sub subdivide_mesh {

  my $input = shift;
  my $npoly = shift;
  my $output = shift;

  my $npoly_input = `print_n_polygons $input`;
  chomp( $npoly_input );
  return if ($npoly_input == $npoly); # do nothing if sizes are same

  if( !CheckFlipOrientation( $input ) ) {
    &run( "subdivide_polygons", $input, $output, $npoly );
  } else {
    # flip right as left first before subdividing, then flip back.
    &run( "param2xfm", '-clobber', '-scales', -1, 1, 1,
          "${tmpdir}/flip.xfm" );
    my $input_flipped = "${tmpdir}/right_flipped.obj";
    &run( "transform_objects", $input,
          "${tmpdir}/flip.xfm", $input_flipped );
    &run( "subdivide_polygons", $input_flipped, $output, $npoly );
    &run( "transform_objects", $output,
          "${tmpdir}/flip.xfm", $output );  # flip.xfm is its own inverse
    unlink( $input_flipped );
    unlink( "${tmpdir}/flip.xfm" );
  }
}

sub simple_chamfer {

  my $wm_mask = shift;
  my $output_chamfer = shift;
  my $tmpdir = shift;
  my $dist = shift;
  my $slope = shift;
  # expect white matter mask to have already been defragmented
  &run( 'mincchamfer', '-quiet', '-max_dist', $dist,
        $wm_mask, "${tmpdir}/chamfer_outside.mnc" );
  &run( 'minccalc', '-quiet', '-clobber', '-expression', '1-A[0]',
        $wm_mask, "${tmpdir}/wm_mask_negated.mnc" );
  &run( 'mincchamfer', '-quiet', '-max_dist', $dist,
        "${tmpdir}/wm_mask_negated.mnc", "${tmpdir}/chamfer_inside.mnc" );
  unlink( "${tmpdir}/wm_mask_negated.mnc" );
  &run( 'minccalc', '-quiet', '-clob', '-expression', "10.0+(A[1]-A[0])*$slope",
        "${tmpdir}/chamfer_inside.mnc", "${tmpdir}/chamfer_outside.mnc",
        $output_chamfer );
  unlink( "${tmpdir}/chamfer_outside.mnc" );
  unlink( "${tmpdir}/chamfer_inside.mnc" );
  return $output_chamfer;
}


# copied from lib/surface-extraction/deform_utils.pl

sub  get_self_intersect( $$$$$ )
{
  my( $self_weight, $self_weight2, $n_selfs, $self_dist, $self_dist2 ) = @_;
  my( $self, $weight, $weight_factor, $s, $dist );

  if( $self_weight > 0.0 ) {
    $self = "";
    $weight = $self_weight;
    $weight_factor = ( $self_weight2 / $self_weight ) **
                      ( 1.0 / $n_selfs );

    for( $s = 0;  $s < $n_selfs;  ++$s ) {
      $dist = $self_dist + ($self_dist2 - $self_dist) *
              $s / $n_selfs;
      $self = $self . " -self_intersect $weight $dist ";
      $weight *= $weight_factor;
    }
    $self = $self . " -self_intersect $self_weight2 $self_dist2 ";
  } else {
    $self = "";
  }
  $self;
}

# Execute a system call.

sub run {
  # print "@_\n";
  system(@_)==0 or die "Command @_ failed with status: $?";
}
