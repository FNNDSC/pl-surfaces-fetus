# Deprecated: Please see https://github.com/FNNDSC/Fetal_Brain_MRI_Surface_Extraction_Pipeline

# pl-surfaces-fetus

## Table of Contents

* [Abstract](#abstract)
* [Description](#description)
* [Publications](#publications)
* [Usage](#usage)
    * [Required Arguments](#required-arguments)
    * [Optional Output Options](#optional-output-options)
    * [Output](#output)
        * [Files](#files)
        * [Visualization](#visualization)
* [Example](#example)
* [Build](#build)
* [TODO](#todo)
    * [Comments](#comments)
* [License](#license)

## Abstract

Extract surfaces (`.obj`) for the *subplate zone*
from pre-segmented fetal brain MRI volume (`.mnc`) using
**CIVET**'s **marching-cubes**/`sphere_mesh` and **ASP**/`surface_fit`.

![brain-view white-matter surface mesh](examples/spectral.png)

## Description

- [OHBM 2020 Poster](https://fnndsc.childrens.harvard.edu/conferences/2020/OHBM/Jennings/Jennings_Zhang_OHBM_2020_Subplace_Surfaces.pdf)
- [Another One](https://fnndsc.childrens.harvard.edu/conferences/2020/OHBM/Jennings/Lana_Vasung_OHBM_2020_Thickness_GE.pdf)


Given an input of a segmented volume (`*.mnc`) for fetal brain MRI between 20-33
weeks of gestational age (GA), this pipeline will produce surface meshes for
vertex-wise quantitative analysis of the subplate zone.

First, a marching-cubes[1] algorithm is used (`sphere_mesh`) to extract the
gray-white boundary as a surface mesh (`wm_81920.obj`). Next, mesh deformation[2]
by the ASP algorithm (`surface_fit`) using a radial distance map shrinks the
surface inwards to find the intermediate zone (IZ)'s superficial boundary
(`iz_81920.obj`), producing two surfaces over the subplate zone (SP) with one-to-one
vertex correspondence. The schedule of parameters used to run `surface_fit` depend
on age. In reality, these parameters should depend on brain size and gyrification index,
which are both correlated with age in normal neurodevelopment of a healthy fetus.
Finally, subplate thickness is calculated as the Euclidean distance between
corresponding vertices using the "tlink"[3] method.

    cortical_thickness -tlink iz_81920.obj wm_81920.obj subplate_thickness.txt

[1] Lepage C. (2017), "Human MR Evaluation of Cortical Thickness Using CIVET v2.1", *OHBM*,
http://mcin.ca/wp-content/uploads/2017/08/HBM2017_civet.png

[2] Kim et. al, (2005). "Automated 3-D extraction and evaluation of the inner and outer
cortical surfaces using a Laplacian map and partial volume effect classification".
*NeuroImage* 27, pp. 210-221.

[3] Lerch J (2005), "In-vivo Analysis of Cortical Thickness using Magnetic Resonance Images",
http://www.bic.mni.mcgill.ca/~jason/jpl-thesis-submitted.pdf

### Publications

Lana Vasung, Caitlin K Rollins, Hyuk Jin Yun, Clemente Velasco-Annis, Jennings Zhang,
Konrad Wagstyl, Alan Evans, Simon K Warfield, Henry A Feldman, P Ellen Grant, Ali Gholipour (2019).
"Quantitative In vivo MRI Assessment of Structural Asymmetries and Sexual Dimorphism of Transient
Fetal Compartments in the Human Brain." *Cerebral Cortex*. https://doi.org/10.1093/cercor/bhz200

Lana Vasung, Caitlin K Rollins, Clemente Velasco-Annis, Hyuk Jin Yun, Jennings Zhang,
Simon K Warfield, Henry A Feldman, Ali Gholipour, P Ellen Grant (2019).
"Spatiotemporal Differences in the Regional Cortical Plate and Subplate Volume Growth
during Fetal Development." *Cerebral Cortex*. https://doi.org/10.1093/cercor/bhaa033

Lana Vasung, Chenying Zhao, Matthew Barkovich, Caitlin K Rollins, Jennings Zhang, Claude Lepage, Teddy Corcoran, Clemente Velasco-Annis, Hyuk Jin Yun, Kiho Im, Simon Keith Warfield, Alan Charles Evans, Hao Huang, Ali Gholipour, Patricia Ellen Grant, Association between Quantitative MR Markers of Cortical Evolving Organization and Gene Expression during Human Prenatal Brain Development, Cerebral Cortex, Volume 31, Issue 8, August 2021, Pages 3610–3621, https://doi.org/10.1093/cercor/bhab035

## Usage

```bash
docker run -u $(id -u) --rm \
    -v <INPUTDIR>:/incoming -v <OUTPUTDIR>:/outgoing \
    surfaces_fetus --side <LEFT|RIGHT> --age <N> /incoming /outgoing
```

### Required Arguments

`<INPUTDIR>` must contain a single file `*.mnc`
of segmentations for a **single brain hemisphere** (either left or right)
with the following labels:


Label | Region
------|-------
3     | gray-white matter boundary
4     | superficial boundary of intermediate zone (IZ)


Also, an estimate of the subject's GA must be given.
The subplate zone changes drastically between 30-32 GA,
hence a different algorithm is used to fit the inner surface mesh.

    [--age <N>]
    Estimate of subject's gestational age (GA), in weeks.

    [--side <left|right>]
    Specify brain hemisphere


### Optional Output Options

Generate additional output which can be helpful for debugging or visualization.

    [--keep-intermediate]
    Keep intermediate files (e.g. *mask.mnc, *chanfer.mnc)

    [--qc]
    Save surface_fit logs and produce vertex-wise quality check files,
    such as fitting distance error, curvature, and triangle area.

### Output

#### Files

Output File Name                        | Purpose
----------------------------------------|--------
`wm_81920.obj`                          | subplate outer surface
`iz_81920.obj`                          | subplate inner surface
`intermediates/wm_mask.mnc`             | subplate outer mask
`intermediates/iz_mask.mnc`             | subplate inner mask
`intermediates/iz_chamfer.mnc`          | distance map to inner surface
`qc/wm_cubes.log`                       | surface extraction log, preprocessing and `surface_fit`
`qc/iz_fit.log`                         | fitting `surface_fit` log
`qc/wm_dist.txt`                        | marching-cubes distance error
`qc/wm_smth.txt`                        | marching-cubes smoothness quality
`qc/wm_area.txt`                        | marching-cubes triangle areas quality
`qc/iz_dist.txt`                        | fitting distance error
`qc/iz_smth.txt`                        | fitting smoothness quality
`qc/iz_area.txt`                        | fitting triangle areas quality
`sp_thickness_tlink.txt`                | `-tlink` thickness between surfaces
`qc/sp_thickness_tnear.txt`             | `-tnear` thickness between surfaces
`qc/tlink_minus_tnear.txt`              | vertex-wise difference between `-tlink` and `-tnear` thicknesses
`diemask.txt`                           | vertex mask for diencephalon
`not_subplate_mask.txt`                 | vertex mask over region where the subplate is discontinuous
`intermediates/mid_81920.obj`           | midpoints between inner and outer surface
`qc/distortion_angles.txt`              | distortion angles between 0-pi

#### Visualization

[MNI Display](http://www.bic.mni.mcgill.ca/software/Display/Display.html)
should be used to view `*.mnc` files. It supports overlay of `*.obj` files.

    Display -spectral labels.mnc wm_81920.obj iz_81920.obj

`brain-view` is used to view `*.obj` surface meshes in 3D, along with overlaying
data or a mask `*.txt`.

    brain-view wm_81920.obj subplate_thickness.txt

Warning: `brain-view` is legacy software, I have spent 30 hours trying to compile it
on a modern GNU/Linux and have made no progress.

## Example

TODO sample data doesn't exist yet, because our lab is stingy about data sharing :(

I found some public stuff here:

- [10.18112/openneuro.ds003105.v1.0.1](https://openneuro.org/datasets/ds003105/versions/1.0.1/) - CC0 license
    - https://arxiv.org/pdf/2009.06275.pdf
- http://neuroimaging.ch/feta

```bash
docker run -u $(id -u) --rm \
    -v $PWD/examples/brain:/incoming:ro $PWD/out:/outgoing:rw \
    fnndsc/pl-surfaces-fetus surfaces_fetus \
    --age 31.29 --side left  \
    --qc --keep-intermediate \
    /incoming /outgoing
```

## Build

Most simply,

```bash
docker build -t pl-surfaces-fetus .
```

Or cross-compile a multi-platform manifest list, and push to Dockerhub

```bash
docker run --rm --privileged aptman/qus -s -- -p ppc64le
docker buildx create --name moc_builder --use
docker buildx build -t fnndsc/pl-surfaces-fetus --platform linux/amd64,linux/ppc64le .
# optional clean up
docker buildx rm
docker run --rm --privileged aptman/qus -- -r
```

## TODO

- `--verbose` option to display helpful messages when running locally
- expansion to *cortical plate*
- `--clobber` option
- forward arguments to `cortical_thickness`
- reuse chamfer for both `surface_fit` and `volume_object_evaluate`

#### Comments

- `parm.py` mathematical model, we're workin' on it
- `edgy.py` and `triangle_aspect.py` not currently used, since they are confounding with `depth_potential -area_simple`

## License

[Civet core](https://github.com/aces/CIVET_Full_Project/blob/master/LICENSE)

- open source
- no commercial use
