import subprocess as sp
from os import mkdir, path
from tempfile import gettempdir
from glob import glob


class UserError(Exception):
    """
    User gave a bad input.
    """
    pass


def get_input_file(in_dir: str, g='*.mnc') -> str:
    loc = path.join(in_dir, g)
    files = glob(loc)
    if len(files) == 0:
        raise UserError(f'no input file - cannot find {loc}')
    elif len(files) > 1:
        raise UserError(f'more than one {g} present in {in_dir}')
    return files[0]


def process(in_dir: str, out_dir: str, side: str, age: float, keep_intermediate: bool, qc: bool):
    segmentation_mnc = get_input_file(in_dir)
    age = str(age)  # will get passed to subprocess.run
    side = side.lower()
    if side not in ('left', 'right'):
        raise ValueError('"--side" must be either "left" or "right"')

    # we don't care about cleanup nor using tempfile's more advanced methods
    # because we assume this script is being run in a stateless container
    intf = gettempdir()
    qcf = gettempdir()

    if keep_intermediate:
        intf = path.join(out_dir, 'intermediate')
        mkdir(intf)

    if qc:
        qcf = path.join(out_dir, 'qc')
        mkdir(qcf)

    layer3_obj = path.join(out_dir, 'wm_81920.obj')
    layer3_log = path.join(qcf, 'wm_cubes.log')
    layer3_mask_mnc = path.join(intf, 'wm_mask.mnc')
    layer3_chamfer_mnc = path.join(intf, 'wm_chamfer.mnc')
    layer3_dist_txt = path.join(qcf, 'wm_dist.txt')
    layer3_smth_txt = path.join(qcf, 'wm_smth.txt')
    layer3_area_txt = path.join(qcf, 'wm_area.txt')
    layer4_obj = path.join(out_dir, 'iz_81920.obj')
    layer4_log = path.join(qcf, 'iz_fit.log')
    layer4_mask_mnc = path.join(intf, 'iz_mask.mnc')
    layer4_chamfer_mnc = path.join(intf, 'iz_chamfer.mnc')
    layer4_dist_txt = path.join(qcf, 'iz_dist.txt')
    layer4_smth_txt = path.join(qcf, 'iz_smth.txt')
    layer4_area_txt = path.join(qcf, 'iz_area.txt')
    thickness_tlink = path.join(out_dir, 'sp_thickness_tlink.txt')
    thickness_tnear = path.join(qcf, 'sp_thickness_tnear.txt')
    tlink_minus_tnear = path.join(qcf, 'tlink_minus_tnear.txt')
    angles_txt = path.join(qcf, 'distortion_angles.txt')
    mid_surface = path.join(intf, 'mid_81920.obj')
    vertexmask = path.join(qcf, 'not_subplate_mask.txt')

    def minccalc_mask(label: int, mask_fname: str):
        """
        Create a binary mask for a given label number from the painted labels volume.
        :param label: desired label
        :param mask_fname: mask file prefix
        """
        sp.run(['minccalc', '-quiet', '-clobber', '-byte', '-unsigned', '-expr',
                f'A[0]>{label - 0.5}', segmentation_mnc, mask_fname], check=True)

    def surface_qc(surface, mask, chamfer, dist_txt, smth_txt, area_txt):
        if not qc:
            return
        sp.run(['chamfer.sh', '-c', '0.0', mask, chamfer], check=True)
        sp.run(['volume_object_evaluate', '-linear', chamfer, surface, dist_txt], check=True)
        sp.run(['smoothness.py', surface, smth_txt], check=True)
        sp.run(['depth_potential', '-area_simple', surface, area_txt], check=True)

    def run_log(*args, logfile_name=None, **kwargs):
        if qc:
            with open(logfile_name, 'w') as log_file:
                sp.run(*args, **kwargs, stderr=sp.STDOUT, stdout=log_file)
        else:
            sp.run(*args, **kwargs, stderr=sp.STDOUT, stdout=sp.DEVNULL)

    minccalc_mask(3, layer3_mask_mnc)
    minccalc_mask(4, layer4_mask_mnc)
    run_log(['marching_cubes_fetus.pl', f'-{side}', '-age', age,
             layer3_mask_mnc, layer3_obj], check=True, logfile_name=layer3_log)
    run_log(['fit_subplate.pl', '-age', age,
             layer4_mask_mnc, layer3_obj, layer4_obj], check=True, logfile_name=layer4_log)
    sp.run(['cortical_thickness', '-tlink', layer4_obj, layer3_obj, thickness_tlink], check=True)

    if not qc:
        return

    surface_qc(layer3_obj, layer3_mask_mnc, layer3_chamfer_mnc,
               layer3_dist_txt, layer3_smth_txt, layer3_area_txt)
    surface_qc(layer4_obj, layer4_mask_mnc, layer4_chamfer_mnc,
               layer4_dist_txt, layer4_smth_txt, layer4_area_txt)
    sp.run(['cortical_thickness', '-tlink',
            layer3_obj, layer4_obj, thickness_tnear], check=True)
    sp.run(['vertstats_math', '-old_style_file', '-sub',
            thickness_tlink, thickness_tnear, tlink_minus_tnear], check=True)
    sp.run(['distortion_angles.py', '-mid', mid_surface,
            layer4_obj, layer3_obj, angles_txt], check=True)
    sp.run(['diemesh.py', '-keep', path.join(intf, 'highlight_uncovered_subplate.mnc'),
            '-boundary', "2", segmentation_mnc, layer3_obj, vertexmask], check=True)
