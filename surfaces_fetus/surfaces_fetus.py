#!/usr/bin/env python                                            
#
# surfaces_fetus ds ChRIS plugin app
#
# (c) 2016-2020 Fetal-Neonatal Neuroimaging & Developmental Science Center
#                   Boston Children's Hospital
#
#              http://childrenshospital.org/FNNDSC/
#                        dev@babyMRI.org
#
import pkg_resources
from chrisapp.base import ChrisApp
from .script import process, UserError

Gstr_title = """
                 __                       __     _             
                / _|                     / _|   | |            
 ___ _   _ _ __| |_ __ _  ___ ___  ___  | |_ ___| |_ _   _ ___ 
/ __| | | | '__|  _/ _` |/ __/ _ \/ __| |  _/ _ \ __| | | / __|
\__ \ |_| | |  | || (_| | (_|  __/\__ \ | ||  __/ |_| |_| \__ \
|___/\__,_|_|  |_| \__,_|\___\___||___/ |_| \___|\__|\__,_|___/
                                    ______                     
                                   |______|                    
"""

Gstr_synopsis = """
    NAME

       surfaces_fetus.py 

    DESCRIPTION

        Extract surfaces (.obj) for the subplate zone from pre-segmented
        fetal brain MRI volume (.mnc) using CIVET marching-cubes and surface_fit
        
    USAGE
    
        python3 -m surfaces_fetus --age <n.n> --side [left|right] <inputDir> <outputDir>

        The input directory should contain a single .mnc file, which is the
        painted labels i.e. segmented volume of a single brain hemisphere (left or right)
        from a fetal brain MR image. The volume should have these labels:
        
            3=white matter
            4=intermediate zone
        
        The script has two required arguments: --age (in gestational weeks)
        and --side (left or right brain hemisphere).
    
    EXAMPLE
    
        A left-side 30.1 GA (gestaional age) brain segmentation can be found at
        /neuro/data/s1_0090/labels.mnc
        This is the input file.
        Desired output directory is /neuro/results/s1_0090/surfaces/
        
        IN=/neuro/data/s1_0090/
        OUT=/neuro/results/s1_0090/surfaces/
        docker run --rm -v $IN:/incoming:ro $OUT:/outgoing fnndsc/pl-surfaces-fetus \
            surfaces_fetus --age 30.1 --side left /incoming /outgoing

    CHRIS PLUGIN OPTIONS

        [-h] [--help]
        If specified, show help message and exit.
        
        [--json]
        If specified, show json representation of app and exit.
        
        [--man]
        If specified, print (this) man page and exit.

        [--meta]
        If specified, print plugin meta data and exit.
        
        [--savejson <DIR>] 
        If specified, save json representation file to DIR and exit. 
        
        [-v <level>] [--verbosity <level>]
        Verbosity level for app. Not used currently.
        
        [--version]
        If specified, print version number and exit. 
"""


class SurfacesFetusWrapper(ChrisApp):
    """
    Extract surfaces (.obj) for the subplate zone from pre-segmented fetal
    brain MRI volume (.mnc) using CIVET's marching-cubes and surface_fit.
    """
    AUTHORS                 = 'Jennings Zhang <Jennings.Zhang@childrens.harvard.edu>'
    SELFPATH                = '/usr/local/bin'
    SELFEXEC                = 'surfaces_fetus'
    EXECSHELL               = 'python'
    TITLE                   = 'Surface extraction for fetal brain MRI'
    CATEGORY                = ''
    TYPE                    = 'ds'
    DESCRIPTION             = 'Extract surfaces (.obj) for the subplate zone from ' \
                              'pre-segmented fetal brain MRI volume (.mnc) using CIVET marching-cubes and surface_fit'
    DOCUMENTATION           = 'https://fnndsc.childrens.harvard.edu/conferences/2020/OHBM/Jennings/' \
                              'Jennings_Zhang_OHBM_2020_Subplace_Surfaces.pdf'
    VERSION                 = pkg_resources.require('surfaces_fetus')[0].version
    ICON                    = '' # url of an icon image
    LICENSE                 = 'Opensource (MIT)'
    MAX_NUMBER_OF_WORKERS   = 1  # Override with integer value
    MIN_NUMBER_OF_WORKERS   = 1  # Override with integer value
    MAX_CPU_LIMIT           = '' # Override with millicore value as string, e.g. '2000m'
    MIN_CPU_LIMIT           = '' # Override with millicore value as string, e.g. '2000m'
    MAX_MEMORY_LIMIT        = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_MEMORY_LIMIT        = '' # Override with string, e.g. '1Gi', '2000Mi'
    MIN_GPU_LIMIT           = 0  # Override with the minimum number of GPUs, as an integer, for your plugin
    MAX_GPU_LIMIT           = 0  # Override with the maximum number of GPUs, as an integer, for your plugin

    # Use this dictionary structure to provide key-value output descriptive information
    # that may be useful for the next downstream plugin. For example:
    #
    # {
    #   "finalOutputFile":  "final/file.out",
    #   "viewer":           "genericTextViewer",
    # }
    #
    # The above dictionary is saved when plugin is called with a ``--saveoutputmeta``
    # flag. Note also that all file paths are relative to the system specified
    # output directory.
    OUTPUT_META_DICT = {}

    def define_parameters(self):
        """
        Define the CLI arguments accepted by this plugin app.
        Use self.add_argument to specify a new app argument.
        """
        self.add_argument('--age', dest='age', type=float, optional=False,
                          help='gestational age estimate in weeks')
        self.add_argument('--side', dest='side', type=str, optional=False,
                          help='brain hemisphere [left, right]')
        self.add_argument('--keep-intermediate', dest='keep', type=bool, default=False, optional=True,
                          help='keep intermediate files (e.g. *mask.mnc, *chanfer.mnc)')
        self.add_argument('--qc', dest='qc', type=bool, default=False, optional=True,
                          help='save surface_fit logs and produce vertex-wise quality check files')

    def run(self, options):
        """
        Define the code to be run by this plugin app.
        """
        try:
            process(options.inputdir, options.outputdir, options.side, options.age, options.keep, options.qc)
        except UserError as e:
            print(e.message)

    def show_man_page(self):
        """
        Print the app's man page.
        """
        print(Gstr_synopsis)
