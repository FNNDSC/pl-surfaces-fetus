from os import path
from glob import glob
from setuptools import setup

# CIVET's container image us an older Ubuntu
# its preferred encoding is ANSI_X3.4-1968, yuck
# gotta use encoding='utf_8' for fancy symbols
with open(path.join(path.abspath(path.dirname(__file__)), 'README.md')) as f:
    readme = f.read()

setup(
    name='surfaces_fetus',
    version='1.0',
    description='Extract surfaces (.obj) for the subplate zone from pre-segmented fetal brain MRI '
                'volume (.mnc) using CIVET marching-cubes and surface_fit',
    long_description=readme,
    author='Jennings Zhang',
    author_email='Jennings.Zhang@childrens.harvard.edu',
    url='https://fnndsc.childrens.harvard.edu/conferences/2020/OHBM/Jennings/'
        'Jennings_Zhang_OHBM_2020_Subplace_Surfaces.pdf',
    packages=['surfaces_fetus'],
    install_requires=['chrisapp~=1.1.6', 'pybicpl==0.1-1'],
    license='MIT',
    zip_safe=False,
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'surfaces_fetus = surfaces_fetus.__main__:main'
        ]
    },
    scripts=glob('scripts/*')
)
