""" Optimal Enzyme Kinetics ""
"Optimal enzyme utilization suggests concentrations and thermodynamics
favor condition-specific saturations and binding mechanisms

.. moduleauthor:: OPEN team

"""

from setuptools import setup, find_packages

version_tag = '0.0.1'


setup(name='open',
      version=version_tag,
      author='OPEN team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/open/',
      download_url='https://github.com/EPFL-LCSB/open/archive/'+version_tag+'.tar.gz',
      install_requires=['jedi==0.17.2', #for ipython in python 3.6
                        'ipython',
                        'ipdb',
                        'lxml',
                        'openpyxl',
                        'tabulate',
                        'sphinx>=1.4',
                        'sphinx-rtd-theme',
                        'pytest',
                        'numpy',
                        'pandas',
                        'tables',
                        'optlang',
                        'matplotlib',
                        'numexpr==2.7.1',
                        'pipdeptree==1.0.0',
                        'cobra==0.17.1',
                        'wheel',
                        'tqdm',
                        'sklearn',
                        ],
      packages = find_packages(),
      python_requires='>=3, <4',
     # python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='Optimal Enzyme Kinetics',
      keywords=['kinetic','models','optimal'],
      license='Apache2',

      # See https://PyPI.python.org/PyPI?%3Aaction=list_classifiers
      classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 3 - Alpha',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry'
            'Environment :: Console',

            # Pick your license as you wish (should match "license" above)
            'License :: OSI Approved :: Apache Software License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 3.6',
      ],
     )
