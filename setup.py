""" GEneralized Elementary Kinetics'

.. moduleauthor:: geek team


"""

from setuptools import setup, find_packages

version_tag = '0.0.1'

setup(name='geek',
      version=version_tag,
      author='geek team',
      author_email='softwares.lcsb@epfl.ch',
      url='https://github.com/EPFL-LCSB/geek/',
      download_url='https://github.com/EPFL-LCSB/geek/archive/'+version_tag+'.tar.gz',
      install_requires=['sympy >= 1.1.',
                        'pytest',
                        'scipy',
                        'pandas',
                        'statsmodels',
                        'openbread',],

      packages = find_packages(),
      python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, <4',
      description='',
      keywords=['crowding','spatial-effects','kinetic','models'],

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
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
      ],
     )
