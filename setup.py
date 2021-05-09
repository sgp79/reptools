"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

#To install:
#   py -3 setup.py sdist
#   pip3 install .

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
from os import path
from io import open
#from reptools import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

#Get the version
def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        version_file,
        re.M,
    )
    if version_match:
        return version_match.group(1)

    raise RuntimeError("Unable to find version string.")
    
    
# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='reptools', 
    version=open("reptools/version.py").readlines()[-1].split()[-1].strip("\"'"),
    
    # https://packagiATR01400 ng.python.org/specifications/core-metadata/#summary
    description='Tools for processing Rep-seq data',
    
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=long_description,
    
    # https://packaging.python.org/specifications/core-metadata/#description-content-type-optional
    long_description_content_type='text/markdown',
    
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    #url='',  # Optional

    author='Stephen Preston', 
    
    author_email='stephen.preston@zoo.ox.ac.uk', 
    
    # For a list of valid classifiers, see https://pypi.org/classifiers/
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',
    
        'Intended Audience :: Immunologists',
        'License :: OSI Approved :: MIT License',
    
        'Programming Language :: Python :: 3',
    ],
    
    # Note that this is a string of words separated by whitespace, not a list.
    #keywords='sample setuptools development',  # Optional

    #
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),  # Required
    
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=['numpy','numba'],
    python_requires='>=3.7',
    
    #extras_require={  # Optional
    #    'dev': ['check-manifest'],
    #    'test': ['coverage'],
    #},
    
    #package_data={  # Optional
    #    'sample': ['package_data.dat'],
    #},
    
    # Although 'package_data' is the preferred approach, in some case you may
    # need to place data files outside of your packages. See:
    # http://docs.python.org/3.4/distutils/setupscript.html#installing-additional-files
    #
    # In this case, 'data_file' will be installed into '<sys.prefix>/my_data'
    #data_files=[('my_data', ['data/data_file'])],  # Optional
    
    # The following provides a command called `reptools` which
    # executes the function `main` from the reptools.cli package when invoked:
    entry_points={ 
        'console_scripts': [
            'reptools=reptools.cli:main',
        ],
    },
    
    # List additional URLs that are relevant to your project as a dict.
    # https://packaging.python.org/specifications/core-metadata/#project-url-multiple-use
    #project_urls={  # Optional
    #    'Bug Reports': 'https://github.com/pypa/sampleproject/issues',
    #    'Funding': 'https://donate.pypi.org',
    #    'Say Thanks!': 'http://saythanks.io/to/example',
    #    'Source': 'https://github.com/pypa/sampleproject/',
    #},
)
