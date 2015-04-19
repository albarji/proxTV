from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build

# We need this workaround to solve the problem of importing prox_tv._prox_tv
# before CFFI is installed. For more details, have a look at:
#   https://caremad.io/2014/11/distributing-a-cffi-project/


def get_ext_modules():
    from prox_tv._prox_tv import _ffi
    return [
        _ffi.verifier.get_extension(),
    ]


class CFFIBuild(build):
    def finalize_options(self):
        self.distribution.ext_modules = get_ext_modules()
        build.finalize_options(self)


class CFFIInstall(install):
    def finalize_options(self):
        self.distribution.ext_modules = get_ext_modules()
        install.finalize_options(self)

setup(
    name="prox_tv",
    version="3.1.0",
    description="Toolbox for fast Total Variation proximity operators",
    long_description="proxTV is a toolbox implementing blazing fast implementations of Total Variation proximity operators. While the core algorithms are implemented in C to achieve high efficiency, Matlab and Python interfaces are provided for ease of use. The library provides efficient solvers for a variety of Total Variation proximity problems, with address input signals of any dimensionality (1d, images, video, ...) and different norms to apply in the Total Variation term.",
    packages=['prox_tv'],
    install_requires=[
        'numpy',
        'cffi',
        'sphinxcontrib-napoleon',
        'sphinx_rtd_theme',
        'matplotlib',
        'scipy',
        'scikit-image'
    ],
    setup_requires=[
        'cffi'
    ],
    cmdclass={
        'build': CFFIBuild,
        'install': CFFIInstall,
    },
    package_data={
        'prox_tv': ['src/*.h', 'src/demos/*']
    },
    author="Alvaro Barbero, Suvrit Sra",
    author_email="alvaro.barbero@uam.es",
    url='https://github.com/albarji/proxTV',
    license='BSD', 
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Mathematics',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.4'
    ],
    keywords='total variation image processing machine learning',
    test_suite="nose.collector",
    zip_safe=False,
)
