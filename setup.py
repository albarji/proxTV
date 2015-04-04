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
    version="3.1",
    packages=['prox_tv'],
    install_requires=[
        'numpy',
        'cffi',
    ],
    setup_requires=[
        'cffi'
    ],
    cmdclass={
        'build': CFFIBuild,
        'install': CFFIInstall,
    },
    author="Alvaro Barbero, Suvrit Sra",
    author_email="alvaro.barbero@uam.es",
    url='https://github.com/albarji/proxTV',
    license='BSD', 
    test_suite="nose.collector",
    zip_safe=False,
)
