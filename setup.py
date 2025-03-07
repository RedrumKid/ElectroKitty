from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages

__version__ = "1.0.11.1"

ext_modules = [
    Pybind11Extension(
        "cpp_ekitty_simulator",
        ["./electrokitty/cpp_code/electrokitty_simulator.cpp", "./electrokitty/cpp_code/electokitty_helper_file.cpp"],
        define_macros=[("VERSION_INFO", __version__)],
    ),
]

setup(
    name="electrokitty",
    version=__version__,
    author="Ožbej Vodeb",
    author_email="ozbej.vodeb@gmail.com",
    url="https://github.com/RedrumKid/ElectroKitty",
    description="A simulation and simple analysis tool for electrochemical data",
    packages = find_packages(),
    ext_modules=ext_modules,
    install_requires=[
        "numpy >= 1.20.0", 
        "scipy >= 1.11.3",
        "matplotlib >= 3.8.0",
        "cma >= 3.3.0"
        ],
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
    license="BSD"
)
