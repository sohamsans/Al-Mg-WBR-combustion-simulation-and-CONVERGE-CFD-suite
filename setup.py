from setuptools import setup, find_packages

setup(
    name="almg-wbr-combustion-suite",
    version="1.0.0",
    description="Al-Mg Composite Solid Propellant Combustion and CONVERGE CFD Simulation Suite",
    author="Soham Sans",
    url="https://github.com/sohamsans/Al-Mg-WBR-combustion-simulation-and-CONVERGE-CFD-suite",
    py_modules=[
        "cea_wbr_final_lab",
        "primary_combustion_model",
        "bates_converge_cfd_exporter",
        "optimize_expulsion",
        "test_physics_and_math"
    ],
    install_requires=[
        "rocketcea>=1.1.25",
        "numpy>=1.26.0",
        "scipy>=1.12.0",
        "pandas>=2.2.0",
        "matplotlib>=3.8.0"
    ],
    entry_points={
        "console_scripts": [
            "wbr-lab-gui=cea_wbr_final_lab:main",
            "wbr-cfd-export=bates_converge_cfd_exporter:simulate_bates_grain_regression",
            "wbr-physics-test=primary_combustion_model:simulate_1d_primary_combustor"
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Operating System :: OS Independent"
    ],
    python_requires=">=3.10",
)
