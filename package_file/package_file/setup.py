from setuptools import setup, find_packages

setup(
    name="NLA_atsh",
    version="0.1.0",
    description="High-Performance Hybrid Python-Julia Numerical Linear Algebra Package",
    author="Atharv Shrivastav",
    # Yeh line setuptools ko batati hai ki saare packages 'src' folder ke andar hain
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    # Aapke package ko chalne ke liye juliacall zaroori hai
    install_requires=[
        "juliacall",
    ],
    python_requires=">=3.7",
) 