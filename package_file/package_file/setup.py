from setuptools import setup, find_packages

setup(
    name="NLA_athsh",
    version="0.2.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "juliacall>=0.9.14",
    ],
    package_data={
        "dumbpy": ["juliapkg.json"],
    },
    python_requires=">=3.7",
)
