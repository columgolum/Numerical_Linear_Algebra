from setuptools import setup, find_packages

setup(
    name="dumbpy", # (or "NLA_athsh" depending on your package name)
    version="0.2.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "juliacall>=0.9.14",
    ],
    package_data={
        # This tells setuptools to include all .jl files and json files in all subfolders
        "": ["*.jl", "*/*.jl", "juliapkg.json"],
    },
    include_package_data=True,
    python_requires=">=3.7",
)
