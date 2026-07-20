import os
import subprocess
import shutil
from setuptools import setup, find_packages
from setuptools.command.install import install

class JuliaInstallCommand(install):
    def run(self):
        # 1. Check if julia exists
        if shutil.which("julia"):
            print("Julia found.")
        else:
            print("Julia not found, downloading Julia...")
            # Logic to download/install Julia for the user
            # subprocess.run(["...download_script..."]) 
        
        # 2. Check/Bridge juliacall
        print("Bridging julia and python")
        # Ensure juliacall is installed/configured
        subprocess.check_call(["pip", "install", "juliacall"])
        
        install.run(self)

setup(
    name="NLA_athsh",
    version="0.2.0",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    cmdclass={'install': JuliaInstallCommand},
    install_requires=["juliacall"],
    python_requires=">=3.7",
)
