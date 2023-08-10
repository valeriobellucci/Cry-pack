from setuptools import setup, find_packages

setup(
    name="crysty",
    version="0.1",
    author="Valerio",
    #author_email="valerio@email.com",  # Replace with your email
    description="Crystalline X-ray Optic package",  # Optional
    long_description=open('README.md').read(),  # Optional
    long_description_content_type='text/markdown',  # Optional
    url="https://github.com/valeriobellucci/crysty",  # Replace with your package's GitHub repo
    packages=find_packages(),
    # List your package's dependencies here
    install_requires=[
        "numpy", "quantities", "scipy", "matplotlib", 
        "scikit-image", "pandas", "pyopencl", "setuptools"
    ],
)
