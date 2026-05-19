# s-SNOM_spectra
An open-source, standalone MATLAB application for calculating the near-field spectrum of bulk homogeneous uniaxial samples based on the spheroid model. The repository also includes a Jupyter notebook and MATLAB script versions for flexible and user-friendly workflows. The model is presented in the paper K.V. Voronin et al. "Quantitative Analytical Spheroid Model for Scattering-type Scanning Near-field Optical Spectroscopy", Advanced Optical Materials (2025). The authors thank Edoardo Vicentini for his help with optimizing the calculations and translating the script from MATLAB to Python.

The folder "MATLAB App project" contains all the files required to compile SNOM_Spectra.exe or s-SNOM Installer.exe. The application s-SNOM Installer.exe is compiled in MATLAB 2025a and does not require a MATLAB license. It allows Users to calculate the near-field spectrum of bulk homogeneous uniaxial samples based on the spheroid model of the tip. Users can set all geometrical parameters and dielectric permittivities of the materials.

Spheroid_model_spectra.m was used for calculations presented in the mentioned paper; spheroid_model_spectra.ipynb is the translation of the MATLAB script to Python.

## Web Calculator

The repository now also contains a first web version of the bulk s-SNOM calculator under `snom_web/`.

Implemented in the first version:

- bulk sample mode;
- Fast and Accurate calculation modes;
- configurable tip, reference, far-field and frequency-range parameters;
- built-in materials from the MATLAB application and uploaded permittivity files;
- uploaded experimental spectrum overlay;
- three interactive Plotly graphs;
- results export to CSV.

Future-oriented structure already exists for:

- `Layered sample` mode;
- `Drude-Lorentz model` material input.

### Run locally

The web app uses only Python's standard HTTP server stack plus `numpy`.

```bash
python3 -m snom_web --host 127.0.0.1 --port 8000
```

Then open `http://127.0.0.1:8000` in a browser.
