This is very much still a work in progress. 

## Section 1: Loading and Parsing StoBe Output Files

### Overview

The data loading pipeline begins with parsing StoBe-deMon output files, which contain the results of DFT calculations for core-excited states. StoBe calculations typically generate several types of output files for each atomic site in the molecule:

Ground State (GND): gnd.out files containing ground state electronic structure
Excited State (EXC): exc.out files with core-excited state information
Transition Potential (TP): tp.out files containing results from transition potential calculation (contains the xray transitions)

### File Structure and Organization

The program expects a specific directory structure and file naming scheme where StoBe output files are organized in folders by calculation type:

```text
project_directory/
├── GND/
│   ├── atom1_gnd.out
│   ├── atom2_gnd.out
│   └── ...
├── EXC/
│   ├── atom1_exc.out
│   ├── atom2_exc.out
│   └── ...
├── TP/
│   ├── atom1_tp.out
│   ├── atom2_tp.out
│   └── ...
```

This can be done via the snippet below. 

```python
from stobeLoader import load_data

# Example usage function
dft_data_folder = r'E:\CuPc\CuPc DFT\CuPc_2005_Isolated\Large Basis Sets\UHF'
width1 = 0.6 / 2.355  # Convert FWHM to sigma
width2 = 12.0 / 2.355
max_energy = 320.0

# Load data
data = load_data(
    directory=dft_data_folder,
    width1=width1,
    width2=width2,
    max_energy=max_energy,
    verbose=True,
    n_jobs=8  # Use n parallel processes
)

if data is not None:
    # Access the different dataframes
    energy_results = data['energy_results']
    xray_transitions = data['xray_transitions']
    orbital_alpha = data['orbital_alpha']
    orbital_beta = data['orbital_beta']
    atomic_coordinates = data['atomic_coordinates']
    basis_sets = data['basis_sets']
else:
    print("Failed to load data")
```

That produces a variety of dataframes with useful information about the current molecule extracted from the StoBe output files. 

The 2 most important ones for now are: 
1. energy_results which contains the ionization potentials and energy corrections for each excitation center
2. xray_transitions which contains the transition info like the energies, oscillator strengths, widths, thetas, tensor components etc. **This is the dataframe you need for the clustering algorithm portion**

For now, save those dataframes into .csv for usage in the next part (I'll allow other file types later like .parquet)

```python
energy_results.to_csv(r'data\CuPc_energy_results.csv')
xray_transitions.to_csv(r'data\CuPc_xray_transitions.csv')
```

## Section 2: Running the Clustering Pipeline

The pipeline loads the xray_transitions and performs the clustering iteratively until the overlap matrix has no elements that exceed the overlap threshold. It also generates the DFT Absorption Step Edge from the DFT Ionization Potentials and scales it to the Molecule Absorption from the Henke data. It also loads and scales the experimental data to the Molecule Absorption.

At this stage, the pipeline takes the energy_results and xray_transitions dataframes as input, as well as the experimental AR-NEXAFS. You'll also enter the OS and OVP thresholds. 

The experimental data should consist of intensity and energy arrays. It should have a column naming scheme like {molName}_{sample_theta[i]} for each NEXAFS intensity array and E_{sample_theta[i]}. In other words, every intensity wave should have its corresponding energy wave.

```python
from nexafs_analyzer import NEXAFSAnalyzer

data_path = r'data\CuPc_xray_transitions2.csv'
df_step_path = r'data\CuPc_energy_results.csv'
exp_nxfs_path1 = r'data\CuPc_CuI_EXP.csv'
exp_nxfs_path2 = r'data\CuPc_Si_EXP.csv'

analyzer = NEXAFSAnalyzer(data_path=data_path)
dft_clustering_results = analyzer.run_complete_analysis(
    os_threshold_percent=2,
    overlap_threshold=50,
    compound_dict={'C': 32, 'H': 16, 'Cu': 1, 'N': 8},
    dft_step_data_path= df_step_path
)
```

You can produce a variety of visualizations as well. 

```python
from plotUtils import plot_iteration_nexafs

# Plot all iterations
plot_iteration_nexafs(dft_clustering_results, 'all')

```

<img width="1189" height="790" alt="image" src="https://github.com/user-attachments/assets/e1dc7418-3518-4f1c-b8ef-edee6ce5d541" />

```python
from plotUtils import plot_overlap_heatmaps

# Plot all iterations
plot_overlap_heatmaps(dft_clustering_results, 'all')
```

<img width="2385" height="495" alt="image" src="https://github.com/user-attachments/assets/c29fa2e5-54ea-4d92-86ed-7638714d201c" />

```python
from plotUtils import plot_gaussian_fits

# Plot all clusters (up to 30)
plot_gaussian_fits(dft_clustering_results, iteration=1, peaks='all')
```

<img width="1991" height="1475" alt="image" src="https://github.com/user-attachments/assets/15a58b12-7092-408f-83c4-35168e423f45" />

```python
from plotUtils import plot_bare_atom_absorption, plot_absorption_vs_step_edge, plot_step_edge_fit_quality

# Plot bare atom absorption
plot_bare_atom_absorption(dft_clustering_results, energy_xlim=(10, 30000))
```

<img width="1189" height="590" alt="image" src="https://github.com/user-attachments/assets/7632712b-26e7-4fa3-9add-885684a54b42" />

<img width="890" height="530" alt="image" src="https://github.com/user-attachments/assets/d2635b6a-100e-40b7-9a6a-9feefbfa8945" />

<img width="989" height="590" alt="image" src="https://github.com/user-attachments/assets/e261ab9b-91e2-421a-88eb-61248e19cd93" />

<img width="989" height="590" alt="image" src="https://github.com/user-attachments/assets/b19467a6-6904-4758-a31f-92e00946d427" />

## Section 3: Running the Simultaneous Fits

(Very much a work in progress still)

To get started with the simultaneous fits, you'll need the experimental data and a .csv with the peak parameters from the clustering routine. 

You can get that data like so: results['clustering']['iteration_results']['cluster_gaussians'][0].to_csv('PARAMETERS.csv')

```python
from multiSpecFitProcs import NEXAFSTensorFitter, analyze_tensor_evolution
# Create and fit data
fitter = NEXAFSTensorFitter("PARAMETERS.csv", "data/CuPc_CuI_EXP.csv")
results = fitter.fit_spectra()

# Analyze results
fitter.print_fit_summary()
fitter.plot_results()
fitter.plot_tensor_analysis()
analyze_tensor_evolution(fitter)

# Export results
fitter.export_results("fitted_parameters.csv")
```

Some initial results below. Again, work in progress...
<img width="1189" height="990" alt="image" src="https://github.com/user-attachments/assets/a1950ed7-8467-46cd-a592-3ed3175057d2" />

