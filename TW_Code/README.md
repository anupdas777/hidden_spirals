### 1. Install the following python requirements and packages
- required python version: 3.11.0
- pip install numpy
- pip install pandas
- pip install matplotlib
- pip install scipy
- pip install scikit-learn
- pip install mpl-tools
- pip install pylablib

### 2. Simulation example 1
- Simulation example 1 namely "wave_spatial_autocorrelation_ica.ipynb" is a jupyterlab notebook which provides an example simulation demonstrating the spatial autocorrelation preserving surrogate data analysis procedure employed in our paper. 
- This example simulates a rotational traveling wave, which ubiquitous in our experimental data.
- It also simulates 500 permutations of circularly shifted (vertical or horizontal or both) traveling waves of the original rotational wave. Note that this circular shuffling preserves the spatial autocorrelation of the rotational traveling wave. 
- These shuffled waves are then passed onto complex ICA to estimate the modes. These are then compared to the original rotational wave to compute a distribution or histogram of error magnitudes. The original rotational wave is similarly passed onto ICA and the empirical error is computed. The error in the later case is significantly lower compared to the error for the circularly shuffled waves.

### 3. Simulation example 2
- Simulation example 2 namely "wave_experiment_stable_epochs_ica.ipynb" is a jupyterlab notebook which provides an example of how the ICA procedure works on stable epochs corresponding to the experimental data. 
- This example first loads the stable epochs UP021_vf_pres.pkl, the location of the electrodes UP021_xyz_grid.pkl, the initialization matrix for ICA UP021_init_mat.pkl, the letter labels UP021_letter_list.pkl, the unique letters UP021_items_upper_unique.pkl, and PCA transformed electrode coordinates UP021_xg.pkl and UP021_yg.pkl.  
- The stable epochs are then passed onto ICA and then we extract and plot the modes extracted from the ICA. 
