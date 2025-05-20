### 1. Install the requirements
- required python version: 3.8 or 3.9
- run `pip install -r requirements.txt`
- for certain parts of the notebook (described in the notebook), a cuda GPU is required (torch with cuda support will be installed in requirements)

### 2. Download the data
- Download the associated data
- to reproduce the plots, 3 tar files are provided:
  - neuronal_responses.tar
  - natural_images.tar
  - model_checkpoints.tar
- copy these 3 .tar files into `.data/` directory
- extract the .tar files in the `.data/` directory, by using for example
  - `tar -xvf neuronal_responses.tar`

### 3. Run the notebook
Full Notebook for recreating the plots can be found under: [**bisc_nat_images.ipynb**](./bisc_nat_images.ipynb)
