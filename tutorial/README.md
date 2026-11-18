# bAIes-IDP-OpenMM

Using bAIes with openMM

## Quickstart

### 0. Run a structure prediction

You can use AlphaFold2, AlphaFold3, Boltz-1, Boltz-2, ESMFold or ESMFold2.

Collect a predicted model (.pdb) and the distograms (.pkl) and deposit them in a subdirectory in the `predictions` directory.

If the predicted model is in mmCIF format (.cif), you can convert it to PDB by running the following on the terminal. Make sure to replace the elements in <> brackets by the desired file paths.

```bash
scripts/cif2pdb.py <path_to_cif_file>.cif <path_to_output_file>.pdb
```

Ensure that your structure predictor outputs the distograms [here](link)

### 1. Prepare the inputs

Run on the terminal:

```bash
./1-prepare_bAIes.bash -pdb <path_to_prediction_pdb>.pdb \
                       -pkl <path_to_distograms>.pkl \
                       -out <path_to_output>
```

### 2. Run bAIes.

Run on the terminal:

```bash
./2-run_bAIes.bash -pdb <path_to_bAIes_pdb>.pdb \
                   -baies <path_to_baies_params>.dat \
                   -out <path_to_output>
```

