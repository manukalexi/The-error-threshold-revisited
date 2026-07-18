Repository for the software described in

# The error threshold revisited: Sustainable information vs. sustainable concentrations

Alexa Iván<sup>1,2,4</sup>; Tamás Czárán<sup>1</sup>; Eörs Szathmáry<sup>1,2,3</sup>; András Szilágyi<sup>1,2</sup>

<sup>1</sup> Institute of Evolution, HUN-REN Centre for Ecological Research, Budapest, Hungary  
<sup>2</sup> Center for the Conceptual Foundations of Science, Parmenides Foundation, Pöcking, Germany  
<sup>3</sup> Department of Plant Systematics, Ecology and Theoretical Biology, ELTE Eötvös Loránd University, Budapest, Hungary  
<sup>4</sup> Doctoral School of Biology, ELTE Eötvös Loránd University, Budapest, Hungary

## Models

In this study three main models were developed and analyzed:

- **Genotype model**
- **Phenotype model**
- **Simplified model**

The Genotype model and the Simplified model are investigated analytically using systems of ordinary differential equations (ODEs), while the Phenotype model is studied through Monte Carlo simulations within an agent-based modeling framework.

The article presents results for the following specific cases of these models:

- **Genotype model: $L=25,  50,  100$**
- **Phenotype model: $L=25$**
  - small loop
  - large loop
  - motif
- **Simplified model: $L=25$**

The corresponding source codes and makefiles for each model are provided in the respective directories.

## Installation and Usage

1. Download the following directories containing the C source codes and makefiles:

- Genotype model
- Phenotype model
- Simplified model

2. On Linux/Ubuntu systems with the GCC toolchain installed, navigate to the directory or subdirectory of the selected model and compile the code using:

```bash
make
```

3. Run the compiled executable corresponding to the selected model:

```bash
./genotype_model
```
```bash
./small_loop
```
```bash
./large_loop
```
```bash
./motif
```
```bash
./simplified_model
```
4. For the Genotype and Simplified models, save the program output to a file using output redirection (`>`):

```bash
./genotype_model > output_filename.txt
```
```bash
./simplified_model > output_filename.txt
```
> [!NOTE]
> For the Phenotype model, specifying an output file is not necessary, as the program creates the output files automatically.

5. Use the raw output files for further data evaluation and visualization.

## Parameters of the models

### Genotype model

- `L`: sequence length (default: 25, 50, 100)
- `N`: number of Hamming classes (default: L+1)
- `mu`: mutation rate
- `sigma_master`: reproductive superiority of the master sequence (default: 10)
- `n_scal`: fitness landscape geometry factor (investigated: 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.3, 1.7, 2.3, 3.0, 6.0)
- `tmax`: maximum integration time (default: 10000)
- `dt_out`: output time interval (default: 0.01)

### Phenotype model
- `L`: genome length (default: 25)
- `N`: population size (default: 100000)
> [!CAUTION]
> If you modify `N`, you must provide a master sequence of length `N` that satisfies the requirements.
- `mu`: mutation rate
- `sigma`: reproductive superiority of the master sequence (default: 10)
- `GEN`: number of generations (default: 10000000)
> [!WARNING]
> If `GEN` is too low, the simulation may terminate before the population reaches its equilibrium state.
- `REP`: number of parallel simulations (default: 50)
- `n`: fitness landscape geometry factor (investigated: 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.3, 1.7, 2.3, 3.0, 6.0)
- `K`: maximum error score, reducing fitness to the baseline level (default: 7.0)
> [!CAUTION]
> Do not change `K`! It will mess up the determination of fitness and alter the population dynamics!

- **Simplified model**
- `L`: sequence length (default: 25)
- `N`: number of Hamming classes (default: 3)
- `mu`: mutation rate
- `sigma`: reproductive superiority of the master sequence (default: 10)
- `lambda`: fitness of the 1-error mutant (investigated: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

## Sample output


## Disclaimer
The code provided in this repository is released for research purposes only and is provided without warranty or liability.

