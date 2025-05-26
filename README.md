# Combinatorial INverse DESigner (CINDES)

This package allows for combinatorial optimizations, in which molecular properties are optimized by introducing substituents onto or core-modifications into the molecular framework.
The CINDES package containing two modules:

*   evaluation: HPC automation of job submission/processing
*   algorithms: Several algorithms to perform combinatorial optimizations. The most important ones are the Best-First-Search algorithm and several Genetic Algorithm implementations (as well as NSGA-II), and Bayesian optimization

## Information

Prof. Dr. Freija De Vleeschouwer

Algemene Chemie (ALGC) – General Chemistry, Vrije Universiteit Brussel (VUB)

Main developer: Jos L. Teunissen

Contributors: Eline Desmedt, David Smets

Email correspondence: freija.de.vleeschouwer@vub.be

## Related publications

*   Jos L. Teunissen, Frank De Proft, Freija De Vleeschouwer, Tuning the HOMO-LUMO Gap of Small Diamondoids Using Inverse Molecular Design. _J. Chem. Theory Comput_. **2017**, _13_, 1351.
*   Jos L. Teunissen, Frank De Proft, Freija De Vleeschouwer, Acceleration of Inverse Molecular Design by Using Predictive Techniques. _J. Chem. Inf. Model_. **2019**, _59_, 2587.
*   Eline Desmedt, Tatiana M. Woller, Jos Teunissen, Freija De Vleeschouwer, Mercedes Alonso, Fine-tuning of Nonlinear Optical Contrasts of Hexaphyrin-based Molecular Switches using Inverse Design. _Front. Chem_. **2021**, _9_, 786036.
*   Eline Desmedt, David Smets, Tatiana Woller, Mercedes Alonso, Freija De Vleeschouwer, Designing hexaphyrins for high-potential NLO switches: The synergy of core-modifications and meso-substitutions. _Phys. Chem. Chem. Phys_. **2023**, _25_, 17128.
*   Eline Desmedt, Léa S. Gimenez, Freija De Vleeschouwer, Mercedes Alonso,	Application of Inverse Design Approaches to the Discovery of Nonlinear Optical Switches. _Molecules_ **2023**, _28_, 7371.
*   Irene Casademont-Reig, Roger Monreal-Corona, Eline Desmedt, Freija De Vleeschouwer, and Mercedes Alonso. Pursuit of Singlet Fission Fulvenes Candidates using Inverse Design. Submitted to Digital Discovery.

## Contents

*   [Getting started](#getting-started)
    *   [Install](#install)
    *   [Requirements](#requirements)
    *   [Usage](#usage)
    *   [Output data](#output)
    *   [Failed CINDES run](#failed)
*   [Examples](#examples)

## Getting started

### Install

Use git to clone this repository into the **Project** directory (located in your working directory).

```
mkdir Project
cd Project
git clone https://gitlab.com/fdevlees/cindes.git
```

### Requirements

Load python module, here _matplotlib/3.9.2-gfbf-2024a_.

```
module load matplotlib/3.9.2-gfbf-2024a
```

Additional packages are required via setup.py.

```
cd cindes/
pip install .

```

### Usage

Before running the script, make sure the following files are in the **Project/cindes/** directory:

*   **INPUT** (input with keywords required for running CINDES) and **ZMAT** (Z-matrix of your molecule). For examples, see **/Examples/**

Make the following scipts executable:

*   **ID_gauss**, your computational software submission script
*   **cindes_submit**, potentially your CINDES submission script, where PYTHONPATH refers to the _Project/cindes/_ directory.

```
chmod +x ID_gauss
chmod +x cindes_submit
```

Also, make the script **qsta** generally available and executable, for example in a bin/ folder in your **$HOME** directory (Note: adapt $HOME accordingly):
```
mkdir  $HOME/bin/
cp qsta $HOME/bin/qsta
chmod +x $HOME/bin/qsta
```

In $HOMO/.bashrc, you add the following:

```
export PATH=$HOME/bin:$PATH
```
 
Also, include the path in **cindes_submit**. 
The **qsta** script collects relevant data from the slurm scheduler, needed to verify whether calculations are still running, have failed or have completed.

Modify **ID_gauss**, the Gaussian software submission script, for your usage.

Run the module.

*   Via command line:

```
python -m CINDES -i INPUT > cindes_output.log
```

*   Via submission script, here using slurm workload manager. Make sure to adjust the paths in the **cindes_submit** script.

```
sbatch cindes_submit
```

### Output data

All input and output files of the computational chemistry software calculations are collected in the newly created directory **CALC**, whereas errors are summarized in **.err** files in the directory **logs**.

A **table.json** database is created containing the data from structures of previous runs. As such, restarts first check the database before submitting computational software calculations.

The file **cyclesinfo** is made available with all molecular data from the inverse design run.

Output file **cindes_output.log**, finally, contains detailed information and data about the integral inverse design process (iteration and steps within) and ends with a listing of the optimum. 

### Failed CINDES run

Gaussian software log files are checked for error terminations or imaginary frequencies (in case of opt+freq), and the user is warned. In these cases, failed calculations need to be restarted manually, after which a restart of the CINDES program can be set up. The restart requires the initial structure and the site sequences. Do not remove any calculations from the **CALC** directory.

## Examples

Examples are given in the **Examples** directory.

