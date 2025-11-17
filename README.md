# KPP - The Kinetics PreProcessor

**Version 2.2.3**

KPP is a software tool that facilitates the numerical simulation of chemical kinetic systems. It automatically generates simulation code from symbolic chemical mechanism descriptions, supporting multiple programming languages and numerical integration methods.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Guide](#usage-guide)
  - [Input File Structure](#input-file-structure)
  - [Running KPP](#running-kpp)
  - [Compiling Generated Code](#compiling-generated-code)
- [Configuration Options](#configuration-options)
- [Examples](#examples)
- [Advanced Features](#advanced-features)
- [Troubleshooting](#troubleshooting)
- [Documentation](#documentation)
- [License](#license)
- [Support](#support)

## Overview

KPP (Kinetics PreProcessor) is a symbolic chemistry preprocessor that:

- Parses chemical mechanism definitions (species and reactions)
- Generates optimized simulation code in your target language
- Provides multiple numerical integration methods
- Supports sensitivity analysis, adjoint models, and stochastic simulations
- Enables efficient simulation of atmospheric chemistry, combustion, and other kinetic systems

**Project Website:** http://www.cs.vt.edu/~asandu/Software/KPP

## Features

### Supported Output Languages
- Fortran 77
- Fortran 90/95
- C
- MATLAB

### Numerical Integrators
- **Rosenbrock methods**: Adaptive, L-stable integrators for stiff systems
- **Runge-Kutta methods**: Explicit and implicit methods
- **SDIRK**: Singly Diagonally Implicit Runge-Kutta
- **LSODE/LSODES**: Livermore Solver for Ordinary Differential Equations
- **RADAU5**: Implicit Runge-Kutta method of order 5
- **SEULEX**: Extrapolation-based method
- **Gillespie/Tau-Leap**: Stochastic simulation algorithms

### Advanced Capabilities
- Sparse and full Jacobian computation
- Hessian calculation
- Forward sensitivity analysis (Tangent Linear Model)
- Adjoint sensitivity analysis
- Stochastic simulation
- Mass balance checking
- Customizable initial conditions and parameters

## Prerequisites

### Required
- **FLEX** (Fast Lexical Analyzer): Version 2.5 or higher
  - Test with: `flex --version`
  - Install on Debian/Ubuntu: `sudo apt-get install flex`
  - Install on macOS: `brew install flex`
  - Download from: http://www.gnu.org/software/flex/

- **C Compiler**: GCC or native C compiler
  - Test with: `gcc --version`

### Optional (depending on target language)
- **Fortran Compiler**: For Fortran output (gfortran, ifort, etc.)
- **MATLAB**: For MATLAB output and execution

## Installation

### Step 1: Verify FLEX Installation

```bash
flex --version
```

If not installed, install FLEX using your package manager or download from the official website.

### Step 2: Locate FLEX Library

Find the FLEX library location (typically `libfl.a` or `libfl.so`):

```bash
# Common locations:
ls /usr/lib/libfl*
ls /usr/lib64/libfl*
ls /usr/lib/x86_64-linux-gnu/libfl*
```

### Step 3: Set Environment Variables

Add these lines to your shell configuration file (`~/.bashrc` for bash, `~/.zshrc` for zsh, `~/.cshrc` for csh):

**For bash/zsh:**
```bash
export KPP_HOME=/path/to/kpp-box
export PATH=$PATH:$KPP_HOME/bin
```

**For csh/tcsh:**
```csh
setenv KPP_HOME /path/to/kpp-box
set path=($path $KPP_HOME/bin)
```

Apply the changes:
```bash
source ~/.bashrc  # or ~/.zshrc, ~/.cshrc, etc.
```

### Step 4: Configure Build Settings

Edit `Makefile.defs` in the KPP_HOME directory:

```bash
cd $KPP_HOME
nano Makefile.defs
```

Configure the following variables:

```makefile
# 1. C Compiler (gcc recommended)
CC=gcc

# 2. FLEX executable
FLEX=flex

# 3. FLEX library directory
FLEX_LIB_DIR=/usr/lib  # Adjust to your system

# 4. Compiler flags
CC_FLAGS=-g -Wall  # Use -O for optimization, -g for debugging

# 5. Include directory
INCLUDE_DIR=/usr/include  # Typically /usr/include or /usr/include/sys
```

### Step 5: Build KPP

```bash
cd $KPP_HOME
make
```

This will compile the KPP executable and place it in `$KPP_HOME/bin/kpp`.

### Step 6: Verify Installation

```bash
which kpp
kpp --version
```

If successful, you should see the path to the kpp executable.

## Quick Start

Here's a simple example to get you started:

### 1. Create a Mechanism Definition

Create a file named `simple.def`:

```
#include simple.spc
#include simple.eqn

#LANGUAGE Fortran90
#INTEGRATOR rosenbrock

#INITVALUES
CFACTOR = 1.0
O3 = 5.0E+11
NO = 1.0E+09

#INLINE F90_INIT
  TSTART = 0.0
  TEND = 86400.0
  DT = 3600.0
  TEMP = 298.0
#ENDINLINE
```

### 2. Create Species File

Create `simple.spc`:

```
#DEFVAR
O3  = O + O + O
NO  = N + O
NO2 = N + O + O

#DEFFIX
O2 = O + O
```

### 3. Create Equations File

Create `simple.eqn`:

```
#EQUATIONS

<R1> NO + O3 = NO2 + O2 : 1.8E-12 * EXP(-1370.0/TEMP);
```

### 4. Run KPP

```bash
kpp simple.def
```

This generates Fortran90 code with the Rosenbrock integrator.

### 5. Compile and Run

```bash
# Generated files will include simple_Main.f90, simple_Model.f90, etc.
gfortran -o simple simple_*.f90
./simple
```

## Usage Guide

### Input File Structure

A KPP mechanism consists of three main files:

#### 1. Definition File (.def)
The master file that includes other files and specifies options:

```
#include mechanism.spc     {Species definitions}
#include mechanism.eqn     {Reaction equations}

#LANGUAGE Fortran90         {Output language: Fortran77, Fortran90, C, Matlab}
#INTEGRATOR rosenbrock      {Numerical integrator}
#DRIVER general             {Driver program type}

#DOUBLE ON                  {Use double precision}
#JACOBIAN SPARSE_LU_ROW     {Sparse Jacobian with LU decomposition}

#LOOKATALL                  {Output all species to file}
#MONITOR O3; NO; NO2;       {Display these species on screen}

#CHECK O; N;                {Check mass balance for these atoms}

#INITVALUES                 {Initial concentrations}
CFACTOR = 1.0
O3 = 5.0E+11
NO = 1.0E+09
NO2 = 2.0E+08

#INLINE F90_INIT            {Inline code for initialization}
  TSTART = 0.0             {Start time (seconds)}
  TEND = 86400.0           {End time (seconds)}
  DT = 3600.0              {Time step (seconds)}
  TEMP = 298.0             {Temperature (K)}
#ENDINLINE
```

#### 2. Species File (.spc)
Defines chemical species with their atomic composition:

```
#include atoms              {Include atomic definitions}

#DEFVAR                     {Variable species (computed)}
O3  = O + O + O            {Ozone}
NO  = N + O                {Nitric oxide}
NO2 = N + O + O            {Nitrogen dioxide}
OH  = O + H                {Hydroxyl radical}

#DEFFIX                     {Fixed species (constant)}
O2  = O + O                {Molecular oxygen}
H2O = H + H + O            {Water}
M   = IGNORE               {Generic air molecule}
```

#### 3. Equations File (.eqn)
Defines chemical reactions with rate coefficients:

```
#EQUATIONS {Reaction Mechanism}

{Photolysis reactions}
<J1> NO2 + hv = NO + O    : 8.0E-3 * SUN;

{Thermal reactions}
<R1> NO + O3 = NO2 + O2   : 1.8E-12 * EXP(-1370.0/TEMP);
<R2> NO2 + O3 = NO3 + O2  : 1.4E-13 * EXP(-2470.0/TEMP);

{Three-body reactions}
<R3> O + O2 + M = O3 + M  : 6.0E-34 * (TEMP/300)^(-2.6);
```

### Running KPP

Basic syntax:

```bash
kpp [options] mechanism.def
```

Common options:
- `-h` or `--help`: Display help information
- `-v` or `--version`: Display version information

### Compiling Generated Code

#### Fortran 90:
```bash
gfortran -o model *.f90
./model
```

#### Fortran 77:
```bash
gfortran -o model *.f
./model
```

#### C:
```bash
gcc -o model *.c -lm
./model
```

#### MATLAB:
Open MATLAB and run the generated driver script.

## Configuration Options

### Language Options
```
#LANGUAGE Fortran77    {or Fortran90, C, Matlab}
```

### Integrator Options

| Integrator | Description | Best For |
|------------|-------------|----------|
| `rosenbrock` | Adaptive Rosenbrock method | Stiff systems, general use |
| `runge_kutta` | Runge-Kutta methods | Non-stiff systems |
| `sdirk` | Singly Diagonally Implicit RK | Moderately stiff systems |
| `kpp_lsode` | LSODE integrator | Large stiff systems |
| `kpp_radau5` | RADAU5 method | Very stiff systems |
| `gillespie` | Gillespie algorithm | Stochastic simulations |
| `tau_leap` | Tau-leaping method | Fast stochastic simulations |

### Driver Options
```
#DRIVER general           {General purpose driver}
#DRIVER general_adj       {Adjoint sensitivity analysis}
#DRIVER general_tlm       {Tangent linear model (forward sensitivity)}
#DRIVER general_stochastic {Stochastic simulation}
```

### Jacobian Options
```
#JACOBIAN SPARSE_LU_ROW   {Sparse, row-oriented}
#JACOBIAN SPARSE_ROW      {Sparse without LU}
#JACOBIAN FULL            {Full Jacobian matrix}
```

### Additional Options
```
#DOUBLE ON                {Double precision (recommended)}
#HESSIAN ON               {Compute Hessian matrix}
#STOICMAT ON              {Generate stoichiometric matrix}
#MEX ON                   {Generate MATLAB MEX files}
```

### Output Control
```
#LOOKATALL                {Output all species}
#LOOKAT species1; species2; {Output specific species}
#MONITOR species1; species2; {Display species during integration}
```

## Examples

KPP includes several example mechanisms in the `demo/` directory:

### Small Stratospheric Mechanism
```bash
cd $KPP_HOME/demo/examples
kpp small_f90.kpp
gfortran -o small small_*.f90
./small
```

### SAPRC-99 Mechanism
```bash
cd $KPP_HOME/demo/examples
kpp saprc_f90.kpp
gfortran -o saprc saprc_*.f90
./saprc
```

### Exploring Available Examples
```bash
ls $KPP_HOME/demo/examples/     # Example .kpp files
ls $KPP_HOME/demo/models/       # Chemical mechanism files (.def, .spc, .eqn)
```

## Advanced Features

### Sensitivity Analysis

#### Forward Sensitivity (Tangent Linear Model)
```
#DRIVER general_tlm
```
Computes sensitivities of species concentrations with respect to parameters.

#### Adjoint Sensitivity
```
#DRIVER general_adj
```
Computes sensitivities of a scalar functional with respect to initial conditions or parameters.

### Stochastic Simulation
```
#DRIVER general_stochastic
#INTEGRATOR gillespie  {or tau_leap}
```

### Custom Rate Functions

Define custom functions in inline code:
```
#INLINE F90_RCONST
  ! Custom rate constant calculation
  K_CUSTOM = 1.0E-12 * EXP(-1500.0/TEMP) * (TEMP/300.0)^2
#ENDINLINE
```

### Including External Files
```
#INLINE F90_GLOBAL
  USE external_module
#ENDINLINE
```

## Troubleshooting

### Common Issues

**1. "FLEX library not found" error**
- Verify FLEX installation: `flex --version`
- Check `FLEX_LIB_DIR` in `Makefile.defs`
- Common locations: `/usr/lib`, `/usr/lib64`, `/usr/lib/x86_64-linux-gnu`

**2. "kpp: command not found"**
- Verify `KPP_HOME` is set: `echo $KPP_HOME`
- Verify PATH includes `$KPP_HOME/bin`: `echo $PATH`
- Re-run: `source ~/.bashrc`

**3. Compilation errors in generated code**
- Ensure you're using a compatible Fortran/C compiler
- Check for proper initialization of all variables in inline code
- Verify species and reactions are correctly defined

**4. Numerical integration fails or produces errors**
- Try a different integrator (e.g., `rosenbrock` instead of `runge_kutta`)
- Adjust integration tolerances in the driver code
- Check initial conditions are physically reasonable
- Reduce the time step (DT)

**5. Mass balance errors**
- Verify atomic compositions in `.spc` file
- Check that all atoms are conserved in reactions
- Use `#CHECK` to identify problematic atoms

### Cleaning the Installation

Remove object files:
```bash
cd $KPP_HOME
make clean
```

Complete clean (removes all binaries):
```bash
cd $KPP_HOME
make distclean
```

## Documentation

- **User Manual**: `doc/kpp-UserManual.pdf` - Comprehensive documentation
- **Integrator Details**: `int/readme` - Information about numerical integrators
- **Example Mechanisms**: `demo/examples/` and `demo/models/`
- **Test Cases**: `testcase/` - Validation test cases

## License

KPP is distributed under the GNU General Public License (GPL).

Copyright (C) 1995-1997, V. Damian & A. Sandu, CGRER, University of Iowa
Copyright (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech

With contributions from:
R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany

See `gpl/` directory for full license text.

## Support

For problems, bug reports, or questions:

**Primary Contact:**
Adrian Sandu
Email: sandu@cs.vt.edu
Computer Science Department
Virginia Polytechnic Institute and State University

**Before Reporting Issues:**
1. Check this README and the user manual (`doc/kpp-UserManual.pdf`)
2. Verify your installation is correct
3. Try the included examples to ensure KPP is working
4. Include the following in your report:
   - KPP version
   - Operating system and version
   - Compiler versions
   - Complete error messages
   - Input files that reproduce the problem

## Citation

If you use KPP in your research, please cite:

```
Damian, V., Sandu, A., Damian, M., Potra, F., and Carmichael, G. R. (2002).
The Kinetic PreProcessor KPP - A Software Environment for Solving Chemical Kinetics.
Computers and Chemical Engineering, 26(11), 1567-1579.
```

---

**Quick Reference Card**

```bash
# Installation
export KPP_HOME=/path/to/kpp-box
export PATH=$PATH:$KPP_HOME/bin
cd $KPP_HOME && make

# Basic usage
kpp mechanism.def          # Generate code
gfortran -o model *.f90    # Compile (Fortran90)
./model                    # Run

# Clean
make clean                 # Remove object files
make distclean            # Remove all generated files
```
