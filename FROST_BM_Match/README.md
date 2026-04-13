# Tools for track matching between WAGASCI/BabyMIND & NINJA at J-PARC (v0.0.1)

This is a set of tools for (muon) track matching between WAGASCI/BabyMIND detectors
and NINJA detectors (scintillation tracker only).

## External dependencies (same as WAGASCI+BabyMIND Geant4 simulation)

Below is a copy of README.md in the MC simulation software written by Giorgio Pintaudi.

### GEANT4

[GEANT4](http://heant4.cern/ch/) is a toolkit for the simulation of
the passage of particles through matter developed by CERN.
Geant4 v10.5.0+ is recommended.

### ROOT

[ROOT](https://root.cern.ch/) is an object-oriented program and
library developed by CERN. ROOT 6.20+ is recommended.

### BOOST

[BOOST](https://www.boost.org/) is a set of libraries for the C++
programming language that provides support for tasks and structures
such as linear algebra, pseudorandom number generation,
multithreading, image processing, regular expressions, and unit
testing. Boost 1.53+ is recommended.

### Wagasci BabyMIND Monte Carlo
[Wagasci BabyMIND Monte Carlo](https://git.t2k.org/wagasci_babymind/wagasci-babymind-monte-carlo) is a developed Monte Carlo simulation software for the WAGASCI-BabyMIND experiment.
It includes necessary libraries for the analyses.
Wagasci BabyMIND Monte Carlo v0.1.12+ is recommended.

### Wagasci BabyMIND event display
[Wagasci BabyMIND event display](https://git.t2k.org/wagasci_babymind/wagasci-babymind-event-display) is an event display tool for the WAGASCI-BabyMIND experiment.
Wagasci BabyMIND event display v0.1.0+ is recommended.

### Wagasce Reconstruction

## Installation

Only the following operative systems are tested and supported:
 - CentOS 7 (KEKCC)
 - MacOS 10.15

To install the software with all its dependencies, use a bash shell script `install.bash` in Wagasci BabyMIND Monte Carlo.
For further info about the software installation, refer to the [WAGASCI BabyMIND Monte Carlo gitlab](https://git.t2k.org/wagasci_babymind/wagasci-babymind-monte-carlo).

To compile this NINJA reconstruction software,
```shell script
mkdir build
cd build
cmake .. -DB2MC_PATH=<path/to/your/wagasci/mc/installation> \
  -DEVENT_DISPLAY_PATH=<path/to/your/wagasci/event/display/installation> \
  -DCMAKE_INSTALL_PREFIX=<path/to/your/ninja/recon/installation>
make
make install
```
also it is useful to export environment variables in .bashrc:
```shell script
export CC=$(command -v gcc) // KEKCC only
export CXX=$(command -v g++) // KEKCC only
export BOOST_ROOT=<path/to/your/boost/installation> // KEKCC only
export B2MCDIR=<path/to/your/wagasci/mc/installation>
export EVENT_DISPLAYDIR=<path/to/your/wagasci/event/display/installation>
```
then you can avoid `B2MC_PATH` and `EVENT_DISPAY_PATH` when you cmake.

### KEKCC

For further info about how to create an account and use the KEKCC
cloud computing servers, refer to the [NINJA
wiki](https://www-he.scphys.kyoto-u.ac.jp/research/Neutrino/ninja/dokuwiki/doku.php?id=software:howto:kekcc).

To recompile the package on KEKCC, move to the `build` folder and type:
```shell script
CC=$(command -v gcc) CXX=$(command -v g++) \
  BOOST_ROOT=<path/to/your/boost/installation> \
  cmake3 -DCMAKE_INSTALL_PREFIX=<path/to/your/ninja/recon/installation> \
  -DB2MC_PATH=<path/to/your/wagasci/mc/installation> \
  -EVENT_DISPLAY_PATH=<path/to/your/wagasci/event/display/installation> ..
```

## Usage

First just type:
```shell script
<path/to/your/program/installation>
```
then you can find discriptions of each variable.

### File Separator

This program is used for separation of the NINJA tracker raw data into a day unit file.
Since WAGASCI-BabyMIND files are provided in a day unit, the tracker file also should be
separated into the unit.

### Hit Converter

This program is used for NINJA tracker raw data conversion to WAGASCI-BabyMIND general data format.
It converts NINJA tracker EASIROC raw data into WAGASCI-BabyMIND data format, B2HitSummary, and
push back each hits into the file.

#### Note: This is only used in real data because simulated data is generated in B2 data format.

### Track Match

This program is used for track matching between NINJA tracker and WAGASCI-BabyMIND detectors.
The merged file created in Hit Converter processes are analyzed and B2TrackSummary with NINJA
tracker hit is created for each spill.