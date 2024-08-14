# Analysis and Plotting Scripts for SVJl Signal

**Abhishek Dey**

This repository documents the work completed during my ATLAS Summer Internship, where I focused on analyzing the Semi-Visible Leptonic Jet (SVJl) signal—a hypothetical dark valley signature [1]. The repository includes scripts and code used to replicate the analysis and generate the relevant plots. Additionally, it provides the `bash` scripts necessary to process background data.

## Table of Contents

- [Analysis and Plotting Scripts for SVJl Signal](#analysis-and-plotting-scripts-for-svjl-signal)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Files in the Repository](#files-in-the-repository)
  - [Sources](#sources)

## Requirements

To use the scripts in this repository, the following tools are required:

- **ROOT**: Essential for processing and analyzing ROOT files. Instructions for installation can be found at [ROOT CERN](https://root.cern/install/).
- **FastJet**: Needed for jet clustering algorithms. Installation details are available at [FastJet](http://fastjet.fr/).
- **bash**: For executing shell scripts.
- It is recommended to use an lxplus environment to ensure compatibility.

## Installation

To install, ensure that ROOT and FastJet are properly set up in your environment. After cloning this repository to your desired directory, make the provided shell scripts executable by running the following commands:

```bash
chmod +x calculate_event_weights.sh
chmod +x run_analyses_using_event_weights.sh
```

## Usage

The analysis workflow begins with calculating event weights using the `calculate_event_weights.sh` script. This script processes ROOT files containing event data, computing weights based on factors such as cross-section, k-factor, filter efficiency, and luminosity. The computed weights are stored in an `event_weights.txt` file. Before running this script, ensure the `base_path` variable within the script is set to the correct directory containing your ROOT files. Here is an example of the tree that is used in this analysis:

```
├── 700323
│   ├── 40385045._000003.output.root
│   ├── 40385045._000004.output.root
│   ├── 40385045._000005.output.root
│   ├── 40385045._000006.output.root
│   ├── 40385045._000007.output.root
│   ├── 40385045._000008.output.root
│   ├── 40385045._000009.output.root
│   ├── 40385045._000010.output.root
│   └── 40385045._000011.output.root
├── 700324
│   ├── 40385047._000001.output.root
│   ├── 40385047._000002.output.root
│   ├── 40385047._000003.output.root
│   ├── 40385047._000004.output.root
│   ├── 40385047._000005.output.root
│   ├── 40385047._000006.output.root
│   ├── 40385047._000007.output.root
│   ├── 40385047._000008.output.root
│   ├── 40385047._000009.output.root
│   ├── 40385047._000010.output.root
│   ├── 40385047._000011.output.root
│   ├── 40385047._000012.output.root
│   └── 40385047._000013.output.root
└── 700325
    ├── 40385048._000001.output.root
    ├── 40385048._000002.output.root
    ├── 40385048._000003.output.root
    ├── 40385048._000004.output.root
    ├── 40385048._000005.output.root
    ├── 40385048._000006.output.root
    ├── 40385048._000007.output.root
    ├── 40385048._000008.output.root
    ├── 40385048._000009.output.root
    ├── 40385048._000010.output.root
    └── 40385048._000011.output.root
```

Each folder in the directory structure is named according to its corresponding DSID, which aligns with the design of the analysis scripts. The `calculate_event_weights.sh` script is tailored to work with this directory structure, automating the event weight calculation process for each DSID and the individual ROOT files contained within them.

The analysis is conducted using a single signal ROOT file, which is processed through the `Analysis_script.C` file. This script is specifically designed to handle multiple DSIDs and ROOT files. The automation provided by `run_analyses_using_event_weights.sh` ensures that the analysis is consistently applied across all datasets.

Upon completing the analysis, plots can be generated using the `plotting_macro.C` script. This macro produces ATLAS-style plots, formatted on a log-y scale. Ensure that all the background DSIDs are merged into a single `.root` file using `hadd` prior to running the plotting script. An example `hadd` command is given below:

```
cd path/to/output_for_this_DSID/
hadd -f all_background.root *.root
```

## Files in the Repository

- **calculate_event_weights.sh**: A script to calculate event weights for various datasets.
- **run_analyses_using_event_weights.sh**: A script to perform analysis using the calculated event weights.
- **Analysis_script.C**: A C++ script used for analyzing the SVJl signal.
- **plotting_macro.C**: A macro script designed to generate ATLAS-style plots.

## Sources

[1]: C. Cazzaniga and A. de Cosa, “Leptons lurking in semi-visible jets at the LHC,” The European Physical Journal C, 82, 10.1140/epjc/s10052-022-10775-2 (2022).

[2]: M. Cacciari, G. P. Salam, and G. Soyez, “Fastjet user manual”, European Physical Journal C 72, 1896 (2012)