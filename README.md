<a id="readme-top"></a>

<!-- ABOUT THE PROJECT -->
## About The Project

### Project Introduction
This project aims to conduct numerical simulations to study the development mechanisms of ground fissures and their potential impact on urban infrastructure. We employ a numerical model based on Peridynamic (PD) theory to simulate ground fissures in the Songzhuang area of Beijing, Xi'an area, and plates with pre-existing cracks, predicting and analyzing the development of ground fissures.

### Technical Approach

* Model Construction: Build a numerical model of geological strata to simulate the behavior of ground fissures under different geological conditions.
* Data-Driven: Calibrate and validate the model using actual drilling data and geological surveys.
* Multi-scenario Analysis: Simulate crack propagation under various geological and environmental conditions to assess risks in different scenarios.

### Key Features

* High-Precision Simulation: The advanced PD theory can accurately capture stress distribution at crack tips and crack propagation paths.
* User-Friendly: Provide a simple and clear user interface and instructions for easy operation and understanding by non-professionals.
* Visualization Output: Generate intuitive displacement and damage distribution diagrams to help users understand simulation results.


<!-- GETTING STARTED -->
## Getting Started


### Prerequisites

* MATLAB 2018a or later.

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/1322403796/Songzhuang_Ground_Fissure.git
   ```
2. Open MATLAB and navigate to the cloned repository directory.

<!-- USAGE EXAMPLES -->
## Usage

Run one of the following three main scripts to start the corresponding simulation:

* plate_with_a_preexisting_crack_v50 and plate_with_a_preexisting_crack_v20 : Simulate plates with pre-existing cracks under varying loads
* XIAN_simulate.m : Simulate ground fissures in the Xi'an area.
* SONGZHUANG_simulation.m : Simulate ground fissures in the Songzhuang area.

Each script will perform the following steps:

* Generate a grid of material points representing geological layers.
* Apply initial conditions and boundary values.
* Execute PD simulation.
* Output displacement and damage fields.

## Visualization

Each script includes commands to generate visualizations of displacement and damage fields.

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<!-- CONTACT -->
## Contact

Lingfei Wang - [@email]wlfzzzzzz@163.com

## Figures and Outputs

The simulation generates several figures and output datasets, including:
* Displacement fields at various time steps.
* Damage fields indicating the extent of fissure propagation.
* Line plots of horizontal and vertical displacements at key points.

<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Gratitude is expressed to the Hebei Cangzhou Groundwater and Land Subsidence National Observation and Research Station, the National Natural Science Foundation of China, and the R&D Program of Beijing Municipal Education Commission for their support.
