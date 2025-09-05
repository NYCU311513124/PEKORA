# Partially Estimated Kernel Optimization and Random Adaptive (PEKORA) Algorithm Applied to Microwave Circuit Design 
This program implements the Partially Estimated Kernel Optimization (PEKO) algorithm for optimizing microwave circuit designs. By combining the particle filter method with transmission line theory, it performs parameter optimization to meet requirements for reflection coefficient and phase response. It can be applied to CRLH phase shifters and impedance matching structures.

## Core Algorithm Overview
1. PEKO Algorithm (Partially Estimated Kernel Optimization) : The goal of this algorithm is to transform the load impedance in a circuit into the desired input impedance, thereby achieving impedance matching.
2. Particle Filter : To overcome the computational complexity of grid search method, this program adopts the particle filter approach. It uses a finite number of randomly sampled points (particles) to estimate results, and by repeatedly resampling and updating particle weights, it efficiently identifies the optimal set of design parameters.

## Supported Circuit Structures
This MATLAB implementation of the PEKORA algorithm can be used to design :
1. CRLH Phase Shifter : Achieves reflection coefficient requirements and satisfies specific phase variation conditions.
2. Impedance Matching Structure : Achieves the specified reflection coefficient requirements.
3. Switchable Mode : Can switch between impedance matching structure and CRLH phase shifter.
4. Result Prioritization : The program ranks results according to user-defined priorities (reflection coefficient, phase variation, or a weighted combination of both).

## Adjustable Parameters
The program allows tuning of multiple parameters to control the calculation process :
1. Basic Specifications : Frequency range, input and load impedances, target reflection coefficient, and center frequency.
2. Structural Parameters : Characteristic impedance and electrical length ranges of transmission lines and stubs.
3. Computation Controls : Circuit order, number of particles, number of recalculations, particle distribution type (uniform or normal), and resampling method.  
Based on these parameters, the program performs iterative calculations using the particle filter and outputs the optimized structural parameters that satisfy the design objectives.

## Usage
1. Open MATLAB or a compatible environment.
2. Load the PEKORA.m file.
3. Set the calculation parameters.
4. Run the program to output the optimal parameter set and the corresponding response results.

## Potential Applications
1. Microwave circuit and antenna matching network design
2. Phase compensation or phase shifter design
3. Impedance optimization in RF and communication systems
