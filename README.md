# Cooperative Positioning Filters
Repository containing the implementation of Bayesian filters used for multi-robot cooperative positioning, including:
- Extended Kalman Filter (EKF),
- Iterative Extended Kalman Filter (IEKF),
- Converted Measurement Extended Kalman Filter, and
- Extended Information Filter.
## Documentation
The documentation for this implementation of the cooperative positioning filters can be found [here](https://danielingham.github.io/Cooperative-Positioning-Filters/).
## Data Handler Submodule
This repository contains the [Cooperative Positioning Data Handler](https://github.com/DanielIngham/Cooperative-Positioning-Data-Handler), as a submodule. Therefore, to clone the repository, the following commands are required:
```bash
git clone https://github.com/DanielIngham/Cooperative-Positioning-Filters.git
cd Cooperative-Positioning-Filters
git submodule update --init --recursive
```
# Running a Filter
To run a one of the aforementioned filters, the user should first build the project using CMake:
```bash
mkdir build && cd build
cmake ..
cmake --build .
```
This compiles a target for each of the filters:
- `./EKF` 
- `./IEKF` 
- `./CMEKF` 
- `./INFO` 

## Selecting a Dataset
To select a dataset, the [DataHandlers](https://github.com/DanielIngham/Cooperative-Positioning-Data-Handler) [ArgumentHandler.h](https://github.com/DanielIngham/Cooperative-Positioning-Data-Handler/blob/master/ArgumentHandler.h) is used. For example, to set the dataset for the `EKF`: 
```bash
./EKF -d <Dataset>
```
`Dataset` can be any value in the range [1,9]. Note that dataset 9 has incorrect groundtruth values, and therefore it is not recommended using it.

# Build Options
## Coupled vs Decoupled
The measurement update can be performed either in a **coupled** or **decoupled** manner. The coupled approach creates an augmented state vector and covariance matrix containing the state and covariance of both the ego vehicle and the agent observed respectively. The decoupled approach adds the covariance of the agent observed to the ego robots measurement noise. Since the state of the agent in marginalized out in the coupled approach, both approaches result in identical inference. However, the decoupled approach is less computationally expensive. 

The decoupled is the default, however, the user can compile the project using the coupled approach during linking by adding the `COUPLED`:
```bash
cmake -DCOUPLED=ON ..
```
## Normal vs Robust 
The measurement update can be performed using the *Huber* loss function for both the EKF and IEKF, which results in improved robustness against outliers. To use the Huber loss function, the user must compile the project using the `ROBUST` flag:
```bash
cmake -DROBUST=ON ..
```
