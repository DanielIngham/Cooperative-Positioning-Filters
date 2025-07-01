# Cooperative Positioning Filters
Repository containing the implementation of Bayesian filters used for multi-robot cooperative positioning, including:
- Extended Kalman Filter (EKF),
- Iterative Extended Kalman Filter (IEKF), and
- Extended Information Filter.

# Submodule
This repository contains the [Cooperative Positioning Data Handler](https://github.com/DanielIngham/Cooperative-Positioning-Data-Handler), as a submodule. Therefore, to clone the repository, the following commands are required:
```bash
git clone https://github.com/DanielIngham/Cooperative-Positioning-Filters.git
cd Cooperative-Positioning-Filters
git submodule update --init --recursive
```

# Measurement Update Options
## Coupled vs Decoupled
The measurement update can be performed either in a **coupled** or **decoupled** manner. The coupled approach creates an augmented state vector and covariance matrix containing the state and covariance of both the ego vehicle and the agent observed respectively. The decoupled approach adds the covariance of the agent observed to the ego robots measurement noise. Since the state of the agent in marginalized out in the coupled approach, both approaches result in identical inference. However, the decoupled approach is less computationally expensive. 

The decoupled is the default, however, the user can compile the the project using the coupled approach during linking by adding:
```
make COUPLED=1
```
## Normal vs Robust 
The measurment update can be performed using the *Huber* loss function for both the EKF and IEKF, which results in improved robustness against outliers. To use the Huber loss function, the user must complie the project using the compile the `ROBUST` flag:
```
make ROBUST=1
```
