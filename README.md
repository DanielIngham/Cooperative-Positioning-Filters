# Cooperative Positioning Filters
Repository containing the implementation of Bayesian filters used for multi-robot cooperative positioning, including:
- Extended Kalman Filter (EKF),
- Iterative Extended Kalman Filter (IEKF), and
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
To run a one of the aforementioned filters, the user can select the filter by setting the `Filter` environment variable before executable:
```bash
make run Filter=<Filter Name>
```
`<Filter Name>` can be either: `EKF`, `IEKF`, or `INFO`. 

To select a dataset, the user can set the `Dataset` environment variable before execution:
```bash
make run Filter=<Filter Name> Dataset=<Dataset Number>
```
`Dataset` can be any value in the range [1,9]. Note that dataset 9 has incorrect groundtruth values, and therefore it is not recomended to use it.

# Build Options
## Coupled vs Decoupled
The measurement update can be performed either in a **coupled** or **decoupled** manner. The coupled approach creates an augmented state vector and covariance matrix containing the state and covariance of both the ego vehicle and the agent observed respectively. The decoupled approach adds the covariance of the agent observed to the ego robots measurement noise. Since the state of the agent in marginalized out in the coupled approach, both approaches result in identical inference. However, the decoupled approach is less computationally expensive. 

The decoupled is the default, however, the user can compile the project using the coupled approach during linking by adding the `COUPLED`:
```bash
make COUPLED=1
```
## Normal vs Robust 
The measurement update can be performed using the *Huber* loss function for both the EKF and IEKF, which results in improved robustness against outliers. To use the Huber loss function, the user must compile the project using the `ROBUST` flag:
```bash
make ROBUST=1
```
