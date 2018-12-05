
## 2D Coupled Physics Simulator
The environment consists of the 2d physics simulator available in this repository

## Installation

Install the [OpenAI gym](https://gym.openai.com/docs/).

Then install this package via

```
pip install -e .
```

## Usage

```
import gym
import gym_sim

env = gym.make('Sim-v0')
```
## The Environment

### Description

The environment simulates the coupled motion of a 2D rigid ball and water using the simulator code available in this repository. A water hose (controller) is the agent which is trained using deep reinforcement learning. The goal is to keep the ball in air was as long as possible. The domain is of size 1 m x 1.4 m and is represented by [0.3, 1.3] x [0.1, 1.5], with origin at (0, 0). The MAC grid is used as a scratch pad for computations and is given by [0, 1.6] x [0, 1.6].

### Actions

Type: Box(3)

| Num  | Action               | Min | Max |
| --- | --- | --- | --- |
| 0    | Acceleration         | -1  | 1   |
| 1    | Angular Acceleration | -1  | 1   |
| 2    | Amount of fluid shot |  0  | 1   |

Note: acceleration and angular acceleration are scaled according in the C++ code directly. This is the action space for the controller (water hose) which shoots fluid inside the domain.

### Observation

Type: Box(401)

| Num  | Action               | Min | Max |
| --- | --- | --- | --- |
| 0    | Ball x-position        | -Inf  | Inf   |
| 1    | Ball y-position | -Inf  | Inf   |
| 2    | Ball angular vel. |  -Inf  | Inf   |
| 3    | Ball x-velocity |  -Inf  | Inf   |
| 4    | Ball y-velocity |  -Inf  | Inf   |
| 5    | Controller x-position |  0.4  | 1.2   |
| 6    | Controller x-velocity |  -Inf  | Inf   |
| 7    | Controller angle |  -60<sup>o</sup>  | 60<sup>o</sup>   |
| 8    | Controller angular vel. |  -Inf  | Inf   |
| 9-400    | Fluid velocity field |  -Inf  | Inf   |

Note: the controller is constrained to move along only x-direction. The values for controller's x-position and angular velocity are clipped to ensure it doesn't go outside the domain boundaries. The fluids velocity field is extracted using a traveling window which is centered at the center of the ball. The size of the window is 14 x 14 and samples the x and y velocities.

### Reward

The following reward function is used,

<img src="https://latex.codecogs.com/gif.latex?r&space;=&space;w_{pos}&space;e^{-\vert&space;\vert&space;x_{t}&space;-&space;x_{start}&space;\vert&space;\vert}&space;&plus;&space;w_{vel}&space;e^{-\vert&space;\vert&space;v_{t}&space;\vert&space;\vert}&space;&plus;&space;w_{shoot}&space;e^{-\vert&space;\vert&space;\delta_{shoot}&space;\vert&space;\vert}" title="r = w_{pos} e^{-\vert \vert x_{t} - x_{start} \vert \vert} + w_{vel} e^{-\vert \vert v_{t} \vert \vert} + w_{shoot} e^{-\vert \vert \delta_{shoot} \vert \vert}" />

Note: the weights are hyperparameters which are to be tuned. For this project, (5, 2.5, 5) were used. <img src="https://latex.codecogs.com/gif.latex?\delta_{shoot}" title="\delta_{shoot}" /> represents the amount of fluid shot by the controller.

### Starting State

The observations from 2-8 are assigned uniform random values between (-0.01, 0.01). The observation 0 is assigned uniform random values between (0.8-0.01, 0.08+0.01) and observation 1 is assigned uniform random values between (1.0-0.01, 1.0+0.01). This is because the starting position of the ball is kept at (0.8, 1.0). The fluids velocity field is assigned a 0 value as its state depends on the controller and I choose to start with a blank slate everytime.

### Episode Termination

The episode ends if the y-position of the ball falls below 0.45 m.
