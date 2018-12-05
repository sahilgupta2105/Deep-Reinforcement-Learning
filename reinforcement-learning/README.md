# Deep Reinforcement Learning: PPO

Uses the [stable baselines](https://github.com/hill-a/stable-baselines) library for training the agent using PPO algorithm. The code implements a multi-processing vectorized environment that drastically speeds up the learning process.

A custom policy network is used for the learning task which is a fully connected neural network with the following architecture: 128 x 64 x 64 x 32.
