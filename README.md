# Fluid directed rigid ball balancing using Deep Reinforcement Learning

This repository contains code for an agent that was trained using deep reinforcement learning (PPO) to balance a ball in air using fluid that is sprayed by a water hose (agent). The problem is based on this [paper](http://gamma.cs.unc.edu/DRL_FluidRigid/).

## Physics Simulator

This folder contains the code the 2D simulator which is used as the backend physics engine for the motion of rigid body and fluid.

## Autoencoder

This folder contains autoencoders which were designed during this project. This approach is used in the original paper, but was abandoned here because it was difficult to integrate with stable baselines in the available time.

## Reinforc
