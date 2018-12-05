import gym
import gym_sim
import numpy as np

from stable_baselines import PPO2

env = gym.make('Sim-v0')

model = PPO2.load('sim_agent')

obs = env.reset()

while(True):
	for i in range(1000):
		action, _states = model.predict(obs)
		obs, rewards, dones, info = env.step(action)
		if dones:
			break

	print('Episode lasted for ',i,' timesteps!')
	
	# because the agent is not robust enough, we force it to keep playing until a successful episode results
	if(i>=950)
		break
