import os
import gym
import gym_sim

from stable_baselines.common.policies import MlpPolicy, FeedForwardPolicy
from stable_baselines.common.vec_env import DummyVecEnv, SubprocVecEnv
from stable_baselines import PPO2

log_dir = '/tmp/gym/'
os.makedirs(log_dir, exist_ok=True)

class CustomPolicy(FeedForwardPolicy):
	def __init__(self, *args, **kwargs):
        super(CustomPolicy, self).__init__(*args, **kwargs,
                                           layers = [128,64,64,32],
                                           feature_extraction="mlp")

env = gym.make('Sim-v0')
env = SubprocVecEnv([lambda: env for i in range(16)])

gamma = 0.95
ent_ceof = 0.1
model = PPO2(MlpPolicy, env, gamma=gamma, ent_ceof=ent_ceof, tensorboard_log=log_dir, verbose=1)
model.learn(total_timesteps=200000)
model.save("sim_agent")
