"""
Simulate the 2D ball and fluid simulation.
Episode ends if ball deviates too much from start position
"""

# core modules
from gym_sim.envs import physics_sim
import math
import random

from gym import spaces
import gym
import numpy as np

class SimEnv(gym.Env):
    """
    Define an environment using Physics simulator.

    The environment defines which actions can be taken at which point and
    when the agent receives which reward.
    """

    def __init__(self):
        self.__version__ = "0.1.0"
        
        self.NX = 30
        self.NY = 30
        self.H = 1.6 / self.NX
        self.TOTAL_TIME_STEPS = 1000
        self.time_step = 0.01
        self.ball_x_start = 0.8
        self.ball_y_start = 1.0
        self.epsilon = 0.55# used to terminate the episode

        self.sim = physics_sim.Simulator(self.NX, self.NY, self.H, 0.2*self.H)

        a_low = np.array([-1,-1,0])
        a_high = np.array([1,1,1])
        self.action_space = spaces.Box(a_low, a_high, dtype=np.float32)

        self.DOF = 401
        low = np.array([-np.inf]*self.DOF)
        high = np.array([np.inf]*self.DOF)
        self.observation_space = spaces.Box(low, high, dtype=np.float32)

        self.state = None
        self.done = False
        self.current_time = 1

    def step(self, action):
        """
        The agent takes a step in the environment.

        Parameters
        ----------
        action : (acc, ang. acc, shot/not)

        Returns
        -------
        ob, reward, episode_over, info : tuple
            ob (object) :
                an environment-specific object representing your observation of
                the environment.
            reward (float) :
                amount of reward achieved by the previous action. The scale
                varies between environments, but the goal is always to increase
                your total reward.
            episode_over (bool) :
                whether it's time to reset the environment again. Most (but not
                all) tasks are divided up into well-defined episodes, and done
                being True indicates the episode has terminated. (For example,
                perhaps the pole tipped too far, or you lost your last life.)
            info (dict) :
                 diagnostic information useful for debugging. It can sometimes
                 be useful for learning (for example, it might contain the raw
                 probabilities behind the environment's last state change).
                 However, official evaluations of your agent are not allowed to
                 use this for learning.
        """
        
        if self.done:
            raise RuntimeError("Episode is done")
        new_obs = list(self.sim.advance(self.time_step, action.astype(float) ))
        reward = self._get_reward(action, new_obs)
        self.state = new_obs
        self._check_ball_status(new_obs)
        self.current_time=self.current_time+1
        if self.current_time > self.TOTAL_TIME_STEPS:
            self.done = True
        return np.array(self.state), reward, self.done, {}

    def _check_ball_status(self, l):
        delta_y = self.ball_y_start - l[4]

        if delta_y >= self.epsilon:
            self.done = True

    def _get_reward(self, a, o):
        """Reward is calculated by observing rigid body."""
        w_pos = 5
        w_vel = 2.5
        w_del = 5
        re = w_pos*(np.exp( -((self.ball_x_start-o[3])**2 + (self.ball_y_start-o[4])**2) )) + \
             w_vel*( np.exp( -(o[0]**2+o[1]**2) ) ) + w_del*( np.exp( -a[2] ) )
        return re

    def reset(self):
        """
        Reset the state of the environment and returns an initial observation.
        NOTE: for fluid initial state is always set at zero because it is dependent
        on the controller's fire
        Returns
        -------
        observation (object): the initial observation of the space.
        """
        s1 = self.sim.reset_simulator();
        s2 = np.zeros((392,))
        self.done = False
        self.current_time = 1
        self.state = np.concatenate((np.array(s1),s2))
        return self.state

    def render(self, mode='human'):
    	raise NotImplementedError

    def seed(self, seed):
        random.seed(seed)
        np.random.seed

    def close(self):
        print('Closing env')
