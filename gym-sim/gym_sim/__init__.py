import logging
from gym.envs.registration import register

logger = logging.getLogger(__name__)

register(
    id='Sim-v0',
    entry_point='gym_sim.envs:SimEnv',
)
