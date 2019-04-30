import math
import numpy as np
import gym
from gym import spaces
from gym.utils import seeding

class RKEnv(gym.Env):
    metadata = {'render.modes': ['air']}

    def __init__(self, min_action, max_action,
    min_position, max_position):
        super(RKEnv, self).__init__()
        self.min_action = min_action
        self.max_action = max_action
        self.min_state = min_position
        self.max_state = max_position
        self.action_spaces= np.array([self.min_action])
        self.observation_space = np.array([self.max_action])
 
        self.seed()


    def step(self, a):
        reward = 0.0
        action = self.action_spaces[a]
        reward += self.ale.act(action)
        ob = self.ale.get_obs()
        return ob, reward, self.ale.game_over(), {"ale.lives": self.ale.lives()}

    
    def reset(self):
        self.state = np.array([self.np_random.uniform(low=-0.6, high=-0.4), 0])
        return np.array(self.state)
    
    def render(self):
        return np.array(self.state)
        
    def seed(self, seed=None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]
        
    def close(self):
        return 0