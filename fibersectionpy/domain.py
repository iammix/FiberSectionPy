from fibersectionpy.analysis import FiberModelDB, Fiber, State
from fibersectionpy.materials import *
import numpy as np
import math
import scipy
from meshpy.geometry import GeometryBuilder, make_circle
import matplotlib.pyplot as plt
from matplotlib.patches import CirclePolygon, Ellipse
from matplotlib.widgets import Slider
import matplotlib.animation as anim
import sys


class FiberModel:
    def __init__(self, figsize=6):
        self.materials = {}
        self.builder = GeometryBuilder()
        self.reinforcement = []
        self.outer_facets = []
        self.mesh_ares = {}
        self.mesh_centroids = {}
        self.mesh = None
        self.ele_mat_primitive = []
        self.ele_mat = {}
        self.load_angle = 0
        self.phi_list = []
        self.M_list = []
        self.state_id = 0
        self.states = {}
        self.state_y_loc = {}
        self.analysis_model = FiberModelDB()
        self.fail = False
        self.figsize = figsize

    @staticmethod
    def conf_pressure(shape, **kwargs):
        if shape == 'circle':
            return conf_pressure_circle(kwargs['fyh'], kwargs['bar'], kwargs['s'], kwargs['D'])
        elif shape == 'rect':
            return conf_pressure_rect(kwargs['fyh'], kwargs['bar'], kwargs['s'], kwargs['b'], kwargs['w'], kwargs['nx'],
                                      kwargs['ny'])
        else:
            raise ValueError(
                "Calculation for confining concrete can be applied only on rectangular and circular secitons")

    @staticmethod
    def color_from_state(state):
        colors = ["Gray", "Purple", "Pink", "Blue", "Black", "Green", "Yellow", "Red", "White"]
        if type(state) == int:
            return colors[state]
        else:
            return state

    def add_material(self, mat_num, mat, **kwargs):
        if mat == 'concrete':
            if 'fple' in kwargs:
                self.materials[mat_num] = Conf