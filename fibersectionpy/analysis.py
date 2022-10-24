import sys


class State:
    def __init__(self, mat_state, max_strain, min_strain, y_loc):
        self.mat_state = mat_state
        self.max_strain = max_strain
        self.min_strain = min_strain
        self.y_loc = y_loc
        self.strains = []
        self.stresses = []


class Fiber:
    def __init__(self, area, coords, mat_id):
        self.area = area
        self.xy = coords
        self.mat_id = mat_id
        self.fail = False

class FiberModel:
    def __init__(self):
        self.fibers = {}
        self.materials = {}
        self.zero_strain_location = 0
        self.phi = 0
        self.P = 0
        self.fail = False
        self.maxy = 0
        self.miny = 0
        self.state = None

    def calc_strain(self, coords, y_intercept):
        return self.phi * (y_intercept - coords[1])
    def force_balance(self, y_intercept):
        self.zero_strain_location = y_intercept
        sum_force = 0
        for i, fiber in self.fibers.items():
            strain = self.calc_strain(fiber.xy, y_intercept)
            stress = self.materials[fiber.mat_id].stress(strain)
            sum_force += fiber.area * stress
        return sum_force + self.P

    def calc_moment(self):
        