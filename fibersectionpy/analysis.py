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
        max_strain = self.phi * (self.zero_strain_location - self.maxy)
        min_strain = self.phi * (self.zero_strain_location - self.miny)

        sum_moment = 0
        mat_states = []
        strains = []
        stresses = []
        for i, fiber in self.fibers.items():
            strain = self.calc_strain(fiber.xy, self.zero_strain_location)
            stress = self.materials[fiber.mat_id].stress(strain)
            strains.append(strain)
            stresses.append(stress)
            sum_moment += fiber.area * stress * (self.zero_strain_location - fiber.xy[1])
            mat_states.append(self.materials[fiber.mat_id].state)
            fiber.fail = self.materials[fiber.mat_id].fail
            if fiber.fail:
                rnd_loc = [round(i, 3) for i in fiber.xy]
                self.fail = fiber.fail + f"\n\tMat_id={fiber.mat_id}\n\tLocation={rnd_loc}"
        sum_moment += self.P * self.zero_strain_location
        self.state = State(mat_states, max_strain, min_strain, self.zero_strain_location)
        self.state.strains = strains
        self.state.stresses = stresses
        return sum_moment
