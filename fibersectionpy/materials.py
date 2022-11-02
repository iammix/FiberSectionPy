import math
import numpy as np
import sys
import utilities

def conf_pressure_circle(fyh, bar_number, d, s):
    bar = Stand


class ReinforcementProperties:
    def __init__(self, points, bar_diameter, mat_id):
        bar = utilities.Bar(bar_diameter)
        self.radius = bar.d / 2
        self.area = bar.area
        self.points = points
        self.mat_id = mat_id


class UnconfConcMat:
    def __init__(self, **kwargs):
        self.fail = False
        self.state = 'White'
        self.fpc = kwargs['fpc']
        ratio = 1000
        if 'Ec' in kwargs:
            self.Ec = kwargs['Ec']
        else:
            self.Ec = math.sqrt(self.fpc * ratio) * 57
        self.fcr = 4 * math.sqrt(self.fpc * ratio) / ratio
        m = 1 + 3600 / (self.fpc * ratio) / ratio
        if 'ecp' in kwargs:
            self.ecp = kwargs['ecp']
        else:
            self.ecp = -self.fpc * m / self.Ec
        self.r_unconf = 2

        if 'ecu' in kwargs:
            self.ecu = kwargs['ecu']
        else:
            self.ecu = -0.005
        if 'tension' in kwargs:
            self.tension = kwargs['tension']
        else:
            self.tension = False
        self.useful_points = [self.ecu, self.ecp * 2, self.ecp, 0]

    def stress(self, ec):
        if ec > 0:
            if self.tension:
                return ec
            else:
                self.state = 'White'
                return 0
        else:
            if ec >= self.ecp * 2:
                r = self.r_unconf
                ecp = self.ecp
                fc = -r * (ec / ecp) / (r - 1 + (ec / ecp) ** r) * self.fpc
                if ec > self.ecp:
                    self.state = 3
                else:
                    self.state = 1
                return fc

            elif ec >= self.ecu:
                slope = self.stress(2 * self.ecp) / (2 * self.ecp - self.ecu)
                b = -slope * self.ecu
                self.state = 7
                return slope * ec + b
            else:
                self.state = 8
                return 0

    def tension(self, ec):
        elastic = ec * self.Ec
        if elastic > self.fcr:
            if ec <= 0.002:
                return 0.7 * self.fcr(1 + math.sqrt(500 * ec))
            else:
                return 0
        else:
            return elastic
        self.state = 8


class SteelBilinMaterial:
    def __init__(self, E, fy, fsu, e_sh, e_su, P):
        self.Es = E
        self.fy = fy
        self.fsu = fsu
        self.e_y = fy / E
        self.e_sh = e_sh
        self.e_su = e_su
        self.P = P
        self.state = 'Black'
        self.fail = False

        self.plot_points = [-self.e_su, -self.e_sh, -self.e_y, 0, self.e_y, self.e_sh, self.e_su]

    def stress(self, e):
        if e < 0:  # Compression
            if e > -self.e_y:
                self.state = 'Yellow'
                return self.Es * e
            elif e > -self.e_sh:
                self.state = 'Pink'
                return -self.fy
            elif e > -self.e_su:
                self.state = 'Red'
                return -(self.fsu - (self.fsu - self.fy) * ((self.e_su + e) / (self.e_su - self.e_sh)) ** self.P)
            else:
                self.state = 'Black'
                self.fail = f'Steel Fracture\nMax Strain = {abs(round(self.e_su, 5))}\n'
                return 0
        else:  # Tension
            if e < self.e_y:
                self.state = 'Yellow'
                return self.Es * e
            elif e < self.e_sh:
                self.state = 'Pink'
                return self.fy
            elif e < self.e_su:
                self.state = 'Red'
                return self.fsu - (self.fsu - self.fy) * ((self.e_su - e) / (self.e_su - self.e_sh)) ** self.P
            else:
                self.state = 'Black'
                self.fail = f'Steel Fracture\n Max Strain = {abs(round(self.e_su, 5))}'
                return 0