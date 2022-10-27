import math
import numpy as np
import sys
import utilities


class ReinforcementProperties:
    def __init__(self, points, bar_diameter, mat_id):
        bar = utilities.Bar(bar_diameter)
        self.radius = bar.d / 2
        self.area = bar.area
        self.points = points
        self.mat_id = mat_id


class UnconfConcMat:
    def __init__(self, kwargs):
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
