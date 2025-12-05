#!/usr/bin/env python
'''PyQt based reciprocal space viewer for the psi diffractometer at 8-ID of APS.

TODO: 
* (Done) target HKL visualization
* visualization of Eward's sphere of area detector
* scan motor coverage

@Author: Hao Zheng
@Email: hao.zheng@anl.gov
@Data: 10/7/2025

@modified psiviewer.py to use hklpy2 package
@modified by Sam Marks 12/04/2025
'''


import os
import sys
import numpy as np
from functools import partial
from time import sleep
from PyQt5.QtWidgets import (
    QMainWindow, 
    QApplication, 
    QLabel, 
    QWidget, 
    QGridLayout, 
    QLineEdit, 
    QPushButton, 
    QVBoxLayout, 
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QComboBox,
)
from PyQt5.QtCore import (
    QTimer,
    Qt,
)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from mpl_toolkits.mplot3d import Axes3D  # Needed to enable 3D plotting
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
# from hkl.geometries import SimulatedE6C
# from hkl.diffract import Constraint
from epics import caget, caput

import hklpy2
sim6c = hklpy2.creator(geometry="E6C", name="sixc")
from hklpy2.blocks.sample import Lattice


DEFAULT_PIX = 94
DEFAULT_ALIGN_TIME = 10000
DEFAULT_PV_TIME = 500 #ms
DEFAULT_PLOT_TIME = 100 #ms
WINDOW_WIDTH = 650
WINDOW_HEIGHT = 1050
DISPLAY_HEIGHT = 35
BUTTON_WIDTH = 70
BUTTON_HEIGHT = 35
ERROR_MSG = 'ERROR'
No_COMM_MSG = 'PV OFFLINE'
SYNCED_MSG = 'Synced'
NOT_SYNCED_MSG = 'Not Synced'

###### MOTOR PVS #######
DELTA_RBV = ' 8ideSoft:CR8-E1:m5.RBV'
NU_RBV = '8ideSoft:CR8-E1:m4.RBV'
MU_RBV = '8ideSoft:CR8-E1:m6.RBV'
ETA_RBV = '8ideSoft:CR8-E1:m7.RBV'
CHI_RBV = '8ideSoft:CR8-E1:m8.RBV'
PHI_RBV = '8ideSoft:CR8-E1:m9.RBV'

###### XRAY PVS #######
ENERGY_RBV = 'S08ID:USID:EnergyM' # this is the upstream gap energy, maybe update to mono
# POLARIZATION_RBV = 'S29ID:ActualModeM'

###### UB PVS ######
UB_MATRIX_RBV = ''
UB_ENERGY_RBV = ''
UB_OR0_RBV = ''
UB_OR1_RBV = ''
LATTICE_RBV = ''
SAMPLE_NAME_RBV = ''
NO_SOLUTION_MSG = 'No solution!'

POLARIZATION_DICT = {
    0: 'CW, RCP',
    1: 'CCW, LCP',
    2: 'V',
    3: 'H',
    4: 'H, Neg',
}

DEFAULT_SAMPLE = {
    'name': 'SrTiO3',
    'lattice': [3.905, 3.905, 3.905, 90, 90, 90],
    'UB_energy': 8000,
    'or0':[0, 0, 1, 98.213, 0, 0, 49.1065, 0, 0],
    'or1':[1, 0, 0, 98.213, 0, 0, 49.1065+90, 0, 0],
    'UB': np.array([[1.609, 0, 0],[0, 1.609, 0],[0, 0, -1.609]]),
}

DEFAULT_GEOMETRY = {
    'delta': 0,
    'nu': 0,
    'mu': 0,
    'eta': 0,
    'chi': 0,
    'phi': 0,
    'energy': 8000,
    'polarization': 0
}

PSI_MODES = [
    'bissector_vertical', 
    'constant_omega_vertical', 
    'constant_chi_vertical', 
    'constant_phi_vertical', 
    'lifting_detector_phi', 
    'lifting_detector_omega', 
    'lifting_detector_mu', 
    'double_diffraction_vertical', 
    'bissector_horizontal', 
    'double_diffraction_horizontal', 
    'psi_constant_vertical', 
    'psi_constant_horizontal', 
    'constant_mu_horizontal'
]

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 
    
    def update(self, x, y, z, dx, dy, dz):
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
    return arrow

setattr(Axes3D, 'arrow3D', _arrow3D)


def rot(v, axis, deg):
    '''
    Returns matrix of rotation about x-, y-, or z-axis
    '''
    rad = np.deg2rad(deg)
    if axis.lower()=='x':
        r = np.array([[1, 0, 0],[0, np.cos(rad), -np.sin(rad)],[0, np.sin(rad), np.cos(rad)]]) 
    elif axis.lower()=='y':
        r = np.array([[np.cos(rad), 0, np.sin(rad)],[0, 1, 0],[-np.sin(rad), 0, np.cos(rad)]]) 
    elif axis.lower()=='z':
        r = np.array([[np.cos(rad), -np.sin(rad), 0],[np.sin(rad), np.cos(rad), 0],[0, 0, 1]])
    elif axis.lower()=='-x':
        r = np.array([[1, 0, 0],[0, np.cos(rad), np.sin(rad)],[0, -np.sin(rad), np.cos(rad)]]) 
    elif axis.lower()=='-y':
        r = np.array([[np.cos(rad), 0, -np.sin(rad)],[0, 1, 0],[np.sin(rad), 0, np.cos(rad)]]) 
    elif axis.lower()=='-z':
        r = np.array([[np.cos(rad), np.sin(rad), 0],[-np.sin(rad), np.cos(rad), 0],[0, 0, 1]])
    else:
        print('axis should be either x, y, or z')
    return np.dot(r, v)


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        self._target = DEFAULT_PIX
        self._background = 0
        self.fig = Figure(figsize=(5, 6))
        self.fig.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)
        self.ax = self.fig.add_subplot(111, projection='3d')
        super().__init__(self.fig)
        self.setParent(parent)

        self._energy = DEFAULT_GEOMETRY['energy']
        self._delta = DEFAULT_GEOMETRY['delta']
        self._nu = DEFAULT_GEOMETRY['nu']
        self._mu = DEFAULT_GEOMETRY['mu']
        self._eta = DEFAULT_GEOMETRY['eta']
        self._chi = DEFAULT_GEOMETRY['chi']
        self._phi = DEFAULT_GEOMETRY['phi']
        self._ub_energy = DEFAULT_SAMPLE['UB_energy']
        self._UB = DEFAULT_SAMPLE['UB']
        self._or0 = DEFAULT_SAMPLE['or0']
        self._or1 = DEFAULT_SAMPLE['or1']
        self._sample_name = DEFAULT_SAMPLE['name']
        self._lattice = DEFAULT_SAMPLE['lattice']
        self.hkl = [0, 0, 0]
        self.plot()
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_plot)
        self.timer.start(DEFAULT_PV_TIME)

    def plot(self):
        self.plot_beam()
        self.plot_lattice()
        self.hkl = self._calc_hkl()
        self.get_default_range()
        self.plt_horizontal_plane()
        self.zoom_reset()
        self.ax.set_title(f'{self._sample_name}\nHKL=({self.hkl[0]:.3f}, {self.hkl[1]:.3f}, {self.hkl[2]:.3f})')
        self.ax.axis('off')
        self.ax.set_aspect('equal')
        self.fig.tight_layout()
        self.draw()

    def update_plot(self):
        self.ax.clear()
        self.plot_beam()
        self.plot_lattice()
        self.hkl = self._calc_hkl()
        self.plt_horizontal_plane()
        self.zoom_set()
        self.ax.set_title(f'{self._sample_name}\nHKL=({self.hkl[0]:.3f}, {self.hkl[1]:.3f}, {self.hkl[2]:.3f})')
        self.ax.axis('off')
        self.ax.set_aspect('equal')
        self.fig.tight_layout()
        self.draw()

    def plot_beam(self):
        wl = 1.2398e4/self._energy 
        k_wl = 2*np.pi/wl # A^-1
        pol_ratio = 0.5
        
        self.incident_beam = np.array([-1, 0, 0, 1, 0, 0])*k_wl
        self.diffracted_beam = np.hstack([[0, 0, 0], rot(rot(np.array([1, 0, 0]), '-y', self._delta)*k_wl, 'z', self._nu)])
        self.q_beam = np.hstack([[0, 0, 0], self.diffracted_beam[3:]-self.incident_beam[3:]])
        beam_vectors = np.stack([self.incident_beam, self.diffracted_beam, self.q_beam])
        colors = ['k', 'k', 'k']
        linestyles = ['dashed', 'dashed', 'solid']
        arrowstyles = ['-|>', '-|>', '-|>']
        labels = ['', '', 'G']
        names = ['k_in', 'k_out', 'G']

        self.beam_range = np.max([np.linalg.norm(self.q_beam[3:]), np.linalg.norm(self.incident_beam[3:])])
        X, Y, Z, U, V, W = zip(*beam_vectors)
        self.beam_obj = {}
        for x, y, z, u, v, w, color, linestyle, arrowstyle, label, name in zip(X, Y, Z, U, V, W, colors, linestyles, arrowstyles, labels, names):
            ob = self.ax.arrow3D(x, y, z, u, v, w, 
                            color = color,
                            mutation_scale=10,
                            arrowstyle=arrowstyle,
                            linestyle=linestyle,
                            label=label)
            self.beam_obj[name] = ob

    def cal_rlattice(self):
        ra = np.dot(self._UB, [1, 0, 0])
        rb = np.dot(self._UB, [0, 1, 0])
        rc = np.dot(self._UB, [0, 0, 1])

        ra_phi = rot(ra, '-z', self._phi)
        ra_chi = rot(ra_phi, 'x', self._chi)
        ra_eta = rot(ra_chi, '-y', self._eta)
        ra_mu = rot(ra_eta, 'z', self._mu)

        rb_phi = rot(rb, '-z', self._phi)
        rb_chi = rot(rb_phi, 'x', self._chi)
        rb_eta = rot(rb_chi, '-y', self._eta)
        rb_mu = rot(rb_eta, 'z', self._mu)

        rc_phi = rot(rc, '-z', self._phi)
        rc_chi = rot(rc_phi, 'x', self._chi)
        rc_eta = rot(rc_chi, '-y', self._eta)
        rc_mu = rot(rc_eta, 'z', self._mu)

        ra_now = np.hstack([[0, 0, 0], ra_mu])
        rb_now = np.hstack([[0, 0, 0], rb_mu])
        rc_now = np.hstack([[0, 0, 0], rc_mu])

        self._latt = np.stack([ra_now, rb_now, rc_now])

    def plot_lattice(self):
        spin_ratio = 0.7
        self.cal_rlattice()
        ra_now, rb_now, rc_now = self._latt
        X, Y, Z, U, V, W = zip(*self._latt)
        colors = ['indianred', 'limegreen', 'b']
        linestyles = ['solid', 'solid', 'solid']
        labels = ['ra', 'rb', 'rc']
        names = ['ra', 'rb', 'rc']
        self.lattice_obj = {}
        for x, y, z, u, v, w, color, linestyle, label, name in zip(X, Y, Z, U, V, W, colors, linestyles, labels, names):
            ob = self.ax.arrow3D(x, y, z, u, v, w, 
                            color = color,
                            mutation_scale=15,
                            arrowstyle="-|>",
                            linestyle=linestyle,
                            label=label)
            self.lattice_obj[name] = ob

        unit_cell = []
        unit_cell.append(np.hstack([self._latt[0, :3]+self._latt[1, 3:], self._latt[0, 3:]]))
        unit_cell.append(np.hstack([self._latt[0, :3]+self._latt[2, 3:], self._latt[0, 3:]]))
        unit_cell.append(np.hstack([self._latt[0, :3]+self._latt[1, 3:]+self._latt[2, 3:], self._latt[0, 3:]]))
        unit_cell.append(np.hstack([self._latt[1, :3]+self._latt[0, 3:], self._latt[1, 3:]]))
        unit_cell.append(np.hstack([self._latt[1, :3]+self._latt[2, 3:], self._latt[1, 3:]]))
        unit_cell.append(np.hstack([self._latt[1, :3]+self._latt[0, 3:]+self._latt[2, 3:], self._latt[1, 3:]]))
        unit_cell.append(np.hstack([self._latt[2, :3]+self._latt[0, 3:], self._latt[2, 3:]]))
        unit_cell.append(np.hstack([self._latt[2, :3]+self._latt[1, 3:], self._latt[2, 3:]]))
        unit_cell.append(np.hstack([self._latt[2, :3]+self._latt[0, 3:]+self._latt[1, 3:], self._latt[2, 3:]]))
        X, Y, Z, U, V, W = zip(*unit_cell)
        colors = ['k']*len(unit_cell)
        linestyles = ['dashed']*len(unit_cell)
        uc_names = [f'uc_line_{i}' for i in range(9)]
        for x, y, z, u, v, w, color, linestyle, name in zip(X, Y, Z, U, V, W, colors, linestyles, uc_names):
            ob = self.ax.arrow3D(x, y, z, u, v, w, 
                            color = color,
                            arrowstyle="-",
                            linestyle=linestyle,
                            lw=0.5)
            self.lattice_obj[name] = ob
        
        or_names = ['or0', 'or1']
        ors = [self._or0, self._or1]
        for o, name in zip(ors, or_names):
            or_h, or_k, or_l = o[:3]
            or_now = or_h*ra_now[3:] + or_k*rb_now[3:] + or_l*rc_now[3:]
            ob = self.ax.scatter(or_now[0], or_now[1], or_now[2], 
                            marker='o', 
                            color='orange', 
                            s=30)
            self.lattice_obj[name] = ob
        
        try:
            hkl_now = self.h_target*ra_now[3:] + self.k_target*rb_now[3:] + self.l_target*rc_now[3:]
            ob = self.ax.scatter(hkl_now[0], hkl_now[1], hkl_now[2], 
                            marker='o', 
                            color='red', 
                            s=20)
            self.lattice_obj['target_hkl'] = ob
        except Exception as e:
            pass

    
    def _calc_hkl(self):
        ra, rb, rc = self._latt
        ra_now = ra[3:]
        rb_now = rb[3:]
        rc_now = rc[3:]
        q_beam = self.q_beam[3:]
        d = np.vstack([ra_now, rb_now, rc_now]).T
        h, k, l = np.dot(np.linalg.inv(d), q_beam)
        return h, k, l
    
    def get_default_range(self):
        ra, rb, rc = self._latt
        ra_norm = np.linalg.norm(ra)
        rb_norm = np.linalg.norm(rb)
        rc_norm = np.linalg.norm(rc)
        self.default_plt_range =  1.2*max(ra_norm, rb_norm, rc_norm)
        self.plt_range = 1.2*max(ra_norm, rb_norm, rc_norm)

    def zoom_in(self):
        self.plt_range *= 0.8

    def zoom_out(self):
        self.plt_range *= 1.25
    
    def zoom_reset(self):
        self.plt_range = self.default_plt_range

    def zoom_set(self):
        self.ax.set_xlim(-self.plt_range, self.plt_range)
        self.ax.set_ylim(-self.plt_range, self.plt_range)
        self.ax.set_zlim(-self.plt_range, self.plt_range)

    def plt_horizontal_plane(self):
        plane_x = np.linspace(-1.2*self.plt_range, 1.2*self.plt_range, 3)
        plane_y = np.linspace(-1.2*self.plt_range, 1.2*self.plt_range, 3)
        plane_X, plane_Y = np.meshgrid(plane_x, plane_y)
        plane_Z = np.zeros_like(plane_X)
        self.ax.plot_surface(plane_X, plane_Y, plane_Z, cmap='viridis', alpha=0.2)


class MonitorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('8-ID Psi-View')
        self.setFixedSize(WINDOW_WIDTH, WINDOW_HEIGHT)
        self.generalLayout = QVBoxLayout()
        centralWidget = QWidget(self)
        centralWidget.setLayout(self.generalLayout)
        self.setCentralWidget(centralWidget)
        self.isMotorOnline = True
        self.isSampleOnline = True
        self.isXrayOnline = True
        self._createQView()
        self._createZoomButtons()
        self._createStatusDisplay()
        self._createMotorDisplay()
        self._createLatticeDisplay()
        self._createReflectionDisplay()
        self._createUBTables()
        self._createModeMenu()
        self._createCalcDisplay()
        self._createButtons()
        # self.diffractometer_constraints = {
        #     # axis: Constraint(lo_limit, hi_limit, value, fit)
        #     "omega": Constraint(-180, 180, 0, True),
        #     "chi": Constraint(-180, 180, 0, True),
        #     "phi": Constraint(-180, 180, 0, True),
        #     "mu": Constraint(-180, 180, 0, True),
        #     "gamma": Constraint(-180, 180, 0, True),
        #     "delta": Constraint(-180, 180, 0, True),
        # }

        sim6c.core.constraints["omega"].limits = (-180, 180)
        sim6c.core.constraints["chi"].limits   = (-90, 90)
        sim6c.core.constraints["phi"].limits   = (-180, 180)
        sim6c.core.constraints["mu"].limits    = (-180, 180)
        sim6c.core.constraints["gamma"].limits = (-180, 180)
        sim6c.core.constraints["delta"].limits = (-180, 180)

        self.diffractometer_constraints = {
            "omega": (-180, 180),
            "chi":   (-90, 90),
            "phi":   (-180, 180),
            "mu":    (-180, 180),
            "gamma": (-180, 180),
            "delta": (-180, 180),
        }
        
        self._psi_sim = hklpy2.creator(geometry="E6C", name="psi_sim")
        self._psi_sim.core.mode = "constant_mu_horizontal"
        self._set_soft_limits()
        # self._psi_sim = SimulatedE6C('', name="psi_sim")
        # self._psi_sim.calc.engine.mode = "constant_mu_horizontal"


        #### possible modes ####
        # bissector_vertical, constant_omega_vertical, constant_chi_vertical, 
        # constant_phi_vertical, lifting_detector_phi, lifting_detector_omega, 
        # lifting_detector_mu, double_diffraction_vertical, bissector_horizontal, 
        # double_diffraction_horizontal, psi_constant_vertical, psi_constant_horizontal, 
        # constant_mu_horizontal

        self._set_soft_limits()
        self.initPV()
        self.delta_target, self.nu_target, self.mu_target, self.eta_target, self.chi_target, self.phi_target =  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        self.updateCalcDisplayLayout()
        # self._createButtons()
        self.pvtimer = QTimer(self)
        self.pvtimer.timeout.connect(self.syncPV)
        self.pvtimer.start(DEFAULT_PV_TIME)

        self._calc_UB()

    # def _set_soft_limits(self):
    #     for axis, (lo, hi) in self.diffractometer_constraints.items():
    #         self._psi_sim.core.constraints[axis].limits = (lo, hi)

    def _set_soft_limits(self):
        for axis, (lo, hi) in self.diffractometer_constraints.items():
            self._psi_sim.core.constraints[axis].limits = (lo, hi)
            sim6c.core.constraints[axis].limits = (lo, hi)
            

    # def _set_soft_limits(self):
    #     self._psi_sim.apply_constraints(self.diffractometer_constraints)
    #     for motor in self.diffractometer_constraints.keys():
    #         lo_limit = self.diffractometer_constraints[motor].low_limit
    #         hi_limit = self.diffractometer_constraints[motor].high_limit
    #         value = self.diffractometer_constraints[motor].value
    #         setattr(getattr(self._psi_sim, motor), '_limits', (lo_limit, hi_limit))
    #         setattr(getattr(self._psi_sim, motor), '_position', value)

    def pause_pvtimer(self):
        if self.pvtimer.isActive():
            self.pvtimer.stop()

    def resume_pvtimer(self):
        if not self.pvtimer.isActive():
            self.pvtimer.start(DEFAULT_PV_TIME)

    def _createQView(self):
        self.qviewLayout = QVBoxLayout()
        self.canvas = MplCanvas(self)
        self.qviewLayout.addWidget(self.canvas)
        self.generalLayout.addLayout(self.qviewLayout)

    def _createZoomButtons(self):
        self.zoomButtonMap = {}
        self.zoomButtonLayout = QGridLayout()
        names = ['zoom_in', 'zoom_out', 'zoom_reset']
        labels = ['zoom in', 'zoom out', 'reset']


        for col, (name, label) in enumerate(zip(names, labels)):
            button = QPushButton(label)
            button.setFixedSize(BUTTON_WIDTH, BUTTON_HEIGHT)
            self.zoomButtonMap[name] = button
            self.zoomButtonLayout.addWidget(button, 0, col)

            if name == 'zoom_in':
                button.setStyleSheet("""
                QPushButton {
                    background-color: #2196F3;
                    color: white;
                    font-weight: bold;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: #1976D2;
                }
                """)
            elif name == 'zoom_out':
                button.setStyleSheet("""
                QPushButton {
                    background-color: #FF9800;
                    color: white;
                    font-weight: bold;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: #F57C00;
                }
                """)
            elif name == 'zoom_reset':
                button.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50;
                    color: white;
                    font-weight: bold;
                    border-radius: 5px;
                }
                QPushButton:hover {
                    background-color: #45a049;
                }
                """)


        self.generalLayout.addLayout(self.zoomButtonLayout)


    def _createStatusDisplay(self):
        self.statusDisplayMap = {}
        self.statusDisplayLayout = QGridLayout()
        self.statusDisplayLayout.setVerticalSpacing(0)
        self.statusDisplayLayout.setHorizontalSpacing(5)
        self.statusDisplayLayout.setContentsMargins(0, 0, 0, 0)

        display_names = ['xray_status', 'motor_status', 'ub_status']
        display_labels = ['X-ray', 'Motor', 'UB']

        for col, label_text in enumerate(display_labels):
            label = QLabel(label_text)
            label.setAlignment(Qt.AlignCenter)
            self.statusDisplayLayout.addWidget(label, 0, col)

        for col, display_name in enumerate(display_names):
            line_edit = QLineEdit(display_name)
            line_edit.setFixedHeight(DISPLAY_HEIGHT)
            line_edit.setAlignment(Qt.AlignCenter)
            self.statusDisplayMap[display_name] = line_edit
            self.statusDisplayLayout.addWidget(line_edit, 1, col)
        self.statusDisplayMap['xray_status'].setReadOnly(True)
        self.statusDisplayMap['motor_status'].setReadOnly(True)
        self.statusDisplayMap['ub_status'].setReadOnly(True)
        self.generalLayout.addLayout(self.statusDisplayLayout)

    def _createMotorDisplay(self):
        self.motorDisplayMap = {}
        self.motorDisplayLayout = QGridLayout()
        self.motorDisplayLayout.setVerticalSpacing(0)
        self.motorDisplayLayout.setHorizontalSpacing(5)
        self.motorDisplayLayout.setContentsMargins(0, 0, 0, 0)

        motor_names = ['energy', 'polarization', 'delta', 'nu', 'mu', 'eta', 'chi', 'phi']
        motor_labels = ['Energy (eV)', 'Polarization', 'delta', 'nu', 'mu', 'eta', 'chi', 'phi']

        for col, label_text in enumerate(motor_labels):
            label = QLabel(label_text)
            label.setAlignment(Qt.AlignCenter)
            self.motorDisplayLayout.addWidget(label, 0, col)

        for col, display_name in enumerate(motor_names):
            line_edit = QLineEdit(display_name)
            line_edit.setFixedHeight(DISPLAY_HEIGHT)
            line_edit.setAlignment(Qt.AlignCenter)
            self.motorDisplayMap[display_name] = line_edit
            self.motorDisplayLayout.addWidget(line_edit, 1, col)

        self.motorDisplayMap['polarization'].setReadOnly(True)
        self.generalLayout.addLayout(self.motorDisplayLayout)

    def _createLatticeDisplay(self):
        self.latticeDisplayMap = {}
        self.latticeDisplayLayout = QGridLayout()
        self.latticeDisplayLayout.setVerticalSpacing(0)
        self.latticeDisplayLayout.setHorizontalSpacing(5)
        self.latticeDisplayLayout.setContentsMargins(0, 0, 0, 0)

        display_names = ['sample', 'a', 'b', 'c', 'alpha', 'beta', 'gamma']
        display_labels = ['Sample', 'a', 'b', 'c', 'alpha', 'beta', 'gamma']

        for col, label_text in enumerate(display_labels):
            label = QLabel(label_text)
            label.setAlignment(Qt.AlignCenter)
            self.latticeDisplayLayout.addWidget(label, 0, col)

        for col, display_name in enumerate(display_names):
            line_edit = QLineEdit(display_name)
            line_edit.setFixedHeight(DISPLAY_HEIGHT)
            line_edit.setAlignment(Qt.AlignCenter)
            if display_name == 'a':
                line_edit.setStyleSheet('color: red;')
            elif display_name == 'b':
                line_edit.setStyleSheet('color: green;')
            elif display_name == 'c':
                line_edit.setStyleSheet('color: blue;')
            self.latticeDisplayMap[display_name] = line_edit
            self.latticeDisplayLayout.addWidget(line_edit, 1, col)

        self.generalLayout.addLayout(self.latticeDisplayLayout)

    def _createReflectionDisplay(self):
        self.reflectionDisplayMap = {}
        self.reflectionDisplayLayout = QGridLayout()
        self.reflectionDisplayLayout.setVerticalSpacing(0)
        self.reflectionDisplayLayout.setHorizontalSpacing(5)
        self.reflectionDisplayLayout.setContentsMargins(0, 0, 0, 0)

        display_labels = ['Reflection', 'H', 'K', 'L', 'delta', 'nu', 'mu', 'eta', 'chi', 'phi']
        display_names = [['or0', 'H0', 'K0', 'L0', 'delta0', 'nu0', 'mu0', 'eta0', 'chi0', 'phi0'], 
                          ['or1', 'H1', 'K1', 'L1', 'delta1', 'nu1', 'mu1', 'eta1', 'chi1', 'phi1']]

        for col, label_text in enumerate(display_labels):
            label = QLabel(label_text)
            label.setAlignment(Qt.AlignCenter)
            self.reflectionDisplayLayout.addWidget(label, 0, col)

        for row, reflection in enumerate(display_names):
            for col, display_name in enumerate(reflection):
                line_edit = QLineEdit(display_name)
                line_edit.setFixedHeight(DISPLAY_HEIGHT)
                line_edit.setAlignment(Qt.AlignCenter)
                self.reflectionDisplayMap[display_name] = line_edit
                self.reflectionDisplayLayout.addWidget(line_edit, row+1, col)

        self.generalLayout.addLayout(self.reflectionDisplayLayout)

    def _createUBTables(self):
        self.tablesLayout = QGridLayout()
        self._u_table = QTableWidget(3, 3)
        self._ub_table = QTableWidget(3, 3)
        self._u_table.setFixedSize(350, 120)
        self._ub_table.setFixedSize(350, 120)
        self._u_label = QLabel('U')
        self._ub_label = QLabel('UB')
        self._u_label.setAlignment(Qt.AlignCenter)
        self._ub_label.setAlignment(Qt.AlignCenter)
        self.tablesLayout.addWidget(self._u_label, 0, 0)
        self.tablesLayout.addWidget(self._ub_label, 0, 1)
        self.tablesLayout.addWidget(self._u_table, 1, 0)
        self.tablesLayout.addWidget(self._ub_table, 1, 1)

        self.generalLayout.addLayout(self.tablesLayout)

    def _createModeMenu(self):
        self.modeMap = {}
        self.modeLayout = QGridLayout()
        label = QLabel('Mode')
        label.setAlignment(Qt.AlignCenter)
        self.modeLayout.addWidget(label, 0, 0)
        combo_box = QComboBox()
        for c in PSI_MODES:
            combo_box.addItem(c)
        self.modeMap['mode'] = combo_box
        self.modeLayout.addWidget(combo_box, 1, 0)
        self.generalLayout.addLayout(self.modeLayout)

    def _createCalcDisplay(self):
        self.calcDisplayMap = {}
        self.calcDisplayLayout = QGridLayout()
        self.calcDisplayLayout.setVerticalSpacing(0)
        self.calcDisplayLayout.setHorizontalSpacing(5)
        self.calcDisplayLayout.setContentsMargins(0, 0, 0, 0)

        labels = ['H', 'K', 'L', 'delta', 'nu', 'mu', 'eta', 'chi', 'phi']
        names = ['h_target', 'k_target', 'l_target', 'delta_target', 'nu_target', 
                 'mu_target', 'eta_target', 'chi_target', 'phi_target']

        for col, label_text in enumerate(labels):
            label = QLabel(label_text)
            label.setAlignment(Qt.AlignCenter)
            self.calcDisplayLayout.addWidget(label, 0, col)

        for col, name in enumerate(names):
            line_edit = QLineEdit(name)
            line_edit.setFixedHeight(DISPLAY_HEIGHT)
            line_edit.setAlignment(Qt.AlignCenter)
            if name in ['delta', 'nu', 'mu', 'eta', 'chi', 'phi']:
                line_edit.setReadOnly(True)
            self.calcDisplayMap[name] = line_edit
            self.calcDisplayLayout.addWidget(line_edit, 1, col)

        self.generalLayout.addLayout(self.calcDisplayLayout)

    def _createButtons(self):
        self.buttonMap = {}
        self.buttonsLayout = QGridLayout()
        keys = ['sync', 'quit']

        for col, key in enumerate(keys):
            button = QPushButton(key)
            button.setFixedSize(BUTTON_WIDTH, BUTTON_HEIGHT)

            # Apply styles based on button type
            if key == 'sync':
                button.setStyleSheet("""
                    QPushButton {
                        background-color: #4CAF50;  /* Green */
                        color: white;
                        font-weight: bold;
                        border-radius: 5px;
                    }
                    QPushButton:hover {
                        background-color: #45a049;
                    }
                """)
            elif key == 'quit':
                button.setStyleSheet("""
                    QPushButton {
                        background-color: #000000;  /* Black */
                        color: white;
                        font-weight: bold;
                        border-radius: 5px;
                    }
                    QPushButton:hover {
                        background-color: #333333;
                    }
                """)

            self.buttonMap[key] = button
            self.buttonsLayout.addWidget(button, 0, col)
        self.generalLayout.addLayout(self.buttonsLayout)

    def updateMotorDisplayLayout(self, pos_lst):
        for name, pos in zip(['delta', 'nu', 'mu', 'eta', 'chi', 'phi'], pos_lst):
            self.motorDisplayMap[name].setText(f'{pos:.3f}')
        self._delta, self._nu, self._mu, self._eta, self._chi, self._phi = pos_lst

    def getMotorDisplayText(self):
        return np.array([
            float(self.motorDisplayMap['delta'].text()),
            float(self.motorDisplayMap['nu'].text()),
            float(self.motorDisplayMap['mu'].text()),
            float(self.motorDisplayMap['eta'].text()),
            float(self.motorDisplayMap['chi'].text()),
            float(self.motorDisplayMap['phi'].text()),
        ])
    
    def updateXrayDisplayLayout(self, name, value):
        if name == 'energy':
            self.motorDisplayMap[name].setText(f'{value:.3f}')
            self._energy = value
            # self._psi_sim.calc.wavelength = 12398.4/self._energy
            self._psi_sim.core.wavelength = 12398.4 / self._energy
        elif name == 'polarization':
            self.motorDisplayMap[name].setText(value)
            self._polarization = value

    def getXrayDisplayLayout(self, name):
        if name == 'energy':
            return float(self.motorDisplayMap[name].text())
        elif name == 'polarization':
            return self.motorDisplayMap[name].text()
    
    def updateStatusDisplayLayout(self, display_name, text):
        self.statusDisplayMap[display_name].setText(str(text))

        if text == SYNCED_MSG:
            self.statusDisplayMap[display_name].setStyleSheet('color: green; font-weight: bold;')
        elif text == NOT_SYNCED_MSG:
            self.statusDisplayMap[display_name].setStyleSheet('color: red; font-weight: bold;')
        else:
            self.statusDisplayMap[display_name].setStyleSheet('')

    def getStatusDisplayText(self, display_name):
        return self.statusDisplayMap[display_name].text()
    
    # def updateLatticeDisplayLayout(self, display_name, lattice):
    #     if display_name == 'lattice':
    #         self._lattice = lattice
    #         # self._psi_sim.calc.sample.lattice = lattice
    #         self._psi_sim.core.sample.lattice = lattice
    #         self._calc_UB()
    #         for n, v in zip(['a', 'b', 'c', 'alpha', 'beta', 'gamma'], lattice):
    #             self.latticeDisplayMap[n].setText(f'{v:.3f}')
    #     elif display_name == 'sample':
    #         self.latticeDisplayMap[display_name].setText(str(lattice))
    #         self._sample = lattice
    
    def updateLatticeDisplayLayout(self, what, lattice):
        if what == "lattice":
            # lattice is [a, b, c, alpha, beta, gamma]
            lat_obj = Lattice(*lattice)
            self._psi_sim.core.sample.lattice = lat_obj

        elif what == "sample":
            self._sample_name = lattice

    # def updateLatticeDisplayLayout(self, what, lattice):
    #     lat_obj = Lattice(*lattice)
    #     self._psi_sim.core.sample.lattice = lat_obj
    #     sim6c.core.sample.lattice = lat_obj

    def getLatticeDisplayText(self, display_name):
        if display_name == 'lattice':
            return np.array([float(self.latticeDisplayMap[n].text()) for n in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']])
        else: 
            return self.latticeDisplayMap[display_name].text()
    
    def updateReflectionDisplayLayout(self, name, reflection):
        if name == 'or0':
            display_names = ['H0', 'K0', 'L0', 'delta0', 'nu0', 'mu0', 'eta0', 'chi0', 'phi0']
            self._or0 = reflection
        elif name == 'or1':
            display_names = ['H1', 'K1', 'L1', 'delta1', 'nu1', 'mu1', 'eta1', 'chi1', 'phi1']
            self._or1 = reflection     
        for n, v, in zip(display_names, reflection):
            self.reflectionDisplayMap[n].setText(f'{v:.3f}')
            # self._psi_sim.core.sample._reflections.clear()        
            # self._psi_sim.calc.sample.clear_reflections()
        for r in [self._or0, self._or1]:
            h, k, l, delta, nu, mu, eta, chi, phi = r
            # self._psi_sim.core.sample.add_reflections(h,k,l,
                                                    # position=self._psi_sim.core.Position(mu=mu,
                                                    #                                         omega=eta,
                                                    #                                         chi=chi,
                                                    #                                         phi=phi,
                                                    #                                         gamma=nu,
                                                    #                                         delta=delta))


    def getReflectionDisplayLayout(self, display_name):
        if display_name == 'or0':
            return np.array([float(self.reflectionDisplayMap[n].text()) for n in ['H0', 'K0', 'L0', 'delta0', 'nu0', 'mu0', 'eta0', 'chi0', 'phi0']])
        elif display_name == 'or1':
            return np.array([float(self.reflectionDisplayMap[n].text()) for n in ['H1', 'K1', 'L1', 'delta1', 'nu1', 'mu1', 'eta1', 'chi1', 'phi1']])
        else:
            return self.reflectionDisplayMap[display_name].text()
        
    def updateSampleDisplay(self, sample):
        name, lattice, or0, or1 = sample
        self.updateLatticeDisplayLayout('sample', name)
        self.updateLatticeDisplayLayout('lattice', lattice)
        self.updateReflectionDisplayLayout('or0', or0)
        self.updateReflectionDisplayLayout('or1', or1)
        
    def getSampleDisplay(self):
        name = self.getLatticeDisplayText('sample')
        lattice = self.getLatticeDisplayText('lattice')
        or0 = self.getReflectionDisplayLayout('or0')
        or1 = self.getReflectionDisplayLayout('or1')
        sample = [name, lattice, or0, or1]
        return sample

    def _calc_UB(self):
        reals0 = self._psi_sim.core.standardize_reals({
            'delta': float(self._or0[3]),
            'gamma': float(self._or0[4]),
            'mu': float(self._or0[5]),
            'omega': float(self._or0[6]),
            'chi': float(self._or0[7]),
            'phi': float(self._or0[8]),
        })
        reals1 = self._psi_sim.core.standardize_reals({
            'delta': float(self._or1[3]),
            'gamma': float(self._or1[4]),
            'mu': float(self._or1[5]),
            'omega': float(self._or1[6]),
            'chi': float(self._or1[7]),
            'phi': float(self._or1[8]),
        })
        or0 = self._psi_sim.add_reflection(
            (self._or0[0], self._or0[1], self._or0[2]),
            reals=reals0,
            name='or0',
            replace=True,
        )
        or1 = self._psi_sim.add_reflection(
            (self._or1[0], self._or1[1], self._or1[2]),
            reals=reals1,
            name='or1',
            replace=True,
        )
        UB_new = self._psi_sim.core.calc_UB(or0, or1)
        self._U = np.asarray(self._psi_sim.sample.U, dtype=float)
        self._UB = np.asarray(self._psi_sim.sample.UB, dtype=float)
        self._UB = np.asarray(UB_new, dtype=float)
        for table, matrix in zip([self._u_table, self._ub_table], [self._U, self._UB]):
            for i in range(3):
                for j in range(3):
                    item = QTableWidgetItem(f'{matrix[i, j]:.3f}')
                    item.setTextAlignment(Qt.AlignCenter)
                    table.setItem(i, j, item)

    
    def _calc_motors(self):
        self._psi_sim.delta.move(self._delta)
        self._psi_sim.gamma.move(self._nu)
        self._psi_sim.mu.move(self._mu)
        self._psi_sim.omega.move(self._eta)
        self._psi_sim.chi.move(self._chi)
        self._psi_sim.phi.move(self._phi)
        self.h_target = self.calcDisplayMap['h_target'].text()
        self.k_target = self.calcDisplayMap['k_target'].text()
        self.l_target = self.calcDisplayMap['l_target'].text()
        # self.phi_target = self.calcDisplayMap['phi_target'].text()
        try:
            self.h_target = float(self.h_target)
            self.k_target = float(self.k_target)
            self.l_target = float(self.l_target)
            self.canvas.h_target = self.h_target
            self.canvas.k_target = self.k_target
            self.canvas.l_target = self.l_target
            self._solution = self._psi_sim.forward((self.h_target,self.k_target,self.l_target))
            self.delta_target, self.nu_target, self.mu_target, self.eta_target, self.chi_target, self.phi_target =  np.round(self._solution,5)
        except:
            self.delta_target, self.nu_target, self.mu_target, self.eta_target, self.chi_target, self.phi_target =  np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        self.updateCalcDisplayLayout()

    def change_psi_mode(self):
        new_mode = self.modeMap['mode'].currentText()
        self._psi_sim.core.mode = new_mode
        self._calc_motors()

    def updateCalcDisplayLayout(self):
        self.calcDisplayMap['delta_target'].setText(f'{self.delta_target:.3f}')
        self.calcDisplayMap['nu_target'].setText(f'{self.nu_target:.3f}')
        self.calcDisplayMap['mu_target'].setText(f'{self.mu_target:.3f}')
        self.calcDisplayMap['eta_target'].setText(f'{self.eta_target:.3f}')
        self.calcDisplayMap['chi_target'].setText(f'{self.chi_target:.3f}')
        self.calcDisplayMap['phi_target'].setText(f'{self.phi_target:.3f}')

    def updateUBDisplayLayout(self, name, reflection):
        if name == 'or0':
            display_names = ['H0', 'K0', 'L0', 'delta0', 'nu0', 'mu0', 'eta0', 'chi0', 'phi0']
            self._or0 = reflection
        elif name == 'or1':
            display_names = ['H1', 'K1', 'L1', 'delta1', 'nu1', 'mu1', 'eta1', 'chi1', 'phi1']
            self._or1 = reflection
        for n, v, in zip(display_names, reflection):
            self.reflectionDisplayMap[n].setText(f'{v:.3f}')
        self._calc_UB()
        self.syncQview()

    def syncMotors(self):
        if self.isMotorOnline:
            try:
                self._delta = float(caget(DELTA_RBV, timeout=0.1))
                self._nu = float(caget(NU_RBV, timeout=0.1))
                self._mu = float(caget(MU_RBV, timeout=0.1))
                self._eta = float(caget(ETA_RBV, timeout=0.1))
                self._chi = float(caget(CHI_RBV, timeout=0.1))
                self._phi = float(caget(PHI_RBV, timeout=0.1))
                self._psi_sim.delta.move(self._delta)
                self._psi_sim.gamma.move(self._nu)
                self._psi_sim.mu.move(self._mu)
                self._psi_sim.omega.move(self._eta)
                self._psi_sim.chi.move(self._chi)
                self._psi_sim.phi.move(self._phi)
                self._motor_synced = True
                self.updateStatusDisplayLayout('motor_status', SYNCED_MSG)
                self.isMotorOnline = True
            except Exception as e:
                self.isMotorOnline = False
                self.loadDefaultMotors()
        else:
            self.loadDefaultMotors()
        self.updateMotorDisplayLayout([self._delta, self._nu, self._mu, self._eta, self._chi, self._phi])

    def loadDefaultMotors(self):
        self.updateStatusDisplayLayout('motor_status', NOT_SYNCED_MSG)
        self._delta = DEFAULT_GEOMETRY['delta']
        self._nu = DEFAULT_GEOMETRY['nu']
        self._mu = DEFAULT_GEOMETRY['mu']
        self._eta = DEFAULT_GEOMETRY['eta']
        self._chi = DEFAULT_GEOMETRY['chi']
        self._phi = DEFAULT_GEOMETRY['phi']
        self._psi_sim.delta.move(self._delta)
        self._psi_sim.gamma.move(self._nu)
        self._psi_sim.mu.move(self._mu)
        self._psi_sim.omega.move(self._eta)
        self._psi_sim.chi.move(self._chi)
        self._psi_sim.phi.move(self._phi)
        self._motor_synced = False

    def syncXray(self):
        if self.isXrayOnline:
            try:
                self._energy = float(caget(ENERGY_RBV, timeout=0.1))
                self._polarization = int(caget(POLARIZATION_RBV, timeout=0.1))
                self._xray_synced= True
                self.updateStatusDisplayLayout('xray_status', SYNCED_MSG)
                self.isXrayOnline = True
            except Exception as e:
                self.isXrayOnline = False
                self.loadDefaultXray()
        else:
            self.loadDefaultXray()
        self.updateXrayDisplayLayout('energy', self._energy)
        self.updateXrayDisplayLayout('polarization', POLARIZATION_DICT[self._polarization])

    def loadDefaultXray(self):
            self._energy = DEFAULT_GEOMETRY['energy']
            self._polarization = DEFAULT_GEOMETRY['polarization']
            self.updateStatusDisplayLayout('xray_status', NOT_SYNCED_MSG)
            self._xray_synced = False

    def syncSample(self):
        if self.isSampleOnline:
            try:
                self._ub_energy = float(caget(UB_ENERGY_RBV, timeout=0.1))
                self._or0 = caget(UB_OR0_RBV, timeout=0.1)
                self._or1 = caget(UB_OR1_RBV, timeout=0.1)
                self._lattice = caget(LATTICE_RBV, timeout=1)
                self._sample_name = caget(SAMPLE_NAME_RBV, timeout=0.1)
                self._UB = caget(UB_MATRIX_RBV, timeout=1).reshape(3, 3)
                self._ub_synced = True
                self.updateStatusDisplayLayout('ub_status', SYNCED_MSG)
                self.isSampleOnline = True
            except Exception as e:
                self.isSampleOnline = False
                self.loadDefaultSample()
        else:
            self.loadDefaultSample()

        self.updateLatticeDisplayLayout('lattice', self._lattice)
        self.updateLatticeDisplayLayout('sample', self._sample_name)
        self.updateReflectionDisplayLayout('or0', self._or0)
        self.updateReflectionDisplayLayout('or1', self._or1)

    def loadDefaultSample(self):
        self._ub_energy = DEFAULT_SAMPLE['UB_energy']
        self._UB = DEFAULT_SAMPLE['UB']
        self._or0 = DEFAULT_SAMPLE['or0']
        self._or1 = DEFAULT_SAMPLE['or1']
        self._lattice = DEFAULT_SAMPLE['lattice']
        self._sample_name = DEFAULT_SAMPLE['name']
        self._ub_synced = False
        self.updateStatusDisplayLayout('ub_status', NOT_SYNCED_MSG)


    def syncQview(self):
        self.canvas._delta = self._delta
        self.canvas._nu = self._nu
        self.canvas._mu = self._mu
        self.canvas._eta = self._eta
        self.canvas._chi = self._chi
        self.canvas._phi = self._phi
        self.canvas._energy = self._energy
        self.canvas._UB = self._UB
        self.canvas._ub_energy = self._ub_energy
        self.canvas._or0 = self._or0
        self.canvas._or1 = self._or1
        self.canvas._sample_name = self._sample_name
        self.canvas._lattice = self._lattice

    def initPV(self):
        self.syncMotors()
        self.syncXray()
        self.syncSample()

    def syncPV(self):
        self.syncMotors()
        self.syncXray()
        self.syncSample()
        self.syncQview()


class Qview:
    def __init__(self, view):
        self._view = view
        self._connectSignalsAndSlots()
        # self._aligntimer = QTimer()
        # self._aligntimer.timeout.connect(partial(self._align_method, target=self._target))

    def _read_positions(self, motor):
        self._delta, self._nu, self._mu, self._eta, self._chi, self._phi = self._view.getMotorDisplayText()
        self._view.updateMotorDisplayLayout([self._delta, self._nu, self._mu, self._eta, self._chi, self._phi])
        self._view.syncQview()
        self._view._psi_sim.delta.move(self._delta)
        self._view._psi_sim.gamma.move(self._nu)
        self._view._psi_sim.mu.move(self._mu)
        self._view._psi_sim.omega.move(self._eta)
        self._view._psi_sim.chi.move(self._chi)
        self._view._psi_sim.phi.move(self._phi)
        self._view.updateStatusDisplayLayout('motor_status', NOT_SYNCED_MSG)

    def _read_xray(self):
        self._energy = self._view.getXrayDisplayLayout('energy')
        self._view.updateXrayDisplayLayout('energy', self._energy)
        self._view.syncQview()
        self._view._calc_motors()
        self._view.updateStatusDisplayLayout('xray_status', NOT_SYNCED_MSG)

    def _read_sample(self):
        self._name, self._lattice, self._or0, self._or1 = self._view.getSampleDisplay()
        self._view.updateSampleDisplay([self._name, self._lattice, self._or0, self._or1])
        self._view.syncQview()
        self._view._calc_motors()
        self._view.updateStatusDisplayLayout('ub_status', NOT_SYNCED_MSG)

    def _sync_all(self):
        self._view.resume_pvtimer()
        self._view._calc_motors()
        self._view.updateStatusDisplayLayout('motor_status', SYNCED_MSG)
        self._view.updateStatusDisplayLayout('xray_status', SYNCED_MSG)
        self._view.updateStatusDisplayLayout('ub_status', SYNCED_MSG)

    def _connectSignalsAndSlots(self):
        self._view.motorDisplayMap['energy'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['energy'].returnPressed.connect(self._read_xray)

        self._view.motorDisplayMap['delta'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['delta'].returnPressed.connect(partial(self._read_positions, motor='delta'))

        self._view.motorDisplayMap['nu'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['nu'].returnPressed.connect(partial(self._read_positions, motor='nu'))

        self._view.motorDisplayMap['mu'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['mu'].returnPressed.connect(partial(self._read_positions, motor='mu'))

        self._view.motorDisplayMap['eta'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['eta'].returnPressed.connect(partial(self._read_positions, motor='eta'))

        self._view.motorDisplayMap['chi'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['chi'].returnPressed.connect(partial(self._read_positions, motor='chi'))

        self._view.motorDisplayMap['phi'].textEdited.connect(self._view.pause_pvtimer)
        self._view.motorDisplayMap['phi'].returnPressed.connect(partial(self._read_positions, motor='phi'))

        self._view.latticeDisplayMap['sample'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['sample'].returnPressed.connect(self._read_sample)

        self._view.latticeDisplayMap['a'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['a'].returnPressed.connect(self._read_sample)

        self._view.latticeDisplayMap['b'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['b'].returnPressed.connect(self._read_sample)

        self._view.latticeDisplayMap['c'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['c'].returnPressed.connect(self._read_sample)

        self._view.latticeDisplayMap['alpha'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['alpha'].returnPressed.connect(self._read_sample)

        self._view.latticeDisplayMap['beta'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['beta'].returnPressed.connect(self._read_sample)

        self._view.latticeDisplayMap['gamma'].textEdited.connect(self._view.pause_pvtimer)
        self._view.latticeDisplayMap['gamma'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['H0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['H0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['K0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['K0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['L0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['L0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['delta0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['delta0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['nu0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['nu0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['mu0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['mu0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['eta0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['eta0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['chi0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['chi0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['phi0'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['phi0'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['H1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['H1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['K1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['K1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['L1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['L1'].returnPressed.connect(self._read_sample)


        self._view.reflectionDisplayMap['delta1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['delta1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['nu1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['nu1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['mu1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['mu1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['eta1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['eta1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['chi1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['chi1'].returnPressed.connect(self._read_sample)

        self._view.reflectionDisplayMap['phi1'].textEdited.connect(self._view.pause_pvtimer)
        self._view.reflectionDisplayMap['phi1'].returnPressed.connect(self._read_sample)


        self._view.zoomButtonMap['zoom_in'].clicked.connect(self._view.canvas.zoom_in)
        self._view.zoomButtonMap['zoom_out'].clicked.connect(self._view.canvas.zoom_out)
        self._view.zoomButtonMap['zoom_reset'].clicked.connect(self._view.canvas.zoom_reset)


        self._view.buttonMap['sync'].clicked.connect(self._sync_all)
        self._view.buttonMap['quit'].clicked.connect(QApplication.quit)

        self._view.modeMap['mode'].currentTextChanged.connect(self._view.change_psi_mode)

        self._view.calcDisplayMap['h_target'].textEdited.connect(self._view.pause_pvtimer)
        self._view.calcDisplayMap['h_target'].returnPressed.connect(self._view._calc_motors)

        self._view.calcDisplayMap['k_target'].textEdited.connect(self._view.pause_pvtimer)      
        self._view.calcDisplayMap['k_target'].returnPressed.connect(self._view._calc_motors)

        self._view.calcDisplayMap['l_target'].textEdited.connect(self._view.pause_pvtimer)
        self._view.calcDisplayMap['l_target'].returnPressed.connect(self._view._calc_motors)

        self._view.calcDisplayMap['phi_target'].textEdited.connect(self._view.pause_pvtimer)
        self._view.calcDisplayMap['phi_target'].returnPressed.connect(self._view._calc_motors)



        

def main():
    qview_app = QApplication([])
    qview_app_window = MonitorWindow()
    mm = Qview(view=qview_app_window)
    qview_app_window.show()
    sys.exit(qview_app.exec())


if __name__ == '__main__':
    main()
