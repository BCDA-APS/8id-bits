"""
Spec Reader (Python) -- live plotting of SPEC data files.

A Python/Qt rewrite of the MATLAB ``specr`` GUI by Zhang Jiang (APS), keeping
the parts that matter for beamline use: scan selection, X/Y column choice, plot
styles, peak/COM/FWHM readout, and a scan monitor that polls the file while a
scan is being written.

Because ``ophyd_scan.py`` appends one complete row per point and closes the
file each time, the monitor can poll safely -- it never holds the file open and
never sees a partial row.

usage::

    python specr.py                                   # then File > Open
    python specr.py /path/to/experiment.spec          # open immediately
    python specr.py /path/to/experiment.spec --watch  # open and start monitoring
"""

import os
import sys

import numpy as np

import matplotlib

matplotlib.use("Qt5Agg")

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavToolbar
from matplotlib.figure import Figure

from PyQt5 import QtCore, QtGui, QtWidgets

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from specfile import SpecDataFile, scan_stats  # noqa: E402

# Matches the MATLAB specr cycling order so plots look familiar.
COLORS = ["b", "r", "g", "c", "m", "k"]
MARKERS = ["o", "s", "d", "^", "v", "<", ">"]
FACES = ["m", "b", "k", "r", "c", "g"]

PLOT_STYLES = ["linear", "logx", "logy", "logxy"]


def _icon(kind, on=False):
    """Draw the toolbar icons, mirroring the MATLAB specr play/eraser buttons.

    Painted rather than loaded from files so the package stays self-contained
    and the icons render identically on any workstation (the MATLAB original
    keeps these in a binary icons.mat).
    """
    pm = QtGui.QPixmap(22, 22)
    pm.fill(QtCore.Qt.transparent)
    p = QtGui.QPainter(pm)
    p.setRenderHint(QtGui.QPainter.Antialiasing, True)
    p.setPen(QtCore.Qt.NoPen)

    if kind == "live":
        if on:  # running -> offer "pause"
            p.setBrush(QtGui.QColor("#c0392b"))
            p.drawRect(6, 4, 4, 14)
            p.drawRect(12, 4, 4, 14)
        else:  # stopped -> offer "play"
            p.setBrush(QtGui.QColor("#27ae60"))
            p.drawPolygon(QtGui.QPolygon([
                QtCore.QPoint(6, 4), QtCore.QPoint(18, 11), QtCore.QPoint(6, 18)]))
    else:  # eraser
        body = QtGui.QColor("#e67e22") if on else QtGui.QColor("#9aa4ad")
        p.save()
        p.translate(11, 11)
        p.rotate(-35)
        p.translate(-11, -11)
        p.setBrush(body)
        p.drawRoundedRect(QtCore.QRectF(4.5, 5.5, 13, 11), 2, 2)
        p.setBrush(QtGui.QColor("#f4f6f7") if on else QtGui.QColor("#dfe4e8"))
        p.drawRect(QtCore.QRect(5, 12, 12, 4))
        p.restore()
    p.end()
    return QtGui.QIcon(pm)


class Settings:
    """Mirrors the ``settings`` struct in specr_CreateFcn."""

    def __init__(self):
        self.monitor_period = 2.0
        self.monitor_auto_period = True
        self.erase_mode = True  # show only the newest scan while monitoring
        self.show_grid = True
        self.show_legend = False


class SettingsDialog(QtWidgets.QDialog):
    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Spec Reader Settings")
        self.settings = settings

        form = QtWidgets.QFormLayout()
        self.period = QtWidgets.QDoubleSpinBox()
        self.period.setRange(0.1, 600.0)
        self.period.setDecimals(2)
        self.period.setSingleStep(0.5)
        self.period.setValue(settings.monitor_period)
        form.addRow("Monitor checking period (s):", self.period)

        self.auto = QtWidgets.QCheckBox("Automatically adjust period for large scans")
        self.auto.setChecked(settings.monitor_auto_period)
        form.addRow(self.auto)

        self.erase = QtWidgets.QCheckBox("Erase previous scans while monitoring")
        self.erase.setChecked(settings.erase_mode)
        form.addRow(self.erase)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def apply(self):
        self.settings.monitor_period = self.period.value()
        self.settings.monitor_auto_period = self.auto.isChecked()
        self.settings.erase_mode = self.erase.isChecked()


class TextDialog(QtWidgets.QDialog):
    """Read-only text viewer used for 'Show Scan' and 'Show Motor Positions'."""

    def __init__(self, title, text, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(800, 460)
        self.view = QtWidgets.QPlainTextEdit()
        self.view.setReadOnly(True)
        self.view.setPlainText(text)
        font = self.view.font()
        font.setFamily("monospace")
        self.view.setFont(font)
        close = QtWidgets.QPushButton("Close")
        close.clicked.connect(self.accept)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.view)
        row = QtWidgets.QHBoxLayout()
        row.addStretch(1)
        row.addWidget(close)
        layout.addLayout(row)

    def set_text(self, text):
        self.view.setPlainText(text)


class ScanSelectDialog(QtWidgets.QDialog):
    """Multi-select scan list, equivalent to MATLAB's listdlg in selectscan.m."""

    def __init__(self, headers, preselect, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Scan")
        self.resize(700, 500)
        self.list = QtWidgets.QListWidget()
        self.list.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        for h in headers:
            self.list.addItem(h)
        for i in preselect:
            if 0 <= i < self.list.count():
                self.list.item(i).setSelected(True)
        if preselect:
            self.list.scrollToItem(self.list.item(preselect[-1]))

        buttons = QtWidgets.QDialogButtonBox()
        buttons.addButton("Plot", QtWidgets.QDialogButtonBox.AcceptRole)
        buttons.addButton(QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.list)
        layout.addWidget(buttons)

    def selection(self):
        return sorted(i.row() for i in self.list.selectedIndexes())


class SpecrWindow(QtWidgets.QMainWindow):
    def __init__(self, path=None, watch=False):
        super().__init__()
        self.setWindowTitle("Spec Reader")
        self.resize(980, 620)

        self.settings = Settings()
        self.specfile = None
        self.selection = []
        self._last_labels = None

        self._build_ui()

        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._on_monitor_tick)

        if path:
            self.open_file(path)
            if watch:
                self.monitor_action.setChecked(True)

    # -- construction ------------------------------------------------------

    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QVBoxLayout(central)

        controls = QtWidgets.QHBoxLayout()
        self.select_btn = QtWidgets.QPushButton("Select Scan")
        self.select_btn.clicked.connect(self.select_scan)
        controls.addWidget(self.select_btn)

        controls.addWidget(QtWidgets.QLabel("X:"))
        self.x_combo = QtWidgets.QComboBox()
        self.x_combo.setMinimumWidth(190)
        self.x_combo.currentIndexChanged.connect(self.replot)
        controls.addWidget(self.x_combo)

        controls.addWidget(QtWidgets.QLabel("Y:"))
        self.y_combo = QtWidgets.QComboBox()
        self.y_combo.setMinimumWidth(190)
        self.y_combo.currentIndexChanged.connect(self.replot)
        controls.addWidget(self.y_combo)

        controls.addWidget(QtWidgets.QLabel("Style:"))
        self.style_combo = QtWidgets.QComboBox()
        self.style_combo.addItems(PLOT_STYLES)
        self.style_combo.currentIndexChanged.connect(self._apply_scale)
        controls.addWidget(self.style_combo)

        self.replot_btn = QtWidgets.QPushButton("Replot")
        self.replot_btn.clicked.connect(self.reload_and_plot)
        controls.addWidget(self.replot_btn)
        controls.addStretch(1)
        layout.addLayout(controls)

        self.figure = Figure(figsize=(7, 4.2), tight_layout=True)
        self.canvas = FigureCanvas(self.figure)
        self.axes = self.figure.add_subplot(111)
        self.axes.grid(True)
        layout.addWidget(NavToolbar(self.canvas, self))
        layout.addWidget(self.canvas, 1)

        self.status = self.statusBar()
        self.status.showMessage("No file loaded.")
        self.canvas.mpl_connect("motion_notify_event", self._on_mouse_move)

        self._build_menus()
        self._build_toolbar()
        self._sync_mode()

    def _build_menus(self):
        filemenu = self.menuBar().addMenu("&File")
        act = filemenu.addAction("&Open Spec File...")
        act.setShortcut("Ctrl+O")
        act.triggered.connect(self.open_dialog)
        act = filemenu.addAction("&Save Displayed Data...")
        act.setShortcut("Ctrl+S")
        act.triggered.connect(self.save_data)
        filemenu.addSeparator()
        act = filemenu.addAction("&Close")
        act.setShortcut("Ctrl+W")
        act.triggered.connect(self.close)

        tools = self.menuBar().addMenu("&Tools")
        act = tools.addAction("Sh&ow Current Scan")
        act.triggered.connect(self.show_scan_text)
        act = tools.addAction("Sho&w Current Motor Positions")
        act.triggered.connect(self.show_motor_positions)
        tools.addSeparator()
        act = tools.addAction("Previous Scan")
        act.setShortcut("Ctrl+Left")
        act.triggered.connect(lambda: self.step_scan(-1))
        act = tools.addAction("Next Scan")
        act.setShortcut("Ctrl+Right")
        act.triggered.connect(lambda: self.step_scan(+1))
        tools.addSeparator()
        act = tools.addAction("Invert &X Axis")
        act.triggered.connect(lambda: self._invert("x"))
        act = tools.addAction("Invert &Y Axis")
        act.triggered.connect(lambda: self._invert("y"))
        tools.addSeparator()
        act = tools.addAction("&Settings...")
        act.triggered.connect(self.edit_settings)

        view = self.menuBar().addMenu("&View")
        self.grid_action = view.addAction("Grid")
        self.grid_action.setCheckable(True)
        self.grid_action.setChecked(self.settings.show_grid)
        self.grid_action.triggered.connect(self.replot)
        self.legend_action = view.addAction("Legend")
        self.legend_action.setCheckable(True)
        self.legend_action.setChecked(self.settings.show_legend)
        self.legend_action.triggered.connect(self.replot)
        view.addSeparator()
        self.erase_action = view.addAction("&Erase mode")
        self.erase_action.setCheckable(True)
        self.erase_action.setChecked(self.settings.erase_mode)
        self.erase_action.toggled.connect(self._on_erase_toggled)
        self.monitor_action = view.addAction("Live mode (Scan &Monitor)")
        self.monitor_action.setCheckable(True)
        self.monitor_action.setShortcut("Ctrl+M")
        self.monitor_action.toggled.connect(self.toggle_monitor)

        helpmenu = self.menuBar().addMenu("&Help")
        helpmenu.addAction("About").triggered.connect(
            lambda: QtWidgets.QMessageBox.about(
                self,
                "About Spec Reader",
                "Spec Reader (Python)\n\n"
                "Live SPEC-file viewer for APS 8-ID.\n"
                "Python/Qt rewrite of the MATLAB 'specr' by Zhang Jiang.",
            )
        )

    def _build_toolbar(self):
        """Play / eraser toggles, as on the MATLAB specr toolbar.

        The same QActions back the View menu, so menu and toolbar can never
        disagree about the current mode.
        """
        bar = self.addToolBar("Mode")
        bar.setIconSize(QtCore.QSize(22, 22))
        bar.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        bar.setMovable(False)

        spacer = QtWidgets.QWidget()
        spacer.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                             QtWidgets.QSizePolicy.Preferred)
        bar.addWidget(spacer)          # push the toggles to the right, as in MATLAB

        self.mode_label = QtWidgets.QLabel()
        self.mode_label.setStyleSheet("QLabel { padding-right: 10px; }")
        bar.addWidget(self.mode_label)
        bar.addAction(self.erase_action)
        bar.addAction(self.monitor_action)

    def _on_erase_toggled(self, on):
        self.settings.erase_mode = on
        self._sync_mode()
        if on and len(self.selection) > 1:
            self.selection = self.selection[-1:]   # collapse to the newest
            self.replot()

    def _sync_mode(self):
        """Keep icons, tooltips and the mode caption in step with the toggles."""
        live = self.monitor_action.isChecked()
        erase = self.erase_action.isChecked()

        self.monitor_action.setIcon(_icon("live", live))
        self.monitor_action.setText("Live" if live else "Display")
        self.monitor_action.setToolTip(
            "Live mode: following the file as scans are written — click to stop "
            "and browse (Ctrl+M)" if live else
            "Display mode: browse and overlay finished scans — click to follow "
            "the file live (Ctrl+M)")

        self.erase_action.setIcon(_icon("erase", erase))
        self.erase_action.setText("Erase")
        self.erase_action.setToolTip(
            "Erase on: only the newest scan is shown" if erase else
            "Erase off: completed scans stay on the plot and accumulate")

        if live and erase:
            caption = "live · current scan only"
        elif live:
            caption = "live · accumulating"
        elif erase:
            caption = "display · single scan"
        else:
            caption = "display · overlay"
        self.mode_label.setText(caption)

    # -- file handling -----------------------------------------------------

    def open_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Select Spec File", "", "SPEC files (*.spec *.spe *.dat);;All files (*)"
        )
        if path:
            self.open_file(path)

    def open_file(self, path):
        specfile = SpecDataFile(path)
        specfile.refresh()
        if len(specfile) == 0:
            QtWidgets.QMessageBox.critical(
                self, "File Open Error", "Invalid file or no scan in this file."
            )
            return
        self.specfile = specfile
        self._last_labels = None
        self.selection = [len(specfile) - 1]
        self.setWindowTitle(f"Spec Reader - {os.path.basename(path)}")
        self._sync_columns(force=True)
        self.replot()

    def reload_and_plot(self):
        if self.specfile:
            self.specfile.refresh()
            self._sync_columns()
            self.replot()

    # -- scan / column selection ------------------------------------------

    def select_scan(self):
        if not self.specfile:
            return
        self.specfile.refresh()
        dialog = ScanSelectDialog(self.specfile.scan_headers(), self.selection, self)
        if dialog.exec_() != QtWidgets.QDialog.Accepted:
            return
        chosen = dialog.selection()
        if not chosen:
            return
        labels = {tuple(self.specfile[i].labels) for i in chosen}
        if len(labels) > 1:
            QtWidgets.QMessageBox.critical(
                self,
                "Select Scan Error",
                "Multi-selections have to be scans of the same type.",
            )
            return
        # Erase off means "keep what is already plotted" -- union with the
        # current selection instead of replacing it.
        if self.erase_action.isChecked():
            self.selection = chosen
        else:
            keep = [i for i in self.selection
                    if i < len(self.specfile)
                    and self.specfile[i].labels == self.specfile[chosen[0]].labels]
            self.selection = sorted(set(keep) | set(chosen))
        self._sync_columns(force=True)
        self.replot()

    def step_scan(self, delta):
        if not self.specfile or not self.selection:
            return
        self.specfile.refresh()
        i = self.selection[-1] if delta > 0 else self.selection[0]
        i = max(0, min(len(self.specfile) - 1, i + delta))
        if self.erase_action.isChecked() or i in self.selection:
            self.selection = [i]
        elif self.specfile[i].labels == self.specfile[self.selection[-1]].labels:
            self.selection = sorted(set(self.selection) | {i})
        else:
            self.selection = [i]
        self._sync_columns(force=True)
        self.replot()

    def _sync_columns(self, force=False):
        """Repopulate X/Y combos, preserving the user's choice when possible.

        Defaults follow the MATLAB reader: X = first column, Y = last column.
        """
        scan = self._primary_scan()
        if scan is None or not scan.labels:
            return
        labels = list(scan.labels)
        if not force and labels == self._last_labels:
            return

        prev_x = self.x_combo.currentText()
        prev_y = self.y_combo.currentText()
        for combo, previous, default in (
            (self.x_combo, prev_x, 0),
            (self.y_combo, prev_y, len(labels) - 1),
        ):
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(labels)
            combo.setCurrentIndex(
                labels.index(previous) if (previous in labels and not force) else default
            )
            combo.blockSignals(False)
        self._last_labels = labels

    def _primary_scan(self):
        if not self.specfile or not self.selection:
            return None
        idx = self.selection[-1]
        if idx >= len(self.specfile):
            return None
        return self.specfile[idx]

    # -- plotting ----------------------------------------------------------

    def replot(self):
        if not self.specfile or not self.selection:
            return
        xlabel = self.x_combo.currentText()
        ylabel = self.y_combo.currentText()
        if not xlabel or not ylabel:
            return

        self.axes.clear()
        plotted = []
        for n, idx in enumerate(self.selection):
            if idx >= len(self.specfile):
                continue
            scan = self.specfile[idx]
            if xlabel not in scan.labels or ylabel not in scan.labels:
                continue
            x, y = scan.column(xlabel), scan.column(ylabel)
            if x.size == 0:
                continue
            size = min(x.size, y.size)
            self.axes.plot(
                x[:size],
                y[:size],
                color=COLORS[n % len(COLORS)],
                marker=MARKERS[n % len(MARKERS)],
                markerfacecolor=FACES[n % len(FACES)],
                markersize=5,
                linestyle="-",
                label=f"Scan {scan.number}",
            )
            plotted.append(scan)

        self.axes.set_ylabel(ylabel)
        self.axes.grid(self.grid_action.isChecked())
        if self.legend_action.isChecked() and plotted:
            self.axes.legend(loc="best", fontsize=8)

        if plotted:
            last = plotted[-1]
            x, y = last.column(xlabel), last.column(ylabel)
            size = min(x.size, y.size)
            peak_x, peak_y, com, fwhm, center = scan_stats(x[:size], y[:size])
            self.axes.set_xlabel(
                f"{xlabel}\n"
                f"Peak {peak_y:g} @ {peak_x:g},   COM {com:g},   "
                f"FWHM {fwhm:g} @ {center:g}"
            )
            numbers = ", ".join(str(s.number) for s in plotted)
            title = f"{os.path.basename(self.specfile.path)},   Scan {numbers}"
            if last.date:
                title += f",   {last.date}"
            self.axes.set_title(f"{title}\n{last.command}", fontsize=9)
            self.status.showMessage(
                f"Scan {last.number}: {last.npoints} points"
                + (f"   [{last.exit_status}]" if last.exit_status else "   [running]")
            )
        else:
            self.axes.set_xlabel(xlabel)

        self._apply_scale(redraw=False)
        self.canvas.draw_idle()

    def _apply_scale(self, redraw=True):
        style = self.style_combo.currentText()
        self.axes.set_xscale("log" if "logx" in style or style == "logxy" else "linear")
        self.axes.set_yscale("log" if "logy" in style or style == "logxy" else "linear")
        if redraw:
            self.canvas.draw_idle()

    def _invert(self, axis):
        for line in self.axes.get_lines():
            if axis == "x":
                line.set_xdata(-np.asarray(line.get_xdata()))
            else:
                y = np.asarray(line.get_ydata())
                line.set_ydata(np.nanmax(y) + np.nanmin(y) - y)
        self.axes.relim()
        self.axes.autoscale_view()
        self.canvas.draw_idle()

    def _on_mouse_move(self, event):
        if event.inaxes:
            self.status.showMessage(f"x = {event.xdata:g}    y = {event.ydata:g}")

    # -- live monitor ------------------------------------------------------

    def toggle_monitor(self, on):
        if on and not self.specfile:
            self.monitor_action.setChecked(False)
            return
        # Scan selection is a display-mode action: while following the file the
        # monitor owns which scan is shown (as in the MATLAB GUI).
        self.select_btn.setEnabled(not on)
        self._sync_mode()
        if on:
            self.specfile.refresh()
            if self.erase_action.isChecked() and len(self.specfile):
                self.selection = [len(self.specfile) - 1]
                self._sync_columns()
                self.replot()
            self.timer.start(int(self.settings.monitor_period * 1000))
            self.status.showMessage("Live — following the file.")
        else:
            self.timer.stop()
            self.status.showMessage("Display mode — monitor stopped.")

    def _on_monitor_tick(self):
        if not self.specfile:
            return
        started = QtCore.QElapsedTimer()
        started.start()

        changed = self.specfile.refresh()
        if changed:
            newest = len(self.specfile) - 1
            if self.erase_action.isChecked():
                self.selection = [newest]
            elif newest not in self.selection:
                # only accumulate scans that share the column layout
                if self.specfile[newest].labels == self.specfile[self.selection[-1]].labels:
                    self.selection.append(newest)
                else:
                    self.selection = [newest]
            self._sync_columns()
            self.replot()

        # Same idea as the MATLAB auto-period: if a refresh takes longer than
        # the poll interval, back the interval off so we never queue up.
        if self.settings.monitor_auto_period:
            elapsed = started.elapsed() / 1000.0
            wanted = max(self.settings.monitor_period, elapsed + 0.1)
            if abs(wanted * 1000 - self.timer.interval()) > 100:
                self.timer.setInterval(int(wanted * 1000))

    # -- tools -------------------------------------------------------------

    def edit_settings(self):
        dialog = SettingsDialog(self.settings, self)
        if dialog.exec_() == QtWidgets.QDialog.Accepted:
            dialog.apply()
            self.erase_action.setChecked(self.settings.erase_mode)
            if self.timer.isActive():
                self.timer.setInterval(int(self.settings.monitor_period * 1000))

    def show_scan_text(self):
        scan = self._primary_scan()
        if scan is None:
            return
        lines = [scan.header, f"#D {scan.date}"]
        lines += [f"#C {c}" for c in scan.comments]
        lines += [f"#MD {k} = {v}" for k, v in scan.meta.items()]
        lines.append("#N " + str(len(scan.labels)))
        lines.append("#L " + "  ".join(scan.labels))
        lines += [" ".join(f"{v:g}" for v in row) for row in scan.rows]
        TextDialog(f"Scan {scan.number}", "\n".join(lines), self).exec_()

    def show_motor_positions(self):
        scan = self._primary_scan()
        if scan is None:
            return
        table = scan.motor_table()
        if not table:
            text = "No #O/#P motor information in this file."
        else:
            width = max(len(n) for n, _ in table)
            text = "\n".join(f"{n:<{width}}  {p:g}" for n, p in table)
        TextDialog(f"Motor positions - scan {scan.number}", text, self).exec_()

    def save_data(self):
        scan = self._primary_scan()
        if scan is None:
            return
        xlabel, ylabel = self.x_combo.currentText(), self.y_combo.currentText()
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save Data As", f"scan{scan.number}.dat", "Data (*.dat);;All files (*)"
        )
        if not path:
            return
        x, y = scan.column(xlabel), scan.column(ylabel)
        size = min(x.size, y.size)
        with open(path, "w") as f:
            f.write(f"# {os.path.basename(self.specfile.path)}  scan {scan.number}\n")
            f.write(f"# {xlabel}\t{ylabel}\n")
            for xi, yi in zip(x[:size], y[:size]):
                f.write(f"{xi:.12g}\t{yi:.12g}\n")
        self.status.showMessage(f"Saved {size} points to {path}")

    def closeEvent(self, event):
        self.timer.stop()
        event.accept()


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    watch = "--watch" in argv
    argv = [a for a in argv if not a.startswith("--")]
    path = argv[0] if argv else None

    app = QtWidgets.QApplication(sys.argv)
    window = SpecrWindow(path=path, watch=watch)
    window.show()
    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
