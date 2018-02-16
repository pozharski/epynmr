from matplotlib.widgets import Button
from matplotlib.figure import Figure

class MPLGUI(Figure):
    def __init__(self, fig, *args, **kwds):
        self.fig = fig
        self.buttons = {}
    def add_button(self, label, fun):
        baxes = self.fig.add_axes([0.901, 0.85-0.05*len(self.buttons), 0.1, 0.048])
        self.buttons[label] = Button(baxes, label)
        self.buttons[label].on_clicked(fun)
    def add_key_button(self, label, key):
        baxes = self.fig.add_axes([0.901, 0.85-0.05*len(self.buttons), 0.1, 0.048])
        self.buttons[label] = Button(baxes, label)
        self.buttons[label].on_clicked(lambda event : self.fig.canvas.key_press_event(key))
    
