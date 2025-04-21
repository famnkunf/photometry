import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from astropy.io import fits
from astropy.visualization import ImageNormalize, LinearStretch
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.style as mplstyle
mplstyle.use('fast')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import patches
# from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
import numpy as np
import regex as re
import time
import pandas as pd
import skimage

class TopWindow(tk.Tk):
    def __init__(self):
        super().__init__()
        ttk.Style(self).theme_use("xpnative")
        self.title("Photometry Tool")
        self.geometry("800x50")
        self.attributes("-topmost", True)
        self.create_widget()
        self.display_windows = []
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.information_window = None
        self.calculator_window = None
        self.graph_window = None
        self.histogram_window = None
        self.objects_window = None
        
    def create_widget(self):
        self.control_frame = ttk.Frame(self)
        self.control_frame.pack(side=tk.TOP, fill=tk.X)
        
        self.open_data_button = ttk.Button(self.control_frame, text="Open FITS File", command=self.open_fits_file)
        self.open_data_button.pack(side=tk.LEFT)
        
        self.information_button = ttk.Button(self.control_frame, text="Information", command=self.show_information)
        self.information_button.pack(side=tk.LEFT)
        
        self.calculator_button = ttk.Button(self.control_frame, text="Calculator", command=self.open_calculator)
        self.calculator_button.pack(side=tk.LEFT)
        
        self.graph_button = ttk.Button(self.control_frame, text="Graph", command=self.show_graph)
        self.graph_button.pack(side=tk.LEFT)
        
        self.histogram_button = ttk.Button(self.control_frame, text="Histogram", command=self.show_histogram)
        self.histogram_button.pack(side=tk.LEFT)
        
        self.objects_button = ttk.Button(self.control_frame, text="Objects", command=self.show_objects)
        self.objects_button.pack(side=tk.LEFT)
        
        
    def open_calculator(self):
        if self.calculator_window:
            self.calculator_window.deiconify()
            self.calculator_window.focus_force()
        else:
            self.calculator_window = Calculator(self)
            self.calculator_window.protocol("WM_DELETE_WINDOW", self.calculator_window.close)
            # self.calculator_window.mainloop()    
    
    def open_fits_file(self):
        file_paths = filedialog.askopenfilenames(filetypes=[("FITS files", "*.fits *.fit")])
        if file_paths:
            for file_path in file_paths:
                try:
                    image, header = fits.getdata(file_path, header=True)
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to open FITS file: {e}")
                self.open_display_window(image, header, file_path)

    def open_display_window(self, image, header, file_path):
        window = DisplayWindow(image, header, file_path, self)
        self.display_windows.append(window)
        self.display_windows[-1].protocol("WM_DELETE_WINDOW", self.display_windows[-1].close)

        if self.objects_window:
            self.objects_window.add_window(window)
        
        # self.display_windows[-1].mainloop()
        
    def close(self):
        if self.information_window:
            self.information_window.close()
        for window in self.display_windows:
            window.close()
        if self.calculator_window:
            self.calculator_window.close()
        if self.graph_window:
            self.graph_window.close()
        if self.histogram_window:
            self.histogram_window.close()
        if self.objects_window:
            self.objects_window.close()
        self.quit()
        self.destroy()
        
    def show_information(self):
        if self.information_window:
            self.information_window.deiconify()
            self.information_window.focus_force()
        else:
            self.information_window = InformationWindow(self)
            self.information_window.protocol("WM_DELETE_WINDOW", self.information_window.close)

    def show_graph(self):
        if self.graph_window:
            self.graph_window.deiconify()
            self.graph_window.focus_force()
        else:
            self.graph_window = Graph(self)
            self.graph_window.protocol("WM_DELETE_WINDOW", self.graph_window.close)

    def show_histogram(self):
        if self.histogram_window:
            self.histogram_window.deiconify()
            self.histogram_window.focus_force()
        else:
            self.histogram_window = Histogram(self)
            self.histogram_window.protocol("WM_DELETE_WINDOW", self.histogram_window.close)
            
    def show_objects(self):
        if self.objects_window:
            self.objects_window.deiconify()
            self.objects_window.focus_force()
        else:
            self.objects_window = ObjectsWindow(self)
            self.objects_window.protocol("WM_DELETE_WINDOW", self.objects_window.close)
class DisplayWindow(tk.Toplevel):
    def __init__(self, image, header, title, parent:TopWindow=None):
        super().__init__()
        self.title(title)
        self.geometry("800x600")
        self.image = image
        self.header = header
        self.parent = parent
        self.header_window = None
        self.norm = ImageNormalize(self.image, stretch=LinearStretch())
        self.pan_start = None
        self.adding_aperture = False
        self.aperture_window = None
        self.added_apertures = []
        self.aperture = ()
        self.aperture_prev_params = ()
        
        self.create_widget()

    def create_widget(self):
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw_idle()
        self.display_image()
        self.drawing_type = False
        self.coords1 = None
        self.coords2 = None
        self.annotate = None
        self.data = None
        self.last_update_time = time.time()
        self.target_interval = 1/10
        
        self.control_frame1 = ttk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = ttk.Frame(self)
        self.control_frame2.pack(side=tk.TOP, fill=tk.X)
        
        self.header_button = ttk.Button(self.control_frame1, text="Show Header", command=self.show_header)
        self.header_button.pack(side=tk.LEFT)
        
        self.save_button = ttk.Button(self.control_frame1, text="Save", command=self.save)
        self.save_button.pack(side=tk.LEFT)
        
        self.add_aperture_button = ttk.Button(self.control_frame1, text="Add Aperture", command=self.add_aperture)
        self.add_aperture_button.pack(side=tk.LEFT)
        
        self.positionX_label = ttk.Label(self.control_frame2, text="X:")
        self.positionX_label.pack(side=tk.LEFT)
        self.positionY_label = ttk.Label(self.control_frame2, text="Y:")
        self.positionY_label.pack(side=tk.LEFT)

        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('button_press_event', self.start_pan)
        self.canvas.mpl_connect('button_release_event', self.stop_pan)
        self.canvas.mpl_connect('motion_notify_event', self.on_pan)
        
        self.bind("<Button-1>", self.on_focus)
        
    def show_header(self):
        if self.header_window is None:
            self.header_window = HeaderWindow(self.header, self)
            self.header_window.protocol("WM_DELETE_WINDOW", self.header_window.close)
        else:
            self.header_window.deiconify()
            self.header_window.focus_force()

    def display_image(self):
        self.ax.imshow(self.image, cmap='gray', norm=self.norm)
        self.ax.set_title("FITS Image")
        self.ax.set_xlim(0, self.image.shape[1])
        self.ax.set_ylim(0, self.image.shape[0])
        self.ax.set_aspect('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.canvas.draw_idle()
        
    def add_aperture(self):
        self.adding_aperture = not self.adding_aperture
        if self.add_aperture:
            self.add_aperture_button.state(['pressed'])
            if len(self.aperture_prev_params) > 0:
                self.aperture_window = Aperture(self, *self.aperture_prev_params)
            else:
                self.aperture_window = Aperture(self)
            self.aperture_window.protocol("WM_DELETE_WINDOW", self.aperture_window.close)
        else:
            self.add_aperture_button.state(['!pressed'])
            self.aperture_window.close()
        
    def draw_aperture(self, x, y):
        a = float(self.aperture_window.aperture_major_var.get())
        b = float(self.aperture_window.aperture_minor_var.get())
        gap = float(self.aperture_window.gap.get())
        angle = float(self.aperture_window.aperture_angle.get())
        background = float(self.aperture_window.background.get())
        self.aperture_prev_params = (a, b, angle, gap, background)
        p = patches.Ellipse((x, y), a, b, angle=angle, color='red', fill=False, lw=1)
        g = patches.Ellipse((x, y), a+gap, b+gap, angle=angle, color='red', fill=False, lw=1, linestyle='dashed')
        bg = patches.Ellipse((x, y), a+gap+background, b+gap+background, angle=angle, color='blue', fill=False, lw=1, linestyle='dashed')
        if len(self.aperture) > 0:
            [i.remove() for i in self.aperture]
        self.aperture = (p, g, bg)
        [self.ax.add_patch(i) for i in self.aperture]
        self.canvas.draw_idle()     
        
    def on_mouse_move(self, event):
        if event.xdata is not None and event.ydata is not None:
            x, y = int(event.xdata), int(event.ydata)
            try:
                pixel_value = self.image[y, x]
            except IndexError:
                pixel_value = None
            if self.parent.information_window:
                self.parent.information_window.pixel_value_label.config(text=f"Pixel Value: {pixel_value}")
            self.positionX_label.config(text=f"X: {x}")
            self.positionY_label.config(text=f"Y: {y}")
            if self.adding_aperture:
                self.draw_aperture(x, y)             
            else:
                if self.drawing_type == "Line":
                    if self.coords1 is not None and self.coords2 is None:
                        self.coords2 = (event.xdata, event.ydata)
                        temp = self.coords1
                        self.drawline()
                        self.coords1 = temp
                elif self.drawing_type == "Horizontal Box":
                    if self.coords1 is not None and self.coords2 is None:
                        self.coords2 = (event.xdata, event.ydata)
                        temp = self.coords1
                        self.drawrect()
                        self.coords1 = temp
                elif self.drawing_type == "Area":
                    if self.coords1 is not None and self.coords2 is None:
                        self.coords2 = (event.xdata, event.ydata)
                        temp = self.coords1
                        self.drawarea()
                        self.coords1 = temp
        
    def on_scroll(self, event):
        base_scale = 1.2
        if event.button == 'up':
            scale_factor = 1 / base_scale
        elif event.button == 'down':
            scale_factor = base_scale
        else:
            return
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        xdata = event.xdata if event.xdata is not None else 0
        ydata = event.ydata if event.ydata is not None else 0
        new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
        new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor
        relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
        rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])
        new_xlim = [xdata - new_width * (1 - relx), xdata + new_width * relx]
        new_ylim = [ydata - new_height * (1 - rely), ydata + new_height * rely]
        if new_xlim[0] < 0:
            new_xlim[0] = 0
        if new_xlim[1] > self.image.shape[1]:
            new_xlim[1] = self.image.shape[1]
        if new_ylim[0] < 0:
            new_ylim[0] = 0
        if new_ylim[1] > self.image.shape[0]:
            new_ylim[1] = self.image.shape[0]
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        self.ax.set_aspect('equal')
        self.canvas.draw_idle()

    def on_focus(self, event):
        if self.parent.graph_window:
            self.drawing_type = self.parent.graph_window.graph_type.get()
        else:
            if self.annotate:
                self.annotate.remove()
                self.annotate = None
            self.drawing_type = False

    def start_pan(self, event):
        if event.button == 1:
            if not self.drawing_type:
                self.pan_start = (event.x, event.y)
        
            
    def stop_pan(self, event):
        if event.button == 1:
            if self.adding_aperture:
                if len(self.aperture) > 0:
                    x, y = self.aperture_window.get_max(self.aperture[0], self.image)
                    for i in self.aperture:
                        i.remove()
                        i.set_color('yellow')
                        i.set_center((x, y))
                        self.ax.add_patch(i)
                    self.added_apertures.append(self.aperture)
                    if self.parent.objects_window:
                        self.parent.objects_window.add(self.aperture, self)
                    self.aperture = ()
                    self.aperture_window.close()
                    self.add_aperture_button.state(['!pressed'])
                    self.adding_aperture = False
                    self.canvas.draw_idle()
            else:
                if not self.drawing_type:
                    self.pan_start = None
                else:
                    if self.coords1 is None and event.xdata is not None and event.ydata is not None:
                        self.coords1 = (event.xdata, event.ydata)
                    elif self.coords2 is None and event.xdata is not None and event.ydata is not None:
                        self.coords2 = (event.xdata, event.ydata)
                        if self.drawing_type == "Horizontal Box":
                            self.drawrect()
                        elif self.drawing_type == "Line":
                            self.drawline()
                        elif self.drawing_type == "Area":
                            self.drawarea()
    
    def drawarea(self):
        current_time = time.time()
        if self.coords1 and self.coords2:
            if self.parent.graph_window:
                if current_time - self.last_update_time > self.target_interval:
                    self.data = self.get_area(int(self.coords1[0]), int(self.coords1[1]), int(self.coords2[0]), int(self.coords2[1]))
                    self.parent.graph_window.update_graph(self.data)
                self.last_update_time = current_time
            x1 = int(self.coords1[0])
            y1 = int(self.coords1[1])
            x2 = int(self.coords2[0])
            y2 = int(self.coords2[1])
            if self.annotate is not None:
                self.annotate.remove()
            self.annotate = patches.Rectangle((x1, y1), x2-x1, y2-y1, linewidth=1, edgecolor='red', facecolor='none')
            self.ax.add_patch(self.annotate)
            self.canvas.draw_idle()
            self.coords1 = None
            self.coords2 = None
    
    def drawrect(self):
        current_time = time.time()
        
        if self.coords1 and self.coords2:
            if self.parent.graph_window:
                if current_time - self.last_update_time > self.target_interval:
                    self.data = self.get_pixels_inside_box(int(self.coords1[0]), int(self.coords1[1]), int(self.coords2[0]), int(self.coords2[1]))
                    self.parent.graph_window.update_graph(self.data)
                self.last_update_time = current_time
            x1 = int(self.coords1[0])
            y1 = int(self.coords1[1])
            x2 = int(self.coords2[0])
            y2 = int(self.coords2[1])
            if self.annotate:
                self.annotate.remove()
            self.annotate = patches.Rectangle((x1, y1), x2-x1, y2-y1, linewidth=1, edgecolor='red', facecolor='none')
            self.ax.add_patch(self.annotate)
            self.canvas.draw_idle()
            self.coords1 = None
            self.coords2 = None

    def drawline(self):
        current_time = time.time()
        
        if self.coords1 and self.coords2:
            if self.parent.graph_window:
                if current_time - self.last_update_time > self.target_interval:
                    self.data = self.get_pixels_along_line(int(self.coords1[0]), int(self.coords1[1]), int(self.coords2[0]), int(self.coords2[1]))
                    self.parent.graph_window.update_graph(self.data)
                self.last_update_time = current_time    
            x = [self.coords1[0], self.coords2[0]]
            y = [self.coords1[1], self.coords2[1]]
            if self.annotate is not None:
                self.annotate.remove()
            self.annotate = self.ax.plot(x, y, color='red', linewidth=1.5)[0]
            self.canvas.draw_idle()
            self.coords1 = None
            self.coords2 = None
        
    def on_pan(self, event):
        if self.pan_start is None or event.button != 1:
            return
        dx = event.x - self.pan_start[0]
        dy = event.y - self.pan_start[1]
        self.pan_start = (event.x, event.y)
        cur_xlim = self.ax.get_xlim()
        cur_ylim = self.ax.get_ylim()
        scale_x = (cur_xlim[1] - cur_xlim[0]) / self.canvas.get_tk_widget().winfo_width()
        scale_y = (cur_ylim[1] - cur_ylim[0]) / self.canvas.get_tk_widget().winfo_height()
        new_xlim = [cur_xlim[0] - dx * scale_x, cur_xlim[1] - dx * scale_x]
        new_ylim = [cur_ylim[0] - dy * scale_y, cur_ylim[1] - dy * scale_y]
        if new_xlim[0] < 0:
            new_xlim[0] = 0
            new_xlim[1] = cur_xlim[1] - cur_xlim[0]
        if new_xlim[1] > self.image.shape[1]:
            new_xlim[1] = self.image.shape[1]
            new_xlim[0] = self.image.shape[1] - (cur_xlim[1] - cur_xlim[0])
        if new_ylim[0] < 0:
            new_ylim[0] = 0
            new_ylim[1] = cur_ylim[1] - cur_ylim[0]
        if new_ylim[1] > self.image.shape[0]:
            new_ylim[1] = self.image.shape[0]
            new_ylim[0] = self.image.shape[0] - (cur_ylim[1] - cur_ylim[0])
        self.ax.set_xlim(new_xlim)
        self.ax.set_ylim(new_ylim)
        self.canvas.draw_idle()
        
    def get_pixels_along_line(self, x1, y1, x2, y2):
        pixels = []
        dx = abs(x2 - x1)
        dy = abs(y2 - y1)
        sx = 1 if x1 < x2 else -1
        sy = 1 if y1 < y2 else -1
        err = dx - dy
        cx, cy = x1, y1

        while True:
            if 0 <= cx < self.image.shape[1] and 0 <= cy < self.image.shape[0]:
                pixels.append(self.image[cy, cx])
            if cx == x2 and cy == y2:
                break
            e2 = 2 * err
            if e2 > -dy:
                err -= dy
                cx += sx
            if e2 < dx:
                err += dx
                cy += sy
        return pixels

    def get_pixels_inside_box(self, x1, y1, x2, y2):
        pixels = []
        x1, x2 = sorted([x1, x2])
        y1, y2 = sorted([y1, y2])
        for x in range(x1, x2):
            temp = []
            for y in range(y1, y2):
                if 0 <= x < self.image.shape[1] and 0 <= y < self.image.shape[0]:
                    temp.append(self.image[y, x])
            if len(temp) > 0:
                pixels.append(np.mean(temp))
        return pixels
    
    def get_area(self, x1, y1, x2, y2):
        area = []
        x1, x2 = sorted([x1, x2])
        y1, y2 = sorted([y1, y2])
        area = self.image[y1:y2, x1:x2]
        return area
        
    def save(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".fit", filetypes=[("FITS Files", "*.fits, *.fit")])
        if file_path:
            fits.writeto(file_path, self.image, header=self.header, overwrite=True)
            messagebox.showinfo("Saved", f"Image saved to {file_path}")
        
    def close(self):
        if self.header_window:
            self.header_window.close()
        if self.aperture_window:
            self.aperture_window.close()
        self.parent.display_windows.remove(self)
        self.destroy()

class InformationWindow(tk.Toplevel):
    def __init__(self, parent: TopWindow=None):
        super().__init__()
        self.parent = parent
        self.title("Information")
        self.geometry("400x300")
        
        self.create_widget()
        
    def create_widget(self):
        self.pixel_value_label = tk.Label(self, text="Pixel Value:")
        self.pixel_value_label.pack(side=tk.TOP, fill=tk.X)
        self.attributes("-topmost", True)
        
    def close(self):
        self.parent.information_window = None
        self.destroy()

class HeaderWindow(tk.Toplevel):
    def __init__(self, header, parent: DisplayWindow=None):
        super().__init__()
        self.parent = parent
        self.title(self.parent.title())
        self.geometry("400x300")
        
        self.create_widget(header)
    def create_widget(self, header):
        self.header_text = tk.Text(self)
        self.header_text.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        
        for key in header:
            self.header_text.insert(tk.END, f"{key}: {header[key]}\n")
        
        self.close_button = tk.Button(self, text="Close", command=self.close)
        self.close_button.pack(side=tk.BOTTOM)

    def close(self):
        self.parent.header_window = None
        self.destroy()

class Calculator(tk.Toplevel):
    def __init__(self, parent: TopWindow=None):
        super().__init__()
        self.parent = parent
        self.title("Calculator")
        self.geometry("400x300")
        self.attributes("-topmost", True)
        
        self.create_widget()
        
        self.variables = {}
        self.n_variables = 0

    def create_widget(self):
        self.control_frame1 = ttk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = ttk.Frame(self)
        self.control_frame2.pack(side=tk.TOP, fill=tk.X)
        self.control_frame3 = ttk.Frame(self)
        self.control_frame3.pack(side=tk.TOP, fill=tk.X)
        self.canvas = tk.Canvas(self)
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.yscrollbar4 = ttk.Scrollbar(self.canvas, orient=tk.VERTICAL, command=self.canvas.yview)
        self.yscrollbar4.pack(side=tk.RIGHT, fill=tk.Y)
        self.canvas.config(yscrollcommand=self.yscrollbar4.set)
        self.control_frame4 = tk.Frame(self.canvas)
        self.control_frame4.pack(side=tk.LEFT, fill=tk.X)

        self.control_frame4_id = self.canvas.create_window((0, 0), window=self.control_frame4, anchor=tk.NW, width=self.winfo_width()-self.yscrollbar4.winfo_width()-5)
        self.canvas.bind("<Configure>", self.on_canvas_resize)
        self.control_frame4.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox("all"))

        
        self.add_new_variable_button = ttk.Button(self.control_frame1, text="Add Variable", command=self.add_new_variable)
        self.add_new_variable_button.pack(side=tk.TOP, fill=tk.X)
        self.formular_label = ttk.Label(self.control_frame2, text="Formula:")
        self.formular_label.pack(side=tk.TOP, fill=tk.X)
        self.formular_entry = tk.Text(self.control_frame2, height=5)
        self.formular_entry.pack(side=tk.TOP, fill=tk.X)
        self.calculate_button = ttk.Button(self.control_frame3, text="Calculate", command=self.calculate)
        self.calculate_button.pack(side=tk.TOP, fill=tk.X)
        
    def on_canvas_resize(self, event):
        self.canvas.itemconfig(self.control_frame4_id, width=event.width-self.yscrollbar4.winfo_width()-5)
        self.canvas.config(scrollregion=self.canvas.bbox("all"))
    
    def add_new_variable(self):
        available_image = [window.title() for window in self.parent.display_windows]
        
        control_frame = ttk.Frame(self.control_frame4)
        control_frame.pack(side=tk.TOP, fill=tk.X)
        select_variable_label = ttk.Label(control_frame, text=f"Variable ${self.n_variables}")
        select_variable_label.pack(side=tk.LEFT)

        select_variable = ttk.Combobox(control_frame, values=available_image)
        select_variable.pack(side=tk.LEFT, fill=tk.X, expand=True)
        select_variable.bind("<ButtonPress-1>", self.update_selections)
        select_variable.bind("<<ComboboxSelected>>", self.update_variable_name)
        
        self.variables[f'${self.n_variables}'] = {}
        self.variables[f'${self.n_variables}']['entry'] = select_variable
        self.n_variables += 1
        
        self.control_frame4.update_idletasks()
        self.canvas.config(scrollregion=self.canvas.bbox("all"))
        
    def update_selections(self, event):
        available_image = [window.title() for window in self.parent.display_windows]
        event.widget['values'] = available_image
        
    def update_variable_name(self, event):
        available_image = [window.title() for window in self.parent.display_windows]
        for i in event.widget.master.winfo_children():
            if isinstance(i, ttk.Label):
                r = re.compile(r'Variable \$(\d+)')
                match = r.search(i.cget("text"))
        if match:
            self.variables[f'${match.group(1)}']['display'] = self.parent.display_windows[available_image.index(event.widget.get())]
    def calculate(self):
        formular = self.formular_entry.get("1.0", tk.END)
        replace_dict = {
            "#sqrt(": "np.sqrt(",
            "#log(": "np.log(",
            "#exp(": "np.exp(",
            "#sin(": "np.sin(",
            "#cos(": "np.cos(",
            "#tan(": "np.tan(",
            "#asin(": "np.arcsin(",
            "#acos(": "np.arccos(",
            "#atan(": "np.arctan(",
        }
        # for key in self.variables.keys():
        #     if key in formular:
        #         formular = formular.replace(key, f"self.variables['{key}']['display'].image")
        for key in replace_dict.keys():
            if key in formular:
                formular = formular.replace(key, replace_dict[key])
        output = re.compile(r'!#(.*)$', re.MULTILINE)
        var = re.compile(r'\$\((\w+)\)')
        formular = var.sub(self.format_var, formular)
        formular = output.sub(self.format, formular)    
        try:
            result = exec(formular)
        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {e, formular}")
            return

    def format_var(self, match):
        return "eval(\"self.variables[('$' + str({w}))]['display'].image\")".format(w=match.group(1))

    def format(self, match):
        return f"self.output({match.group(1)})"

    def output(self, result):
        self.parent.open_display_window(result, None, "Result")

    def close(self):
        self.parent.calculator_window = None
        self.destroy()
    
class Graph(tk.Toplevel):
    def __init__(self, parent: TopWindow=None):
        super().__init__()
        self.parent = parent
        self.title("Graph")
        self.geometry("400x300")
        self.attributes("-topmost", True)        
        self.data = None
        self.create_widget()
        
    def create_widget(self):
        options = ['Line', 'Horizontal Box', 'Area']
        self.control_frame1 = ttk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)

        self.control_frame2 = ttk.Frame(self)
        self.control_frame2.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.graph_type_label = ttk.Label(self.control_frame1, text="Graph Type:")
        self.graph_type_label.pack(side=tk.LEFT, fill=tk.X)

        self.graph_type = ttk.Combobox(self.control_frame1, values=options)
        self.figure, self.ax = plt.subplots()
        self.graph_type.pack(side=tk.LEFT, expand=True)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw_idle()
        
        self.save_button = ttk.Button(self.control_frame2, text="Save", command=self.save)
        self.save_button.pack(side=tk.LEFT)

    def update_graph(self, data):
        if self.graph_type.get() == "Line" or self.graph_type.get() == "Horizontal Box":
            self.figure.clear()
            self.data = data
            self.ax = self.figure.add_subplot(111)
            self.ax.plot(data)
            self.ax.set_title("Graph")
            self.ax.set_xlabel("X-axis")
            self.ax.set_ylabel("Y-axis")
            self.canvas.draw_idle()
        elif self.graph_type.get() == "Area":
            self.figure.clear()
            self.data = np.array(data)
            self.ax = self.figure.add_subplot(111)
            if data.shape[0] > 1 and data.shape[1] > 1:
                t = self.ax.imshow(data, cmap='viridis', interpolation='nearest')
                self.figure.colorbar(t)
            self.ax.set_title("Graph")
            self.ax.set_xlabel("X-axis")
            self.ax.set_ylabel("Y-axis")
            self.canvas.draw_idle()
        
    def close(self):
        self.parent.graph_window = None
        self.destroy()
    
    def save(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if file_path:
            df = pd.DataFrame(self.data)
            df.to_csv(file_path, index=False)
            messagebox.showinfo("Saved", f"Data saved to {file_path}")

class Histogram(tk.Toplevel):
    def __init__(self, parent: TopWindow=None):
        super().__init__()
        self.parent = parent
        self.title("Histogram")
        self.geometry("400x300")
        self.parent.histogram_window = self
        
        self.mean = None
        self.sigma = None
        self.vmin = tk.IntVar()
        self.vmax = tk.IntVar()
        self.selected_window = None
        self.scatter_plot = None
        
        self.create_widget()
        
    def create_widget(self):
        self.control_frame1 = ttk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = ttk.Frame(self)
        self.control_frame2.pack(side=tk.BOTTOM, fill=tk.X)    
        self.control_frame3 = ttk.Frame(self)
        self.control_frame3.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.display_window_select = ttk.Combobox(self.control_frame1, values=[window.title() for window in self.parent.display_windows])
        self.display_window_select.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.display_window_select.bind("<ButtonPress-1>", self.update_selections)
        self.display_window_select.bind("<<ComboboxSelected>>", self.update_selected_window)
        
        self.compute_button = ttk.Button(self.control_frame1, text="Compute", command=self.compute)
        self.compute_button.config(state=tk.DISABLED)
        self.compute_button.pack(side=tk.LEFT)
        
        
        self.select_range_slider1_label = ttk.Label(self.control_frame2, text="VMin:")
        self.select_range_slider1_label.pack(side=tk.LEFT, fill=tk.X)
        self.select_range_slider1_entry = ttk.Entry(self.control_frame3, textvariable=self.vmin)
        self.select_range_slider1_entry.bind("<Return>", self.update_norm)
        self.select_range_slider1_entry.pack(side=tk.LEFT, expand=True)
        self.select_range_slider1 = ttk.Scale(self.control_frame2, from_=0, to=100, orient=tk.HORIZONTAL, variable=self.vmin, command=self.update_range)
        self.select_range_slider1.pack(side=tk.LEFT, fill=tk.X, expand=True)
        # self.select_range_slider1.set(0)
        self.select_range_slider1.config(state=tk.DISABLED)
        self.select_range_slider1.bind("<ButtonRelease-1>", self.update_norm)
        
        self.select_range_slider2_label = ttk.Label(self.control_frame2, text="VMax:")
        self.select_range_slider2_label.pack(side=tk.LEFT, fill=tk.X)
        self.select_range_slider2_entry = ttk.Entry(self.control_frame3, textvariable=self.vmax)
        self.select_range_slider2_entry.bind("<Return>", self.update_norm)  
        self.select_range_slider2_entry.pack(side=tk.LEFT, expand=True)
        self.select_range_slider2 = ttk.Scale(self.control_frame2, from_=0, to=100, orient=tk.HORIZONTAL, variable=self.vmax, command=self.update_range)
        self.select_range_slider2.pack(side=tk.LEFT, fill=tk.X, expand=True)
        # self.select_range_slider2.set(0)
        self.select_range_slider2.config(state=tk.DISABLED)
        self.select_range_slider2.bind("<ButtonRelease-1>", self.update_norm)
        
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw_idle()
        self.attributes("-topmost", True)
        
    def update_selections(self, event):
        available_image = [window.title() for window in self.parent.display_windows]
        event.widget['values'] = available_image
        event.widget.set(event.widget.get())

    def update_selected_window(self, event):
        selected_window = event.widget.get()
        available_image = [window.title() for window in self.parent.display_windows]
        self.selected_window = self.parent.display_windows[available_image.index(selected_window)]
        self.compute_button.config(state=tk.NORMAL)
        
        self.ax.clear()
    
    def update_range(self, event=None):
        if self.vmin.get() > self.vmax.get():
            self.vmin.set(self.vmax.get())
        if self.scatter_plot is not None:
            self.scatter_plot[0].remove()
            self.scatter_plot[1].remove()
        self.scatter_plot = (self.ax.scatter(self.vmin.get(), 0, color='red', label='VMin', marker='x'), self.ax.scatter(self.vmax.get(), 0, color='green', label='VMax', marker='x'))
        self.canvas.draw_idle()
    
    def update_norm(self, event):
        if self.vmin.get() > self.vmax.get():
            self.vmin.set(self.vmax.get())
        self.selected_window.norm.vmin = self.vmin.get()
        self.selected_window.norm.vmax = self.vmax.get()
        self.selected_window.canvas.draw_idle()
    
    def compute(self):
        self.data = np.nan_to_num(self.selected_window.image.flatten(), -99999)
        self.sigma = np.std(self.data)
        self.mean = np.mean(self.data)
        self.min_value = self.mean - self.sigma
        self.vmin.set(self.min_value)
        self.max_value = self.mean + self.sigma
        self.vmax.set(self.max_value)
        self.select_range_slider1.config(state=tk.NORMAL)
        self.select_range_slider1.config(from_=self.min_value, to=self.max_value)
        self.select_range_slider1.set(self.mean)
        self.select_range_slider2.config(state=tk.NORMAL)
        self.select_range_slider2.config(from_=self.min_value, to=self.max_value)
        self.select_range_slider2.set(self.mean)
        self.plot_histogram(self.data)        
    
    def plot_histogram(self, data):
        self.ax.clear()
        self.ax.hist(data, bins=256, color='blue', alpha=0.7, range=(self.min_value, self.max_value))
        self.scatter_plot = (self.ax.scatter(self.vmin.get(), 0, color='red', label='VMin', marker='x'), self.ax.scatter(self.vmax.get(), 0, color='green', label='VMax', marker='x'))
        self.ax.set_title("Histogram")
        self.ax.set_xlabel("Pixel Value")
        self.ax.set_ylabel("Frequency")
        self.canvas.draw_idle()
    
    def close(self):
        self.parent.histogram_window = None
        self.destroy()
        
class Aperture(tk.Toplevel):
    def __init__(self, parent: DisplayWindow=None, a=20, b=20, angle=0, gap=10, background=10):
        super().__init__()
        self.parent = parent
        self.title(self.parent.title())
        self.parent.aperture_window = self
        self.geometry("400x300")
        self.attributes("-topmost", True)
        
        self.enable_circle_aperture = tk.BooleanVar()
        self.create_widget(a, b, angle, gap, background)
        
    def create_widget(self, a, b, angle, gap, background):
        self.control_frame1 = ttk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = ttk.Frame(self)
        self.control_frame2.pack(side=tk.TOP, fill=tk.X)
        self.control_frame3 = ttk.Frame(self)
        self.control_frame3.pack(side=tk.TOP, fill=tk.X)
        self.control_frame4 = ttk.Frame(self)
        self.control_frame4.pack(side=tk.TOP, fill=tk.X)
        self.control_frame5 = ttk.Frame(self)
        self.control_frame5.pack(side=tk.TOP, fill=tk.X)
        self.control_frame6 = ttk.Frame(self)
        self.control_frame6.pack(side=tk.TOP, fill=tk.X)
        
        self.aperture_major_var = tk.StringVar()
        self.aperture_major_label = ttk.Label(self.control_frame1, text="Aperture Major axis:")
        self.aperture_major_label.pack(side=tk.TOP, fill=tk.X)
        self.aperture_major = ttk.Entry(self.control_frame1, textvariable=self.aperture_major_var)
        self.aperture_major.pack(side=tk.TOP, fill=tk.X)
        self.aperture_major.insert(0, str(a))
        self.aperture_major.bind("<Up>", self.increase_axis)
        self.aperture_major.bind("<Down>", self.decrease_axis)

        self.aperture_minor_var = tk.StringVar()
        self.aperture_minor_label = ttk.Label(self.control_frame2, text="Aperture Minor axis:")        
        self.aperture_minor_label.pack(side=tk.TOP, fill=tk.X)
        self.aperture_minor = ttk.Entry(self.control_frame2, textvariable=self.aperture_minor_var)
        self.aperture_minor.pack(side=tk.TOP, fill=tk.X)
        self.aperture_minor.insert(0, str(b))
        self.aperture_minor.bind("<Up>", self.increase_axis)
        self.aperture_minor.bind("<Down>", self.decrease_axis)
        
        self.aperture_major_var.trace_add("write", lambda *args: self.aperture_major_onchange())
        self.aperture_minor_var.trace_add("write", lambda *args: self.aperture_minor_onchange())
        
        self.aperture_angle_label = ttk.Label(self.control_frame3, text="Aperture Angle:")
        self.aperture_angle_label.pack(side=tk.TOP, fill=tk.X)
        self.aperture_angle = ttk.Entry(self.control_frame3)
        self.aperture_angle.pack(side=tk.TOP, fill=tk.X)        
        self.aperture_angle.insert(0, str(angle))
        self.aperture_angle.bind("<Up>", self.increase)
        self.aperture_angle.bind("<Down>", self.decrease)
        
        self.gap_label = ttk.Label(self.control_frame4, text="Gap:")
        self.gap_label.pack(side=tk.TOP, fill=tk.X)
        self.gap = ttk.Entry(self.control_frame4)
        self.gap.pack(side=tk.TOP, fill=tk.X)
        self.gap.insert(0, str(gap))
        self.gap.bind("<Up>", self.increase)
        self.gap.bind("<Down>", self.decrease)
        
        self.background_label = ttk.Label(self.control_frame5, text="Background:")
        self.background_label.pack(side=tk.TOP, fill=tk.X)
        self.background = ttk.Entry(self.control_frame5)
        self.background.pack(side=tk.TOP, fill=tk.X)
        self.background.insert(0, str(background))
        self.background.bind("<Up>", self.increase)
        self.background.bind("<Down>", self.decrease)
        
        self.enable_circle_aperture_checkbox = ttk.Checkbutton(self.control_frame6, text="Enable Circle Aperture?", variable=self.enable_circle_aperture)
        self.enable_circle_aperture.set(True)
        self.enable_circle_aperture_checkbox.pack(side=tk.RIGHT)

    def aperture_major_onchange(self):
        if self.enable_circle_aperture.get():
            if self.aperture_major_var.get() == "":
                return
            self.aperture_minor_var.set(self.aperture_major_var.get())
            x, y = self.parent.aperture[0].get_center()
            self.parent.draw_aperture(x, y)
            
    def aperture_minor_onchange(self):
        if self.enable_circle_aperture.get():
            if self.aperture_minor_var.get() == "":
                return
            self.aperture_major_var.set(self.aperture_minor_var.get())
            x, y = self.parent.aperture[0].get_center()
            self.parent.draw_aperture(x, y)

    def increase_axis(self, event):
        if self.enable_circle_aperture.get():
            old = event.widget.get()
            if old == "":
                old = 0
            else:
                old = float(old)
            self.aperture_major.delete(0, tk.END)
            self.aperture_major.insert(0, str(old+1))
            self.aperture_minor.delete(0, tk.END)
            self.aperture_minor.insert(0, str(old+1))
            x, y = self.parent.aperture[0].get_center()
            self.parent.draw_aperture(x, y)
        else:
            old = event.widget.get()
            if old == "":
                old = 0
            else:
                old = float(old)
            event.widget.delete(0, tk.END)
            event.widget.insert(0, str(old+1))
            x, y = self.parent.aperture[0].get_center()
            self.parent.draw_aperture(x, y)
            
    def decrease_axis(self, event):
        if self.enable_circle_aperture.get():
            old = event.widget.get()
            if old == "":
                old = 0
            else:
                old = float(old)
            self.aperture_major.delete(0, tk.END)
            self.aperture_major.insert(0, str(old-1))
            self.aperture_minor.delete(0, tk.END)
            self.aperture_minor.insert(0, str(old-1))
            x, y = self.parent.aperture[0].get_center()
            self.parent.draw_aperture(x, y)
        else:
            old = event.widget.get()
            if old == "":
                old = 0
            else:
                old = float(old)
            event.widget.delete(0, tk.END)
            event.widget.insert(0, str(old-1))
            x, y = self.parent.aperture[0].get_center()
            self.parent.draw_aperture(x, y)

    def increase(self, event):
        old = event.widget.get()
        if old == "":
            old = 0
        else:
            old = float(old)
        event.widget.delete(0, tk.END)
        event.widget.insert(0, str(old+1))
        x, y = self.parent.aperture[0].get_center()
        self.parent.draw_aperture(x, y)
    
    def decrease(self, event):
        old = event.widget.get()
        if old == "":
            old = 0
        else:
            old = float(old)
        event.widget.delete(0, tk.END)
        event.widget.insert(0, str(old-1))
        x, y = self.parent.aperture[0].get_center()
        self.parent.draw_aperture(x, y)
        
    def close(self):
        self.parent.aperture_window = None
        self.parent.adding_aperture = False
        if len(self.parent.aperture) > 0:
            [i.remove() for i in self.parent.aperture]
            self.parent.aperture = ()
            self.parent.canvas.draw_idle()
            self.parent.add_aperture_button.state(['!pressed'])
        self.destroy()
        
    def get_max(self, aperture: patches.Ellipse, image):
        rr, cc = skimage.draw.ellipse(aperture.center[0], aperture.center[1], aperture.width/2, aperture.height/2, shape=(image.shape[1], image.shape[0]), rotation=np.deg2rad(aperture.angle))
        if rr.shape[0] == 0 or cc.shape[0] == 0:
            return (aperture.center[0], aperture.center[1])
        rr_min = np.min(rr)
        rr_max = np.max(rr)
        cc_min = np.min(cc)
        cc_max = np.max(cc)
        image = image[cc_min:cc_max, rr_min:rr_max]
        coords = np.unravel_index(np.argmax(image, axis=None), image.shape)
        return (coords[1] + rr_min, coords[0] + cc_min)
    
class ObjectsWindow(tk.Toplevel):
    def __init__(self, parent: TopWindow=None):
        super().__init__()
        self.parent = parent
        self.title(self.parent.title())
        self.geometry("600x300")
        self.attributes("-topmost", True)
        self.parent.objects_window = self
        self.current_selection = None
        
        self.create_widget()
        
    def create_widget(self):
        self.control_frame1 = ttk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        
        self.object_table = ttk.Treeview(self.control_frame1, columns=("#0","Index", "Name", "X", "Y", "Intensity", "Note"), show="headings")
        for i in self.object_table["columns"]:
            if i == "#0":
                self.object_table.heading(i, text="#")
                self.object_table.column(i, width=50)
                continue
            self.object_table.heading(i, text=i, anchor=tk.CENTER)
            self.object_table.column(i, width=50)
        self.object_table.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.object_table.bind("<ButtonRelease-1>", self.on_left_click)
        self.object_table.bind("<Button-3>", self.on_right_click)
        self.object_table.bind("<Double-1>", self.on_double_click)
        self.object_table.bind("<Delete>", self.delete)
        self.update_table()
        
        self.control_frame2 = ttk.Frame(self)
        self.control_frame2.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.save_button = ttk.Button(self.control_frame2, text="Save", command=self.save)
        self.save_button.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
    def update_table(self):
        self.object_table.delete(*self.object_table.get_children())
        for display_window in self.parent.display_windows:
            iid = self.object_table.insert("", "end", values=(display_window.title().split("/")[-1], "", "", "", "", "", display_window.title()), open=True)
            n = 1
            for aperture in display_window.added_apertures:
                x, y = aperture[0].center
                intensity = self.get_intensity(aperture[0], aperture[1], aperture[2], display_window.image)
                self.object_table.insert(iid, "end", values=("", n, "", x, y, intensity, ""))
                n+=1
        
    def add(self, aperture: list|tuple, window: DisplayWindow):
        for w in self.object_table.get_children():
            if self.object_table.item(w, "values")[6] == window.title():
                if len(self.object_table.get_children(w)) == 0:
                    n = 1
                else:
                    last_item = self.object_table.get_children(w)[-1]
                    n = int(self.object_table.item(last_item, "values")[1]) + 1
                intensity = self.get_intensity(aperture[0], aperture[1], aperture[2], window.image)
                self.object_table.insert(w, "end", values=("", n, "", aperture[0].center[0], aperture[0].center[1], intensity, ""))
                break
            
    def add_window(self, window: DisplayWindow):
        self.object_table.insert("", "end", values=(window.title().split("/")[-1], "", "", "", "", "", window.title()), open=True)
                
    def delete(self, event: tk.Event):
        if len(self.object_table.selection()) > 0:
            item = self.object_table.selection()[0]
            if self.object_table.item(item, "values")[1] == "":
                for window in self.parent.display_windows:
                    if window.title() == self.object_table.item(item, "values")[6]:
                        for i in window.added_apertures:
                            for j in i:
                                j.remove()
                        window.added_apertures = []
                        window.canvas.draw_idle()
                        window.close()
                        self.object_table.delete(item)
                        return
            parent = self.object_table.parent(item)
            for window in self.parent.display_windows:
                if window.title() == self.object_table.item(parent, "values")[6]:
                    index = int(self.object_table.item(item, "values")[1])-1
                    for i in window.added_apertures[index]:
                        i.remove()
                    window.added_apertures.pop(index)
                    break
            self.object_table.delete(item)
            if self.current_selection is not None:
                self.current_selection[0].remove()
                self.current_selection = None
            window.canvas.draw_idle()
            for i, c in enumerate(self.object_table.get_children(parent)):
                self.object_table.set(c, column=1, value=i+1)
                      
    def on_double_click(self, event):
        col_n = self.object_table.identify_column(event.x)
        if col_n == "#3":
            selected_item = self.object_table.selection()[0]
            name = self.object_table.item(selected_item, "values")[2]
            name = tk.simpledialog.askstring("Edit Name", "Enter new name:", initialvalue=name)
            if name:
                self.object_table.set(selected_item, column=2, value=name)
        elif col_n == "#7":
            selected_item = self.object_table.selection()[0]
            note = self.object_table.item(selected_item, "values")[6]
            note = tk.simpledialog.askstring("Edit Note", "Enter new note:", initialvalue=note)
            if note:
                self.object_table.set(selected_item, column=6, value=note)
                      
    def on_left_click(self, event):
        if len(self.object_table.selection()) == 0:
            return
        selected_item = self.object_table.selection()[0]
        if self.object_table.item(selected_item, "values")[1] == "":
            return
        parent_title = self.object_table.item(self.object_table.parent(selected_item), "values")[6]
        display_window = None
        if self.current_selection is not None:
            self.current_selection[0].remove()
        for window in self.parent.display_windows:
            window.canvas.draw_idle()
            if window.title() == parent_title:
                display_window = window
        if display_window is None:
            return
        _, index, name, x, y, intensity, note = self.object_table.item(selected_item, "values")
        x = int(x)
        y = int(y)
        intensity = float(intensity)
        self.current_selection = display_window.ax.plot(x, y, 'ro')
        display_window.deiconify()
        display_window.lift()
        self.object_table.focus_force()
        display_window.ax.set_aspect('equal')
        display_window.canvas.draw_idle()
        
    def on_right_click(self, event):
        if self.current_selection is not None:
            self.current_selection[0].remove()
            for window in self.parent.display_windows:
                window.canvas.draw_idle()
            self.current_selection = None
        
    def get_intensity(self, inner_aperture, gap_aperture, outer_aperture, image):
        inner_rr, inner_cc = skimage.draw.ellipse(inner_aperture.center[0], inner_aperture.center[1], inner_aperture.width/2, inner_aperture.height/2, shape=(image.shape[1], image.shape[0]), rotation=np.deg2rad(inner_aperture.angle))
        gap_rr, gap_cc = skimage.draw.ellipse(gap_aperture.center[0], gap_aperture.center[1], gap_aperture.width/2, gap_aperture.height/2, shape=(image.shape[1], image.shape[0]), rotation=np.deg2rad(gap_aperture.angle))
        outer_rr, outer_cc = skimage.draw.ellipse(outer_aperture.center[0], outer_aperture.center[1], outer_aperture.width/2, outer_aperture.height/2, shape=(image.shape[1], image.shape[0]), rotation=np.deg2rad(outer_aperture.angle))
        inner_mask = np.zeros(image.shape, dtype=bool)
        inner_mask[inner_cc, inner_rr] = True
        gap_mask = np.zeros(image.shape, dtype=bool)
        gap_mask[gap_cc, gap_rr] = True
        outer_mask = np.zeros(image.shape, dtype=bool)
        outer_mask[outer_cc, outer_rr] = True
        background_mask = outer_mask & ~inner_mask & ~gap_mask
        background = np.mean(image[background_mask])
        intensity = np.sum(image[inner_mask]-background)
        return intensity
    
    def close(self):
        self.parent.objects_window = None
        self.destroy()
    
    def save(self):
        df = pd.DataFrame(columns=["#0", "Name", "X", "Y", "Intensity", "Note"])
        current = ""
        for group in self.object_table.get_children():
            df.loc[len(df)] = [self.object_table.item(group, "values")[0], "", "", "", "", ""]
            for item in self.object_table.get_children(group):
                _, index, name, x, y, intensity, note = self.object_table.item(item, "values")
                x = int(x)
                y = int(y)
                intensity = float(intensity)
                df.loc[len(df)] = ["", name, x, y, intensity, note]
        file_path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if file_path:
            df.to_csv(file_path, index=False)
            messagebox.showinfo("Saved", f"Data saved to {file_path}")
    
if __name__ == "__main__":
    # mpl.rcParams['path.simplify'] = True
    # mpl.rcParams['path.simplify_threshold'] = 1.0
    app = TopWindow()
    app.mainloop()