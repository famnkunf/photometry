import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from astropy.io import fits
from astropy.visualization import ImageNormalize, LinearStretch
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import patches
# from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
import numpy as np
import regex as re
import time
import threading
import pandas as pd

class TopWindow(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Photometry Tool")
        self.geometry("800x50")
        self.attributes("-topmost", True)
        self.create_widget()
        self.display_windows = []
        self.protocol("WM_DELETE_WINDOW", self.close)
        self.information_window = None
        self.calculator_window = None
        self.graph_window = None
        
    def create_widget(self):
        self.control_frame = tk.Frame(self)
        self.control_frame.pack(side=tk.TOP, fill=tk.X)
        
        self.open_data_button = tk.Button(self.control_frame, text="Open FITS File", command=self.open_fits_file)
        self.open_data_button.pack(side=tk.LEFT)
        
        self.information_button = tk.Button(self.control_frame, text="Information", command=self.show_information)
        self.information_button.pack(side=tk.LEFT)
        
        self.calculator_button = tk.Button(self.control_frame, text="Calculator", command=self.open_calculator)
        self.calculator_button.pack(side=tk.LEFT)
        
        self.graph_button = tk.Button(self.control_frame, text="Graph", command=self.show_graph)
        self.graph_button.pack(side=tk.LEFT)
        
        self.histogram_button = tk.Button(self.control_frame, text="Histogram", command=self.show_histogram)
        self.histogram_button.pack(side=tk.LEFT)
        
        
    def open_calculator(self):
        if self.calculator_window:
            self.calculator_window.focus()
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
        self.display_windows.append(DisplayWindow(image, header, file_path, self))
        self.display_windows[-1].protocol("WM_DELETE_WINDOW", self.display_windows[-1].close)
        # self.display_windows[-1].mainloop()
        
    def close(self):
        if self.information_window:
            self.information_window.close()
        for window in self.display_windows:
            window.close()
        self.quit()
        self.destroy()
        
    def show_information(self):
        if self.information_window:
            self.information_window.focus()
        else:
            self.information_window = InformationWindow(self)
            self.information_window.protocol("WM_DELETE_WINDOW", self.information_window.close)

    def show_graph(self):
        if self.graph_window:
            self.graph_window.focus()
        else:
            self.graph_window = Graph(self)
            self.graph_window.protocol("WM_DELETE_WINDOW", self.graph_window.close)

    def show_histogram(self):
        if self.graph_window:
            self.graph_window.focus()
        else:
            self.histogram_window = Histogram(self)
            self.histogram_window.protocol("WM_DELETE_WINDOW", self.histogram_window.close)
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
                
        self.create_widget()

    def create_widget(self):
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw()
        self.display_image()
        self.drawing_type = False
        self.coords1 = None
        self.coords2 = None
        self.annotate = None
        self.data = None
        self.last_update_time = time.time()
        self.target_interval = 1/10
        
        self.control_frame1 = tk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = tk.Frame(self)
        self.control_frame2.pack(side=tk.TOP, fill=tk.X)
        
        self.header_button = tk.Button(self.control_frame1, text="Show Header", command=self.show_header)
        self.header_button.pack(side=tk.LEFT)
        
        self.save_button = tk.Button(self.control_frame1, text="Save", command=self.save)
        self.save_button.pack(side=tk.LEFT)
        
        self.positionX_label = tk.Label(self.control_frame2, text="X:")
        self.positionX_label.pack(side=tk.LEFT)
        self.positionY_label = tk.Label(self.control_frame2, text="Y:")
        self.positionY_label.pack(side=tk.LEFT)

        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('button_press_event', self.start_pan)
        self.canvas.mpl_connect('button_release_event', self.stop_pan)
        self.canvas.mpl_connect('motion_notify_event', self.on_pan)
        
        self.bind("<Button-1>", self.on_focus)
        
    def show_header(self):
        if self.header is not None:
            self.header_window = HeaderWindow(self.header, self)
            self.header_window.protocol("WM_DELETE_WINDOW", self.header_window.close)
        else:
            messagebox.showerror("Error", "No header information available.")
        # self.header_window.mainloop()

    def display_image(self):
        self.ax.imshow(self.image, cmap='gray', norm=self.norm)
        self.ax.set_title("FITS Image")
        self.ax.set_xlim(0, self.image.shape[1])
        self.ax.set_ylim(0, self.image.shape[0])
        self.ax.set_aspect('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.canvas.draw()
        
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
        self.canvas.draw()

    def on_focus(self, event):
        if self.parent.graph_window:
            self.drawing_type = self.parent.graph_window.graph_type.get()
        else:
            self.drawing_type = False

    def start_pan(self, event):
        if event.button == 1:
            if not self.drawing_type:
                self.pan_start = (event.x, event.y)
            
    def stop_pan(self, event):
        if event.button == 1:
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
            self.canvas.draw()
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
            self.canvas.draw()
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
        self.canvas.draw()
        
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
        
    def save(self):
        file_path = filedialog.asksaveasfilename(defaultextension=".fit", filetypes=[("FITS Files", "*.fits, *.fit")])
        if file_path:
            fits.writeto(file_path, self.image, header=self.header, overwrite=True)
            messagebox.showinfo("Saved", f"Image saved to {file_path}")
        
    def close(self):
        if self.header_window:
            self.header_window.close()
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
        self.add_new_variable_button = tk.Button(self, text="Add Variable", command=self.add_new_variable)
        self.add_new_variable_button.pack(side=tk.TOP, fill=tk.X)
        self.formular_label = tk.Label(self, text="Formula:")
        self.formular_label.pack(side=tk.TOP, fill=tk.X)
        self.formular_entry = tk.Entry(self)
        self.formular_entry.pack(side=tk.TOP, fill=tk.X)
        self.calculate_button = tk.Button(self, text="Calculate", command=self.calculate)
        self.calculate_button.pack(side=tk.TOP, fill=tk.X)
    
    def add_new_variable(self):
        available_image = [window.title() for window in self.parent.display_windows]
        control_frame = tk.Frame(self)
        control_frame.pack(side=tk.TOP, fill=tk.X)
        
        select_variable_label = tk.Label(control_frame, text=f"Variable ${self.n_variables}")
        select_variable_label.pack(side=tk.LEFT)

        select_variable = ttk.Combobox(control_frame, values=available_image)
        select_variable.pack(side=tk.LEFT, fill=tk.X, expand=True)
        select_variable.bind("<ButtonPress-1>", self.update_selections)
        select_variable.bind("<<ComboboxSelected>>", self.update_variable_name)
        
        self.variables[f'${self.n_variables}'] = {}
        self.variables[f'${self.n_variables}']['entry'] = select_variable
        self.n_variables += 1
        
    def update_selections(self, event):
        available_image = [window.title() for window in self.parent.display_windows]
        event.widget['values'] = available_image
        
    def update_variable_name(self, event):
        available_image = [window.title() for window in self.parent.display_windows]
        for i in event.widget.master.winfo_children():
            if isinstance(i, tk.Label):
                r = re.compile(r'Variable \$(\d+)')
                match = r.search(i.cget("text"))
        if match:
            self.variables[f'${match.group(1)}']['display'] = self.parent.display_windows[available_image.index(event.widget.get())]
    def calculate(self):
        formular = self.formular_entry.get()
        replace_dict = {
            "#sqrt(": "np.sqrt(",
            "#log(": "np.log(",
            "#exp(": "np.exp(",
            "#sin(": "np.sin(",
            "#cos(": "np.cos(",
            "#tan(": "np.tan(",
        }
        for key in self.variables.keys():
            if key in formular:
                formular = formular.replace(key, f"self.variables['{key}']['display'].image")
        for key in replace_dict.keys():
            if key in formular:
                formular = formular.replace(key, replace_dict[key])
        try:
            result = eval(formular)
        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {e, formular}")
            return
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
        options = ['Line', 'Horizontal Box']
        self.control_frame1 = tk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)

        self.control_frame2 = tk.Frame(self)
        self.control_frame2.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.graph_type_label = tk.Label(self.control_frame1, text="Graph Type:")
        self.graph_type_label.pack(side=tk.LEFT, fill=tk.X)


        self.graph_type = ttk.Combobox(self.control_frame1, values=options)
        self.figure, self.ax = plt.subplots()
        self.graph_type.pack(side=tk.LEFT, expand=True)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw()
        
        self.save_button = tk.Button(self.control_frame2, text="Save", command=self.save)
        self.save_button.pack(side=tk.LEFT)

    def update_graph(self, data):
        self.figure.clear()
        self.data = data
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(data)
        self.ax.set_title("Graph")
        self.ax.set_xlabel("X-axis")
        self.ax.set_ylabel("Y-axis")
        self.canvas.draw()
        
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
        self.control_frame1 = tk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = tk.Frame(self)
        self.control_frame2.pack(side=tk.BOTTOM, fill=tk.X)    
        self.control_frame3 = tk.Frame(self)
        self.control_frame3.pack(side=tk.BOTTOM, fill=tk.X)
        
        self.display_window_select = ttk.Combobox(self.control_frame1, values=[window.title() for window in self.parent.display_windows])
        self.display_window_select.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.display_window_select.bind("<ButtonPress-1>", self.update_selections)
        self.display_window_select.bind("<<ComboboxSelected>>", self.update_selected_window)
        
        self.compute_button = tk.Button(self.control_frame1, text="Compute", command=self.compute)
        self.compute_button.config(state=tk.DISABLED)
        self.compute_button.pack(side=tk.LEFT)
        
        self.select_range_slider1 = tk.Scale(self.control_frame2, from_=0, to=100, orient=tk.HORIZONTAL, label="Vmin", resolution=0.1, variable=self.vmin, repeatinterval=1000, command=self.update_range)
        self.select_range_slider1.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.select_range_slider1.set(0)
        self.select_range_slider1.config(state=tk.DISABLED)
        self.select_range_slider1.bind("<ButtonRelease-1>", self.update_norm)
        
        self.select_range_slider2 = tk.Scale(self.control_frame2, from_=0, to=100, orient=tk.HORIZONTAL, label="Vmax", resolution=0.1, variable=self.vmax, repeatinterval=1000, command=self.update_range)
        self.select_range_slider2.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.select_range_slider2.set(0)
        self.select_range_slider2.config(state=tk.DISABLED)
        self.select_range_slider2.bind("<ButtonRelease-1>", self.update_norm)
        
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw()
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
    
    def update_range(self, event):
        if self.vmin.get() > self.vmax.get():
            self.vmin.set(self.vmax.get())
        self.scatter_plot[0].remove()
        self.scatter_plot[1].remove()
        self.scatter_plot = (self.ax.scatter(self.vmin.get(), 0, color='red', label='VMin', marker='x'), self.ax.scatter(self.vmax.get(), 0, color='green', label='VMax', marker='x'))
        self.canvas.draw()
    
    def update_norm(self, event):
        if self.vmin.get() > self.vmax.get():
            self.vmin.set(self.vmax.get())
        self.selected_window.norm.vmin = self.vmin.get()
        self.selected_window.norm.vmax = self.vmax.get()
        self.selected_window.canvas.draw()
    
    def compute(self):
        self.data = self.selected_window.image.flatten()
        self.sigma = np.std(self.selected_window.image.flatten())
        self.mean = np.mean(self.selected_window.image.flatten())
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
        self.canvas.draw()
    
    def close(self):
        self.parent.histogram_window = None
        self.destroy()
    
if __name__ == "__main__":
    # mpl.rcParams['path.simplify'] = True
    # mpl.rcParams['path.simplify_threshold'] = 1.0
    app = TopWindow()
    app.mainloop()