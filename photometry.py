import tkinter as tk
from tkinter import ttk
from tkinter import filedialog, messagebox
from astropy.io import fits
from astropy.visualization import ImageNormalize, LinearStretch
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
import numpy as np
import regex as re

class TopWindow(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Photometry Tool")
        self.geometry("800x100")
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
        
        self.open_data_button = tk.Button(self, text="Open FITS File", command=self.open_fits_file)
        self.open_data_button.pack(side=tk.LEFT)
        
        self.information_button = tk.Button(self, text="Information", command=self.show_information)
        self.information_button.pack(side=tk.LEFT)
        
        self.calculator_button = tk.Button(self, text="Calculator", command=self.open_calculator)
        self.calculator_button.pack(side=tk.LEFT)
        
        self.graph_button = tk.Button(self, text="Graph", command=self.show_graph)
        self.graph_button.pack(side=tk.LEFT)
        
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

class DisplayWindow(tk.Toplevel):
    def __init__(self, image, header, title, parent=None):
        super().__init__()
        self.title(title)
        self.geometry("800x600")
        self.image = image
        self.header = header
        self.parent = parent
        self.header_window = None
        
        self.pan_start = None
        
        self.create_widget()

    def create_widget(self):
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.draw()
        self.display_image()
        
        self.control_frame1 = tk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.control_frame2 = tk.Frame(self)
        self.control_frame2.pack(side=tk.TOP, fill=tk.X)
        
        self.header_button = tk.Button(self.control_frame1, text="Show Header", command=self.show_header)
        self.header_button.pack(side=tk.LEFT)
        
        self.positionX_label = tk.Label(self.control_frame2, text="X:")
        self.positionX_label.pack(side=tk.LEFT)
        self.positionY_label = tk.Label(self.control_frame2, text="Y:")
        self.positionY_label.pack(side=tk.LEFT)
        
        self.canvas.mpl_connect("motion_notify_event", self.on_mouse_move)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.canvas.mpl_connect('button_press_event', self.start_pan)
        self.canvas.mpl_connect('button_release_event', self.stop_pan)
        self.canvas.mpl_connect('motion_notify_event', self.on_pan)
        
    def show_header(self):
        if self.header is not None:
            self.header_window = HeaderWindow(self.header, self)
            self.header_window.protocol("WM_DELETE_WINDOW", self.header_window.close)
        else:
            messagebox.showerror("Error", "No header information available.")
        # self.header_window.mainloop()

    def display_image(self):
        norm = ImageNormalize(self.image, stretch=LinearStretch())
        self.ax.imshow(self.image, cmap='gray', norm=norm)
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
            pixel_value = self.image[y, x]
            if self.parent.information_window:
                self.parent.information_window.pixel_value_label.config(text=f"Pixel Value: {pixel_value}")
            self.positionX_label.config(text=f"X: {x}")
            self.positionY_label.config(text=f"Y: {y}")
        
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

    def start_pan(self, event):
        if event.button == 1:
            self.pan_start = (event.x, event.y)
    
    def stop_pan(self, event):
        if event.button == 1:
            self.pan_start = None
    
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
        
    def close(self):
        if self.header_window:
            self.header_window.close()
        self.parent.display_windows.remove(self)
        self.destroy()

class InformationWindow(tk.Toplevel):
    def __init__(self, parent=None):
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
    def __init__(self, header, parent=None):
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
    def __init__(self, parent=None):
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
                print(i.cget("text"))
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
    def __init__(self, parent=None):
        super().__init__()
        self.parent = parent
        self.title("Graph")
        self.geometry("400x300")
        self.attributes("-topmost", True)        
        self.create_widget()
        
    def create_widget(self):
        options = ['Line']
        self.control_frame1 = tk.Frame(self)
        self.control_frame1.pack(side=tk.TOP, fill=tk.X)
        self.graph_type_label = tk.Label(self.control_frame1, text="Graph Type:")
        self.graph_type_label.pack(side=tk.LEFT, fill=tk.X)
        self.graph_type = ttk.Combobox(self.control_frame1, values=options)
        self.graph_type.pack(side=tk.LEFT, expand=True)
        
    def close(self):
        self.parent.graph_window = None
        self.destroy()
if __name__ == "__main__":
    app = TopWindow()
    app.mainloop()