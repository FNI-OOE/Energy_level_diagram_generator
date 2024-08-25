import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import json
from rdkit import Chem
from rdkit.Chem import Draw
import os
from PIL import Image, ImageTk
import matplotlib.ticker as ticker

class EnergyDiagramApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Energy Diagram Generator")
        self.root.state('zoomed')  # Start the program in a maximized window
        self.molecules = []

        # Top Frame for image display
        self.top_frame = tk.Frame(root)
        self.top_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas = tk.Canvas(self.top_frame, bg='white')
        self.canvas.pack(fill=tk.BOTH, expand=True)

        # Bottom Frame for parameters and controls
        self.bottom_frame = tk.Frame(root)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

        self.param_frame = tk.LabelFrame(self.bottom_frame, text="Molecules & Parameters")
        self.param_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Control buttons on the bottom right
        self.control_frame = tk.Frame(self.bottom_frame)
        self.control_frame.pack(side=tk.RIGHT, fill=tk.Y)

        button_width = 10
        button_height = 2
        pad_x = 5
        pad_y = 5

        self.generate_button = tk.Button(self.control_frame, text="Generate", bg="lightgrey", fg="black", font=("Arial", 10, "bold"), command=self.generate_image,
                                         width=button_width, height=button_height)
        self.generate_button.pack(side=tk.TOP, padx=pad_x, pady=pad_y)

        self.save_image_button = tk.Button(self.control_frame, text="Save Image", bg="lightgrey", fg="black", font=("Arial", 10, "bold"), command=self.save_image,
                                           width=button_width, height=button_height)
        self.save_image_button.pack(side=tk.TOP, padx=pad_x, pady=pad_y)

        self.save_config_button = tk.Button(self.control_frame, text="Save Config", bg="lightgrey", fg="black", font=("Arial", 10, "bold"), command=self.save_config,
                                            width=button_width, height=button_height)
        self.save_config_button.pack(side=tk.TOP, padx=pad_x, pady=pad_y)

        self.load_config_button = tk.Button(self.control_frame, text="Load Config", bg="lightgrey", fg="black", font=("Arial", 10, "bold"), command=self.load_config,
                                            width=button_width, height=button_height)
        self.load_config_button.pack(side=tk.TOP, padx=pad_x, pady=pad_y)


        # Add the initial molecule frame
        self.add_molecule()

    def add_molecule(self):
        mol_frame = tk.Frame(self.param_frame, bd=2, relief=tk.RIDGE)
        mol_frame.grid(row=0, column=len(self.molecules) + 1, padx=5, pady=5)

        mol = {}

        mol['frame'] = mol_frame
        mol['name'] = tk.StringVar()
        mol['S1'] = tk.DoubleVar()
        mol['T1'] = tk.DoubleVar()
        mol['S2'] = tk.DoubleVar()
        mol['T2'] = tk.DoubleVar()
        mol['S1_color'] = tk.StringVar(value="#000000")
        mol['T1_color'] = tk.StringVar(value="#000000")
        mol['S2_color'] = tk.StringVar(value="#000000")
        mol['T2_color'] = tk.StringVar(value="#000000")
        mol['mol_file'] = None

        self.molecules.append(mol)
    
        #Width/height of input boxes, color boxes, add and delete molecule buttons
        self.button_h = 1
        self.button_w = 2
        self.box_h = 1
        self.box_w = 2
        self.input_energy = 5
        self.input_molecule = 9

        # Molecule Name and File Loader
        tk.Label(mol_frame, text="        Molecule:").grid(row=0, column=0)
        entry = tk.Entry(mol_frame, width=self.input_molecule, textvariable=mol['name'])
        entry.grid(row=0, column=1, columnspan=3, sticky='w')
        tk.Button(mol_frame, text="Load", bg="lightgrey", fg="black", font=("Arial", 10, "bold"), 
                  padx=1, width=4, height=1,
                  command=lambda: self.load_mol_file(mol)).grid(row=0, column=4)

        # Add molecule button
        add_button = tk.Button(mol_frame, text="+", bg="green", fg="white", font=("Arial", 10, "bold"),
                               width=self.button_w, height=self.button_h, command=self.add_molecule)
        add_button.grid(row=1, column=4, sticky="ne", padx=5, pady=5)

        # S1 Energy and Color
        tk.Label(mol_frame, text="S1 Energy (eV):").grid(row=1, column=0)
        tk.Entry(mol_frame, width=self.input_energy, textvariable=mol['S1']).grid(row=1, column=1)
        mol['S1_color_box'] = tk.Label(mol_frame, bg=mol['S1_color'].get(), width=self.box_w, height=self.box_h, relief=tk.RIDGE)
        mol['S1_color_box'].grid(row=1, column=2, padx=5, pady=5)
        mol['S1_color_box'].bind("<Button-1>", lambda e: self.choose_color(mol, 'S1'))

        # T1 Energy and Color
        tk.Label(mol_frame, text="T1 Energy (eV):").grid(row=2, column=0)
        tk.Entry(mol_frame, width=self.input_energy, textvariable=mol['T1']).grid(row=2, column=1)
        mol['T1_color_box'] = tk.Label(mol_frame, bg=mol['T1_color'].get(), width=self.box_w, height=self.box_h, relief=tk.RIDGE)
        mol['T1_color_box'].grid(row=2, column=2, padx=5, pady=5)
        mol['T1_color_box'].bind("<Button-1>", lambda e: self.choose_color(mol, 'T1'))

        # S2 Energy and Color
        tk.Label(mol_frame, text="S2 Energy (eV):").grid(row=3, column=0)
        tk.Entry(mol_frame, width=self.input_energy, textvariable=mol['S2']).grid(row=3, column=1)
        mol['S2_color_box'] = tk.Label(mol_frame, bg=mol['S2_color'].get(), width=self.box_w, height=self.box_h, relief=tk.RIDGE)
        mol['S2_color_box'].grid(row=3, column=2, padx=5, pady=5)
        mol['S2_color_box'].bind("<Button-1>", lambda e: self.choose_color(mol, 'S2'))

        # T2 Energy and Color
        tk.Label(mol_frame, text="T2 Energy (eV):").grid(row=4, column=0)
        tk.Entry(mol_frame, width=self.input_energy, textvariable=mol['T2']).grid(row=4, column=1)
        mol['T2_color_box'] = tk.Label(mol_frame, bg=mol['T2_color'].get(), width=self.box_w, height=self.box_h, relief=tk.RIDGE)
        mol['T2_color_box'].grid(row=4, column=2, padx=5, pady=5)
        mol['T2_color_box'].bind("<Button-1>", lambda e: self.choose_color(mol, 'T2'))

        # Red X button to remove the molecule
        remove_button = tk.Button(mol_frame, text="X", bg="red", fg="white", font=("Arial", 10, "bold"),
                                  width=self.button_w, height=self.button_h, command=lambda: self.remove_molecule(mol))
        remove_button.grid(row=4, column=4, sticky="se", padx=5, pady=5)

    def remove_molecule(self, mol):
        if len(self.molecules) == 1:
            return
        mol['frame'].destroy()
        self.molecules.remove(mol)
        self.rearrange_molecules()

    def rearrange_molecules(self):
        for i, mol in enumerate(self.molecules):
            mol['frame'].grid(row=0, column=i + 1, padx=5, pady=5)

    def load_mol_file(self, mol):
        mol_file_path = filedialog.askopenfilename(filetypes=[("Molecule Files", "*.mol")])
        if mol_file_path:
            mol['mol_file'] = Chem.MolFromMolFile(mol_file_path)
            mol_name = os.path.splitext(os.path.basename(mol_file_path))[0]  # Get filename minus extension
            mol['name'].set(mol_name)

    def choose_color(self, mol, level):
        color_code = colorchooser.askcolor(title="Choose color")[1]
        if level == 'S1':
            mol['S1_color'].set(color_code)
            mol['S1_color_box'].config(bg=color_code)
        elif level == 'T1':
            mol['T1_color'].set(color_code)
            mol['T1_color_box'].config(bg=color_code)
        elif level == 'S2':
            mol['S2_color'].set(color_code)
            mol['S2_color_box'].config(bg=color_code)
        elif level == 'T2':
            mol['T2_color'].set(color_code)
            mol['T2_color_box'].config(bg=color_code)

    def generate_image(self):
        fig = self.plot_diagram()  # Get the figure object
        self.display_image(fig)  # Display it in the top section

    def plot_diagram(self, filename=None):
        fig, ax = plt.subplots(1, len(self.molecules), sharey=True, figsize=(3 * len(self.molecules), 5))
        plt.subplots_adjust(wspace=0.5)

        if len(self.molecules) == 1:
            ax = [ax]  # Ensure ax is iterable even with one subplot

        # Initialize a list to keep track of maximum y-values
        all_y_values = []

        for i, mol in enumerate(self.molecules):
            S0 = 0
            S1 = mol['S1'].get()
            T1 = mol['T1'].get()
            S2 = mol['S2'].get()
            T2 = mol['T2'].get()
            S1_color = mol['S1_color'].get()
            T1_color = mol['T1_color'].get()
            S2_color = mol['S2_color'].get()
            T2_color = mol['T2_color'].get()
            name = mol['name'].get()

            # Add energy levels to the list for later max calculation
            all_y_values.extend([S1, T1, S2, T2])

            ax[i].set_xlim(0, 1.25)  # Set x-axis limits

            # Plot S0 level (always visible)
            ax[i].hlines(S0, 0, 1, colors='black', linewidth=2)
            ax[i].text(-0.08, S0, 'S0', ha='center', va='center', fontsize=10)

            # Plot S1 level if not 0
            if S1 != 0:
                ax[i].hlines(S1, 0, 0.45, colors=S1_color, linewidth=2)
                ax[i].text(0.225, S1+0.05, f"{S1:.2f} eV", ha='center', va='bottom')
                ax[i].text(-0.08, S1, 'S1', ha='center', va='center', fontsize=10)

            # Plot T1 level if not 0
            if T1 != 0:
                ax[i].hlines(T1, 0.55, 1, colors=T1_color, linewidth=2)
                ax[i].text(0.775, T1+0.05, f"{T1:.2f} eV", ha='center', va='bottom')
                ax[i].text(1.08, T1, 'T1', ha='center', va='center', fontsize=10)

            # Plot S2 level if not 0
            if S2 != 0:
                ax[i].hlines(S2, 0, 0.45, colors=S2_color, linewidth=2)
                ax[i].text(0.225, S2+0.05, f"{S2:.2f} eV", ha='center', va='bottom')
                ax[i].text(-0.08, S2, 'S2', ha='center', va='center', fontsize=10)

            # Plot T2 level if not 0
            if T2 != 0:
                ax[i].hlines(T2, 0.55, 1, colors=T2_color, linewidth=2)
                ax[i].text(0.775, T2+0.05, f"{T2:.2f} eV", ha='center', va='bottom')
                ax[i].text(1.08, T2, 'T2', ha='center', va='center', fontsize=10)

            # Add a vertical dashed line between molecules
            if i < len(self.molecules) - 1:
                ax[i].axvline(x=1.2, color='grey', linestyle='--', linewidth=2)

            if i != 0:
                ax[i].axis('off')  # Remove plot frame and ticks except for first plot

            # Add molecule name below the subplot
            ax[i].text(0.5, S0-0.3, name, ha='center', va='top', fontsize=16, fontname='Arial')

            # Add molecule structure if available
            if mol['mol_file']:
                img = Draw.MolToImage(mol['mol_file'], size=(100, 100))
                imagebox = OffsetImage(img, zoom=0.5)
                ab = AnnotationBbox(imagebox, (0.5, S0-1), frameon=False, box_alignment=(0.5, 0))
                ax[i].add_artist(ab)

        # Determine the global y-axis limit
        y_max = max(all_y_values)

        # Add scale to the first subplot (ax[0])
        ax[0].yaxis.set_major_locator(ticker.MultipleLocator(1))
        ax[0].yaxis.set_minor_locator(ticker.MultipleLocator(0.5))

        # Customize the left spine
        ax[0].spines['left'].set_color('grey')
        ax[0].spines['left'].set_linewidth(2)
        ax[0].spines['left'].set_position(('outward', 30))  # Move the left spine outward
        ax[0].spines['left'].set_visible(True)  # Make sure the left spine is visible
        ax[0].spines['left'].set_bounds(0, y_max + 0.5)  # Set the bounds for the left spine

        # Hide the top, bottom, and right spines
        ax[0].spines['top'].set_visible(False)
        ax[0].spines['bottom'].set_visible(False)
        ax[0].spines['right'].set_visible(False)

        # Customize tick parameters
        ax[0].tick_params(axis='both', which='major', color='grey', direction='in', width=2, length=8)  # Major ticks
        ax[0].tick_params(axis='both', which='minor', color='grey', direction='in', width=2, length=4)   # Minor ticks

        # Hide the ticks for the top, bottom, and right
        ax[0].xaxis.set_ticks([])
        ax[0].yaxis.set_ticks_position('left')  # Keep ticks only on the left side

        # Add arrow at the top of the scale
        ax[0].annotate('', xy=(0, y_max + 0.5), xytext=(0, y_max + 0.5),
                       arrowprops=dict(arrowstyle="->", facecolor='grey', linewidth=2))

        # Add scale label "E (eV)" at the top of the scale
        ax[0].text(-0.18, y_max + 0.6, 'E (eV)', ha='center', va='bottom', fontsize=12)

        plt.tight_layout()

        if filename:
            plt.savefig(filename)

        return fig  # Return the figure object

    def display_image(self, fig):
        # Clear the canvas before drawing the new image
        self.canvas.delete("all")

        # Draw the figure on the canvas
        canvas = FigureCanvasAgg(fig)
        canvas.draw()

        # Convert canvas to image
        buf = canvas.buffer_rgba()
        width, height = fig.get_size_inches() * fig.get_dpi()
        image = Image.frombuffer('RGBA', (int(width), int(height)), buf, 'raw', 'RGBA', 0, 1)
        image_tk = ImageTk.PhotoImage(image)

        # Display the image on the canvas
        self.canvas.create_image(0, 0, image=image_tk, anchor=tk.NW)
        self.canvas.image = image_tk  # Keep a reference to avoid garbage collection

    def save_image(self):
        filename = filedialog.asksaveasfilename(defaultextension=".png")
        if filename:
            fig = self.plot_diagram(filename=filename)

    def save_config(self):
        filename = filedialog.asksaveasfilename(defaultextension=".json")
        if filename:
            config = []
            for mol in self.molecules:
                mol_config = {
                    'name': mol['name'].get(),
                    'S1': mol['S1'].get(),
                    'T1': mol['T1'].get(),
                    'S2': mol['S2'].get(),
                    'T2': mol['T2'].get(),
                    'S1_color': mol['S1_color'].get(),
                    'T1_color': mol['T1_color'].get(),
                    'S2_color': mol['S2_color'].get(),
                    'T2_color': mol['T2_color'].get(),
                    'mol_file': mol['mol_file']  # Store the molecule file path, not the RDKit object
                }
                config.append(mol_config)
            with open(filename, 'w') as file:
                json.dump(config, file)

    def load_config(self):
        filename = filedialog.askopenfilename(defaultextension=".json")
        if filename:
            with open(filename, 'r') as file:
                config = json.load(file)
            
            # Clear existing molecules
            for mol in self.molecules:
                for widget in mol['frame'].winfo_children():
                    widget.destroy()
                mol['frame'].destroy()
            self.molecules.clear()
            
            # Load new molecules from config
            for mol in config:
                self.add_molecule()
                last_mol = self.molecules[-1]
                last_mol['name'].set(mol['name'])
                last_mol['S1'].set(mol['S1'])
                last_mol['T1'].set(mol['T1'])
                last_mol['S2'].set(mol['S2'])
                last_mol['T2'].set(mol['T2'])
                last_mol['S1_color'].set(mol['S1_color'])
                last_mol['T1_color'].set(mol['T1_color'])
                last_mol['S2_color'].set(mol['S2_color'])
                last_mol['T2_color'].set(mol['T2_color'])
                if mol['mol_file']:
                    last_mol['mol_file'] = Chem.MolFromMolFile(mol['mol_file'])
            
            # Refresh the display
            self.generate_image()

# Initialize Tkinter root and the EnergyDiagramApp
root = tk.Tk()
app = EnergyDiagramApp(root)
root.mainloop()
