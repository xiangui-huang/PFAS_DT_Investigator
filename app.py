"""
This program is designed to investigate the distribution and transformation of PFAS within systems consisting of single or multiple compartments.

@author: Xiangui Huang, ZIWR, Ben-Gurion University.

OCT. 07, 2025, Midreshet Ben_Gurion, Israel.
"""

import toga
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, CENTER
import pandas as pd
from pathlib import Path
import json

from . import FunctionM_1
from . import FunctionM_2
from importlib import resources #used in native app.





CONFIG_FILE = Path.home() / ".my_toga_app_config.json"


class MyApp(toga.App):
    def startup(self):
        self.main_window = toga.MainWindow(title="PFAS_DT_Investigator")

        
        
        # Define a manu command for About author.
        about_cmd = toga.Command(
            self.show_about,            # function to run
            text="About the author",           # text shown in menu
            tooltip="About the author",
            group=toga.Group.HELP,
            icon=None,                  # optional icon path
        )

        # Add command to the menu bar
        self.commands.add(about_cmd)

        #self.main_window.show()



        #-------------------------
        #step1: Input File box
        #-------------------------
        self.input_box = toga.Box(style=Pack(direction=COLUMN, margin=10))
        
        # Load saved configuration
        self.config = self.load_config()
               
        
        # --- Widgets ---
        #set up the input port for the inputdata file for DT analysis.
        self.path_input = toga.TextInput(value=self.config.get("inputdata_file_path"), placeholder='Enter inputdata CSV file path for DT analysis...', style=Pack(flex=1))
        self.load_button = toga.Button('Loading input file (CSV) for DT analysis', id='inputdata', on_press=self.load_csv, style=Pack(margin_left=5))
        self.flag_label = toga.Label('Try to check the received input file path...', style=Pack(margin_top=5))
        self.rows_label = toga.Label('Rows: -', style=Pack(margin_top=5))
        self.cols_label = toga.Label('Columns: -', style=Pack(margin_top=5))
        # --- Layouts input data file---
        path_box = toga.Box(children=[self.path_input, self.load_button], style=Pack(direction=ROW, margin=5))
        info_box = toga.Box(children=[self.flag_label, self.rows_label, self.cols_label], style=Pack(direction=ROW, margin=5, align_items=CENTER))



        #set up the input port for the inputdata file for simulation and calibration.
        self.path_input_1 = toga.TextInput(value=self.config.get("inputdata_file_path_1"), placeholder='Enter inputdata CSV file path for calibration...', style=Pack(flex=1))
        self.load_button_1 = toga.Button('Loading input file (CSV) for calibration', id='inputdata_1', on_press=self.load_csv, style=Pack(margin_left=5))
        self.flag_label_1 = toga.Label('Try to check the received input file path...', style=Pack(margin_top=5))
        self.rows_label_1 = toga.Label('Rows: -', style=Pack(margin_top=5))
        self.cols_label_1 = toga.Label('Columns: -', style=Pack(margin_top=5))
        # --- Layouts input data file---
        path_box_1 = toga.Box(children=[self.path_input_1, self.load_button_1], style=Pack(direction=ROW, margin=5))
        info_box_1 = toga.Box(children=[self.flag_label_1, self.rows_label_1, self.cols_label_1], style=Pack(direction=ROW, margin=5, align_items=CENTER))



        
        #set up the input port for the input PFAS transformation lib file 
        self.libfile_path_input = toga.TextInput(placeholder='Enter PFAS_trans_lib CSV file path... (default file is Pre_Trans_Lib.csv)', style=Pack(flex=1))
        self.libfile_load_button = toga.Button('Loading precursor trans lib file (CSV)', id = "libfile", on_press=self.load_csv, style=Pack(margin_left=5))
        self.flag_label_2 = toga.Label('Try to check the received input file path...', style=Pack(margin_top=5))
        self.rows_label_2 = toga.Label('Rows: -', style=Pack(margin_top=5))
        self.cols_label_2 = toga.Label('Columns: -', style=Pack(margin_top=5))       
        # --- Layouts input lib file---
        libfile_path_box = toga.Box(children=[self.libfile_path_input, self.libfile_load_button], style=Pack(direction=ROW, margin=5))
        libfile_info_box = toga.Box(children=[self.flag_label_2, self.rows_label_2, self.cols_label_2], style=Pack(direction=ROW, margin=5, align_items=CENTER))



        #set up the output folder.
        self.output_path_input = toga.TextInput(value=self.config.get("output_folder_path"), placeholder='Enter the output folder path...', style=Pack(flex=1))
        self.output_check_button = toga.Button('Checking output folder dir', on_press=self.check_folder, style=Pack(margin_left=5))
        self.foldercheck_label = toga.Label('Try to check the received folder path...', style=Pack(margin_top=5))
        # --- Layouts output folder---
        output_path_box = toga.Box(children=[self.output_path_input, self.output_check_button], style=Pack(direction=ROW, margin=5))
        output_info_box = toga.Box(children=[self.foldercheck_label], style=Pack(direction=ROW, margin=5))
        
        
        #a button, allowing the user to save received input and out file/folder path as default.
        save_file_button=toga.Button('Save the inputdata file path and the output folder path', on_press=self.save_file_folder, style=Pack(flex=1, color="green", font_weight="bold"))
        
             
        self.input_box.add(path_box)
        self.input_box.add(info_box)
        
        self.input_box.add(path_box_1)
        self.input_box.add(info_box_1)
        
        self.input_box.add(libfile_path_box)
        self.input_box.add(libfile_info_box)
        
        self.input_box.add(output_path_box)
        self.input_box.add(output_info_box)
        
        self.input_box.add(save_file_button)
        
       


        #----------------------------------------------------------------------
        #the content of the rest steps are defined in the following functions.
        #-----------------------------------------------------------------------
        
        #define a global variable for collecting the PFAS species selected by the user in System_Infor, which will be used in the next two steps.
        self.selected_PFAS=[]
        #define a global variable for containing the optimized parameters obtained in the section - Parameter Calibration, and will be used in the Section: DT_Simulation.
        self.theta_opt_dict={}
        
        #loading png file for step 3 and 4.
        mass_flow_png = resources.files("pfas_dt_investigator.resources").joinpath("mass_flow.png") #used in native app
        DT_analysis_png=resources.files("pfas_dt_investigator.resources").joinpath("DT_analysis.png")
        
        mass_flow_image = toga.Image(src=mass_flow_png) #path or src
        self.mass_flow_image_view = toga.ImageView(image=mass_flow_image, style=Pack(width=400, height=120, margin=10))
        
        DT_analysis_image = toga.Image(src=DT_analysis_png) #path or src
        self.DT_analysis_image_view = toga.ImageView(image=DT_analysis_image, style=Pack(width=400, height=120, margin=10))
        
        
        
       
        # --- Header with navigation buttons ---
        self.header = toga.Box(style=Pack(direction=ROW, margin=5, background_color="white"))
        self.header.add(toga.Button("Input and Output File Setup", on_press=self.file_path_setup))
        self.header.add(toga.Button("System Infor", on_press=self.show_system_infor))
        self.header.add(toga.Button("Mass Flow", on_press=self.Mass_Flow))
        self.header.add(toga.Button("DT_Analysis", on_press=self.DT_analysis))
        self.header.add(toga.Button("Parameter Calibration", on_press=self.Para_Calibration))
        self.header.add(toga.Button("DT_Simulation", on_press=self.DT_Simulation))
        
        




        # Combine header and content
        self.container_step1 = toga.Box(style=Pack(direction=COLUMN))
        #self.container.clear()
        self.container_step1.add(self.header)
        
        #the default content is the first step.
        self.container_step1.add(self.input_box) 
        
        self.main_window.content = self.container_step1
        self.main_window.show()





#----------------------------------------------
#define functions for the self-defined menu
#----------------------------------------------
    def show_about(self, widget):
        self.main_window.info_dialog(
            "About the author",
            "XIANGUI HUANG, a postdoc at Zuckerberg Institute for Water Research, Ben-Gurion University, Israel since APR-2026,\nunder the supervision of Dr. Avner Ronen.\n"
            )











        
#---------------------------------
#Define functions for step 1
#--------------------------------------
    async def on_mass_balance(self, widget):
        await self.main_window.dialog(
            toga.info_dialog(
            "Mass balance calculation triggered", "Hi there!",
            )
        )


    async def error_infor_step1(self, widget):
        await self.main_window.dialog(
            toga.info_dialog(
            "error!, please check the inputdata file path.",
            )
        )
        
        
    
    def load_csv(self, widget):
        if widget.id == 'inputdata':
            path = self.path_input.value.strip()
            try:
               df = pd.read_csv(path, header=0, index_col=0)
               rows, cols = df.shape
               self.rows_label.text = f"Rows: {rows}"
               self.cols_label.text = f"Columns: {cols}"
               self.flag_label.text='A valid input file path is received.'
              
            except Exception as e:
               self.rows_label.text = "Error loading file"
               self.cols_label.text = str(e)
        
        if widget.id == 'inputdata_1':
            path = self.path_input_1.value.strip()
            try:
               df = pd.read_csv(path, header=0, index_col=0)
               rows, cols = df.shape
               self.rows_label_1.text = f"Rows: {rows}"
               self.cols_label_1.text = f"Columns: {cols}"
               self.flag_label_1.text='A valid input file path is received.'
              
            except Exception as e:
               self.rows_label_1.text = "Error loading file"
               self.cols_label_1.text = str(e)
        
        
            
        if widget.id == 'libfile':
            path = self.libfile_path_input.value.strip()
            try:
               df = pd.read_csv(path, header=0, index_col=0)
               rows, cols = df.shape
               self.rows_label_2.text = f"Rows: {rows}"
               self.cols_label_2.text = f"Columns: {cols}"
               self.flag_label_2.text='A valid input file path is received.'
            except Exception as e:
               self.rows_label_2.text = "Error loading file"
               self.cols_label_2.text = str(e)
        
           

    def check_folder(self, widget):
        folder_path=self.output_path_input.value.strip()
        folder_path=Path(folder_path)
        if folder_path.exists() and folder_path.is_dir():
            self.foldercheck_label.text='A valid folder path is received.'
            
        else:
            self.foldercheck_label.text='A invalid folder path is received.'






    def save_file_folder(self, widget):
        inputdata_file_path = self.path_input.value.strip()
        inputdata_file_path=Path(inputdata_file_path)
        
        inputdata_file_path_1=self.path_input_1.value.strip()
        inputdata_file_path_1=Path(inputdata_file_path_1)
        
        folder_path=self.output_path_input.value.strip()
        folder_path=Path(folder_path)
        
        data={'inputdata_file_path': str(inputdata_file_path), 'inputdata_file_path_1': str(inputdata_file_path_1), 'output_folder_path':str(folder_path)}
        with open(CONFIG_FILE, "w") as f:
            json.dump(data, f)
        
        self.main_window.info_dialog("Successfully saved: \n", f"{inputdata_file_path}, \n {inputdata_file_path_1}, \n and {folder_path}")


    
    def load_config(self):
        if CONFIG_FILE.exists():
            try:
                with open(CONFIG_FILE, "r") as f:
                    return json.load(f)
            except json.JSONDecodeError:
                return {}
        return {}




    def file_path_setup(self, widget):
       #step 1. enter input and output file path/folder path.
       self.container_step1 = toga.Box(style=Pack(direction=COLUMN))
       
       self.container_step1.add(self.header)
       
       self.container_step1.add(self.input_box)
       
       self.main_window.content = self.container_step1
       self.main_window.show()




#--------------------------------
#define function for step 2:
#-------------------------------- 

    def show_system_infor(self, widget):
        #step 2 system information dispaly.
         self.system_infor_box=toga.SplitContainer(style=Pack(flex=1, margin=10))
         #since all the file path is valid, let's extract the infor for step2 in this program.
         
         #define the content in the left panel.
         input_file= self.path_input.value.strip() #receive the input file path from step 1.
         file_path=Path(input_file)
         if file_path.exists():
             #-------------------------------
             #Left: show PFAS species list
             #----------------------------------
                          
             PFAS_species_dir_2=FunctionM_1.PFAS_species_reader(input_file)
             
             self.PFAS_species_list = []
             index=0
             for key_i in list(PFAS_species_dir_2.keys()):
                 for specie_i in PFAS_species_dir_2[key_i]:
                     index +=1
                     self.PFAS_species_list += [ [index, specie_i, key_i] ]
             
            #PFAS species table
             #table = toga.Table( headings=['No', "PFAS", "type"], data=PFAS_species_list, style=Pack(flex=3, margin_right=10))
     
            
            # add Checkbox list ---
             checkbox_box = toga.Box(style=Pack(direction=COLUMN, margin=5, flex=1))
             checkbox_box.add(toga.Label("Select species for analysis:", style=Pack(margin_bottom=5, color="blue", font_weight="bold")))
     
             self.switch_list = []
             for i in self.PFAS_species_list:
                  species_i=str(i[0])+'_'+i[1]+'_'+i[2]
                  cb = toga.Switch(species_i, value=False, style=Pack(margin=3))
                  self.switch_list.append(cb)
                  checkbox_box.add(cb)
                  
                  
             #add select all and selection confirm button.
             self.select_all_button = toga.Button("Select all species", on_press=self.select_all, style=Pack(margin=10))
             checkbox_box.add(self.select_all_button)
             
             self.confirm_button = toga.Button("Confirm Selection", on_press=self.PFAS_selection_confirm, style=Pack(margin=10, color='red', font_weight="bold"))
             checkbox_box.add(self.confirm_button)
             
             
             left_container = toga.ScrollContainer(horizontal=False)   
             left_container.content=checkbox_box
             
             
             
             #---------------------------------------------------------------
             #Right: show workflow table and basic infor of each compartment.
             #----------------------------------------------------------------
             #define the content in the right panel, which consists of two boxes, one shows the system workflow;
             #another one display the basic infor of a selected compartment.
             right_box=toga.Box(style=Pack(direction=COLUMN, margin=10))
             
             
             #Right-1: First to show the system workflow.
             self.system_infor_dir_2=FunctionM_1.system_infor_reader(input_file)
             system_diagram_dir=FunctionM_1.PFAS_massflow_network(self.system_infor_dir_2, [self.PFAS_species_list[0][1]], input_file)
             system_diagram_dir_keys=list(system_diagram_dir.keys())
             
             from_list=system_diagram_dir[system_diagram_dir_keys[0]]
             to_list=system_diagram_dir[system_diagram_dir_keys[1]]
             massflow_list=system_diagram_dir[system_diagram_dir_keys[2]]
             massflow_SD_list=system_diagram_dir[system_diagram_dir_keys[3]]
             
             chain_dict=FunctionM_1.chain_seperation(from_list, to_list, massflow_list, massflow_SD_list)
             chain_list=list(chain_dict.keys())
             
             system_workflow_box = toga.Box(style=Pack(direction=ROW, margin=10, flex=1))
             for chain_i in chain_list:
                 chain_i_box= toga.Box(style=Pack(direction=COLUMN, margin=10, flex=1))
                 chain_i_dict=chain_dict[chain_i]
                 chain_i_from_list=chain_i_dict['from']
                 chain_i_to_list=chain_i_dict['to']
                              
                 from_to_list=[]
                 index=0
                 for j in range(len(chain_i_from_list)):
                     index +=1
                     from_to_list +=[ [index, chain_i_from_list[j], chain_i_to_list[j]]]
             
                 #system workflow table
                 chain_i_table = toga.Table( headings=['Edge_No', "From", "To"], data=from_to_list, style=Pack(flex=3, margin_right=10, font_weight="bold"))
                 chain_i_box.add(chain_i_table)
                 system_workflow_box.add(chain_i_box)
            
            
             #-----------------------------------------------------------
             #add the content of Right-2 into the right_box.
             #--------------------------------------------------------------
             right_box.add(system_workflow_box)
             
             
            
             #Right 2: define the content for the bottom box in the right container.
             compt_infor_box = toga.Box(style=Pack(direction=COLUMN, margin=10, flex=1))
             compt_infor_box.add(toga.Label("Check the basic information of each compartment here.", style=Pack(margin_bottom=5, font_weight="bold")))
             
             compt_button_box=toga.Box(style=Pack(direction=ROW, margin=2, height=50))
             compt_list=list(self.system_infor_dir_2.keys())
             for compt_i in compt_list:
                 compt_i_button = toga.Button(compt_i, id= compt_i, on_press=self.show_compt_infor, style=Pack(margin=10, font_weight="bold"))
                 compt_button_box.add(compt_i_button)
                 
             compt_infor_box.add(compt_button_box)
             
            
             compt_para_box=toga.Box(style=Pack(direction=COLUMN, margin=5, flex=1))
             self.para_list=list(self.system_infor_dir_2[compt_list[0]].keys())
             
             # A dict to store Label widgets, keyed by name
             self.label_dict_2 = {}
             for para_i in self.para_list:
                 label_i = toga.Label(f"{para_i}: N/A", style=Pack(margin=(5, 0)))
                 self.label_dict_2[para_i] = label_i  # store label reference
                 compt_para_box.add(label_i)
            
            
             compt_infor_box.add(compt_para_box) 
                             
             
             #-----------------------------------------------------------
             #add the content of Right-2 into the right_box.
             #--------------------------------------------------------------
             right_box.add(compt_infor_box)
             
                     
                      
             #add the right and left panels into the main box.
             self.system_infor_box.content = [(left_container, 1), (right_box, 3)]
             

             self.container_step2 = toga.Box(style=Pack(direction=COLUMN))
             
             self.container_step2.add(self.header)
             self.container_step2.add(self.system_infor_box)
             
             # self.container.clear()
             # self.container.add(self.header)
             # self.container.add(self.system_infor_box)
             
             
             self.main_window.content = self.container_step2
             self.main_window.show()
        
        
         else:
             self.main_window.info_dialog("Error", "An invalid inputdata file path was assigned.")  
         
         
         
         




    def select_all(self, widget):
        self.selected_PFAS=[list_i[1] for list_i in self.PFAS_species_list]
        number_species=len(self.selected_PFAS)
        self.main_window.info_dialog("Selected PFAS", "You selected: {0} species.".format(number_species))



    def PFAS_selection_confirm(self, widget):
        # Get selected PFAS
        self.selected_PFAS = [sw.text for sw in self.switch_list if sw.value]
        if self.selected_PFAS:
            number_species=len(self.selected_PFAS)
            self.main_window.info_dialog("Selected PFAS", "You selected: {0} species.".format(number_species))
        else:
            self.main_window.info_dialog("Selected PFAS", "No PFAS selected.")
    
    
    def show_compt_infor(self, widget):
        compt_i =widget.id
        compt_i_infor_dict=self.system_infor_dir_2[compt_i]
             
        for key, value in compt_i_infor_dict.items():
            if key in self.label_dict_2:
                self.label_dict_2[key].text = f"{key}:  {value}"
        




#-----------------------------------------------------------------
#defing the functions and content for Step 3.
#show the PFAS mass flow at each edges in the system workflow.
#-----------------------------------------------------------------  
    def Mass_Flow(self, widget):
        #step3: Mass Balance box
        #create the box to contain the content for PFAS DT analysis.
        self.mass_flow_analysis_box=toga.Box(style=Pack(direction=COLUMN, margin=5, flex=1))
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        file_path=Path(input_file)
        
        
        if file_path.exists():
            
            #step3-1: add a select widget allowing user to select a species for the following analysis.
            species_selection_box=toga.Box(style=Pack(direction=ROW, margin=5, margin_bottom=10, align_items=CENTER))
            
            self.PFAS_species_dir_3=FunctionM_1.PFAS_species_reader(input_file)
            PFAS_type_list=list(self.PFAS_species_dir_3.keys())
            self.selected_species_3=self.PFAS_species_dir_3[PFAS_type_list[0]][0] #select the first species as the default.
            
            
            self.PFAS_selection_3 = toga.Selection(items=self.selected_PFAS + ['SUM'], style=Pack(margin=10))
            self.PFAS_selection_3.on_change = self.update_massflow_table
            
            
            heading_label=toga.Label("Select a PFAS species: ", style=Pack(margin=5, color="black", font_weight="bold"))
            self.selection_confrim_label_3=toga.Label(self.selected_species_3, style=Pack(margin=5, color="black", font_weight="bold"))
            
            species_selection_box.add(heading_label)
            species_selection_box.add(self.PFAS_selection_3)
            species_selection_box.add(self.selection_confrim_label_3)
            
            
            self.mass_flow_analysis_box.add(species_selection_box)
            self.mass_flow_analysis_box.add(self.mass_flow_image_view)
            
            
            
            #step3-2: show the PFAS mass flow for each edge at each tree in a box.
            mass_flow_box=  toga.Box(style=Pack(direction=ROW, background_color="white", margin=5, flex=0, align_items=CENTER))
            
            system_workflow_box1=toga.Box(style=Pack(direction=COLUMN, margin=5, flex=0, alignment="left"))
            system_workflow_box2=toga.Box(style=Pack(direction=COLUMN, margin=5, flex=0, alignment="left"))
            system_workflow_box3=toga.Box(style=Pack(direction=COLUMN, margin=5, flex=0, alignment="left"))
            system_workflow_box4=toga.Box(style=Pack(direction=COLUMN, margin=5, flex=0, alignment="left"))
            
            PFAS_mass_box=toga.Box(style=Pack(direction=COLUMN, margin=5, alignment="left"))
            PFAS_mass_SD_box=toga.Box(style=Pack(direction=COLUMN, margin=5, alignment="left"))
            
            self.system_infor_dir_3=FunctionM_1.system_infor_reader(input_file)
            system_diagram_dir=FunctionM_1.PFAS_massflow_network(self.system_infor_dir_3, [self.selected_species_3], input_file)
            system_diagram_dir_keys=list(system_diagram_dir.keys())
            
            from_list=system_diagram_dir[system_diagram_dir_keys[0]]
            to_list=system_diagram_dir[system_diagram_dir_keys[1]]
            massflow_list=system_diagram_dir[system_diagram_dir_keys[2]]
            massflow_SD_list=system_diagram_dir[system_diagram_dir_keys[3]]
            
            
            
            chain_dict=FunctionM_1.chain_seperation(from_list, to_list, massflow_list, massflow_SD_list)
            chain_list=list(chain_dict.keys())
            from_to_list=[['Tree No', 'Edge_No', "From", "To"]]
            mass_flow_list=[[f'{system_diagram_dir_keys[2]}', 'SD']]
            for chain_i in chain_list:
                             
                chain_i_dict=chain_dict[chain_i]
                chain_i_from_list=chain_i_dict['from']
                chain_i_to_list=chain_i_dict['to']
                chain_i_massflow_list=chain_i_dict['PFAS_massflow']
                chain_i_massflow_SD_list=chain_i_dict['PFAS_massflow_SD']
                             
                
                index=0
                for j in range(len(chain_i_from_list)):
                    index +=1
                    from_to_list +=[ (chain_i, index, chain_i_from_list[j], chain_i_to_list[j])]
                    mass_flow_list += [[chain_i_massflow_list[j], chain_i_massflow_SD_list[j]] ]
                
                from_to_list +=[("", "", "", "")]
                mass_flow_list +=[["", ""]]
                
           
            
           #creat labels to show the system workflow edges.
            for i, sublist in enumerate(from_to_list):
               Tree_label=toga.Label(sublist[0], style=Pack(margin=6, color="black", font_weight="bold"))
               edge_label=toga.Label(sublist[1], style=Pack(margin=6, color="black", font_weight="bold"))
               from_label=toga.Label(sublist[2], style=Pack(margin=6, color="black", font_weight="bold"))
               to_label=toga.Label(sublist[3], style=Pack(margin=6, color="black", font_weight="bold"))
               
               system_workflow_box1.add(Tree_label)
               system_workflow_box2.add(edge_label)
               system_workflow_box3.add(from_label)
               system_workflow_box4.add(to_label)
               
              
            #creat labels to show the mass flow and SD.
            self.mass_flow_dir={}
            for j, sublist_j in enumerate(mass_flow_list):
                key_massflow=f"massflow_{j}"
                key_massflow_SD=f"massflow_SD_{j}"
                mass_flow_j=sublist_j[0]
                massflow_SD_j=sublist_j[1]
                
                massflow_label=toga.Label(mass_flow_j, style=Pack(margin=6, color="blue", font_weight="bold"))
                massflow_SD_label=toga.Label(massflow_SD_j, style=Pack(margin=6, color="blue", font_weight="bold"))
                
                self.mass_flow_dir[key_massflow]=massflow_label
                self.mass_flow_dir[key_massflow_SD]=massflow_SD_label
                
                PFAS_mass_box.add(massflow_label)
                PFAS_mass_SD_box.add(massflow_SD_label)
                
            
            mass_flow_box.add(system_workflow_box1)
            mass_flow_box.add(system_workflow_box2)
            mass_flow_box.add(system_workflow_box3)
            mass_flow_box.add(system_workflow_box4)
            
            mass_flow_box.add(PFAS_mass_box)
            mass_flow_box.add(PFAS_mass_SD_box)
            
            # creat a Outer box (blue border box) to creat a boundary line for the massflow box.
            outer_box = toga.Box(children=[mass_flow_box], style=Pack(margin=3, background_color="gray", alignment=CENTER, direction=ROW, flex=0) ) 
            
            
            self.mass_flow_analysis_box.add(outer_box)


            
            #step 3-3: add a button to allow the user to save all the mass flow including all the PFAS species and all the edges in the system.
            save_file_box=toga.Box(style=Pack(direction=ROW, margin=10))
            self.massflow_file_name=toga.TextInput(placeholder='Enter the csv file name ...', style=Pack(width=200, margin=10))
            save_all_massflow_analysis = toga.Button('Save all the mass flow analysis results into a csv file', on_press=self.save_massflow_analysis, style=Pack(margin=10, color= "black", font_weight="bold"))
            
            save_file_box.add(self.massflow_file_name)
            save_file_box.add(save_all_massflow_analysis)
            
            self.mass_flow_analysis_box.add(save_file_box)

    
        
        
        
        else:
            self.main_window.info_dialog("Error", "An invalid inputdata file path was assigned.")
                
        
        self.container_step3 = toga.Box(style=Pack(direction=COLUMN))
        
        self.container_step3.add(self.header)
        self.container_step3.add(self.mass_flow_analysis_box)
        
        # self.container.clear()
        # self.container.add(self.header)
        # self.container.add( self.mass_flow_analysis_box)
                
        self.main_window.content = self.container_step3
        self.main_window.show()
       
        
        
        
        
    def update_massflow_table(self, widget):
        selected = widget.value
        self.selected_species_3 = selected
        self.selection_confrim_label_3.text = f"Selected species: {selected}"
        
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        if self.selected_species_3 == 'SUM':
            system_diagram_dir=FunctionM_1.PFAS_massflow_network(self.system_infor_dir_3, self.selected_PFAS, input_file)
        else:
            system_diagram_dir=FunctionM_1.PFAS_massflow_network(self.system_infor_dir_3, [self.selected_species_3], input_file)
                
        #updated the mass flow table.                
        system_diagram_dir_keys=list(system_diagram_dir.keys())
        
        from_list=system_diagram_dir[system_diagram_dir_keys[0]]
        to_list=system_diagram_dir[system_diagram_dir_keys[1]]
        massflow_list=system_diagram_dir[system_diagram_dir_keys[2]]
        massflow_SD_list=system_diagram_dir[system_diagram_dir_keys[3]]
        
        chain_dict=FunctionM_1.chain_seperation(from_list, to_list, massflow_list, massflow_SD_list)
        chain_list=list(chain_dict.keys())
        
        mass_flow_list=[[f'{system_diagram_dir_keys[2]}', 'SD']]
        for chain_i in chain_list:
                         
            chain_i_dict=chain_dict[chain_i]
            
            chain_i_massflow_list=chain_i_dict['PFAS_massflow']
            chain_i_massflow_SD_list=chain_i_dict['PFAS_massflow_SD']
                                                 
            for j in range(len(chain_i_massflow_list)):
                mass_flow_list += [[chain_i_massflow_list[j], chain_i_massflow_SD_list[j]] ]
                        
            mass_flow_list +=[["", ""]]
        
        for j, sublist_j in enumerate(mass_flow_list):
            key_massflow=f"massflow_{j}"
            key_massflow_SD=f"massflow_SD_{j}"
            mass_flow_j=sublist_j[0]
            massflow_SD_j=sublist_j[1]
            
            self.mass_flow_dir[key_massflow].text=mass_flow_j
            self.mass_flow_dir[key_massflow_SD].text=massflow_SD_j
        
        
    def save_massflow_analysis(self, widget):
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        
        all_massflow_result_dict={}
        PFAS_list=[]
        from_list=[]
        to_list=[]
        massflow_list=[]
        massflow_SD_list=[]
        
        key_list=[]
        
        for i in self.selected_PFAS+['SUM']:
            
            if i != 'SUM':
                PFAS_i=i
                system_diagram_dir=FunctionM_1.PFAS_massflow_network(self.system_infor_dir_3, [PFAS_i], input_file)
            else:
                system_diagram_dir=FunctionM_1.PFAS_massflow_network(self.system_infor_dir_3, self.selected_PFAS, input_file)
                
            
            
            system_diagram_dir_keys=list(system_diagram_dir.keys())
            key_list=system_diagram_dir_keys
            
            from_list_i=system_diagram_dir[system_diagram_dir_keys[0]]
            to_list_i=system_diagram_dir[system_diagram_dir_keys[1]]
            massflow_list_i=system_diagram_dir[system_diagram_dir_keys[2]]
            massflow_SD_list_i=system_diagram_dir[system_diagram_dir_keys[3]]
            
            if i != 'SUM': PFAS_list_i=[i]*len(from_list_i)
            else: PFAS_list_i=["SUM"]*len(from_list_i)
            
            PFAS_list +=PFAS_list_i
            from_list +=from_list_i 
            to_list += to_list_i
            massflow_list += massflow_list_i
            massflow_SD_list +=massflow_SD_list_i
      
            
        all_massflow_result_dict['PFAS'] = PFAS_list
        all_massflow_result_dict[key_list[0]] = from_list 
        all_massflow_result_dict[key_list[1]] = to_list      
        all_massflow_result_dict[key_list[2]] = massflow_list
        all_massflow_result_dict[key_list[3]] = massflow_SD_list
        
        all_massflow_result_df=pd.DataFrame(all_massflow_result_dict, index=list(range(len(PFAS_list))) )
        
        output_folder_path=Path( self.output_path_input.value.strip() )
        output_filename=self.massflow_file_name.value.strip()
        output_file_path=output_folder_path/output_filename
        
        all_massflow_result_df.to_csv(output_file_path, index=True)
        
        self.main_window.info_dialog("Successfully saved", f"Saved PFAS DT results into: {output_file_path}")
            
            
            
        
        
        

        
       
       
   
#----------------------------------------------------------------
#define the actions for step 4:
    #step4: PFAS distribution and transformation (DT) analysis
#---------------------------------------------------------------        
   
    def DT_analysis(self, widget):
        #create the box to contain the content for PFAS DT analysis.
        self.DT_analysis_box=toga.Box(style=Pack(direction=COLUMN, margin=5))
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        file_path=Path(input_file)
        if file_path.exists():
            #step4-1: add a select widget allowing user to select a species for the following analysis.
            species_selection_box=toga.Box(style=Pack(direction=ROW, margin=5, margin_bottom=10, align_items=CENTER))
            
            self.PFAS_species_dir_4=FunctionM_1.PFAS_species_reader(input_file)
            PFAS_type_list=list(self.PFAS_species_dir_4.keys())
            self.selected_species=self.PFAS_species_dir_4[PFAS_type_list[0]][0] #select the first species as the default.
            
            self.PFAS_selection = toga.Selection(items=self.selected_PFAS, style=Pack(margin=10))
            self.PFAS_selection.on_change = self.on_select_species
            
            
            self.heading_label=toga.Label("Select a PFAS species: ", style=Pack(margin=5, color="black", font_weight="bold"))
            self.selection_confrim_label=toga.Label(self.selected_species, style=Pack(margin=5, color="black", font_weight="bold"))
            
            species_selection_box.add(self.heading_label)
            species_selection_box.add(self.PFAS_selection)
            species_selection_box.add(self.selection_confrim_label)
            
            
            #step4-2: add a select widget allowing user decide wheather consider standard deviation (SD) in the following analysis.
            self.consider_SD_selection = toga.Selection(items=[True, False], style=Pack(margin=10))
            self.consider_SD_selection.on_change = self.SD_consideration
            self.heading_label_1=toga.Label("Consider PFAS analysis SD: ", style=Pack(margin_left=20, color="Green", font_weight="bold"))
            self.SD_confrim_label=toga.Label(str(FunctionM_1.consider_SD_PFASmass_CalC), style=Pack(margin=5, color="Green", font_weight="bold"))
            
            species_selection_box.add(self.heading_label_1)
            species_selection_box.add(self.consider_SD_selection)
            species_selection_box.add(self.SD_confrim_label)
            
            
            self.DT_analysis_box.add(species_selection_box)
            self.DT_analysis_box.add(self.DT_analysis_image_view)
            
            
            
            
            
            #step4-3: perform the DT analysis based on the selected PFAS species.
            # creat a Outer box (blue border box) to creat a boundary line for the massflow box.
            PFAS_DT_inner_box = toga.Box(style=Pack(margin=5, background_color="white", direction=COLUMN, flex=0) )
        
            self.system_infor_dir_4=FunctionM_1.system_infor_reader(input_file)
            self.compt_list=list(self.system_infor_dir_4.keys())
           
            self.label_dict_4 = {}
            for compt_i in self.compt_list:
                compt_i_box=toga.Box(style=Pack(direction=ROW, margin=1, background_color="white", align_items=CENTER))
                
                PFAS_DT_compt_i=FunctionM_1.DT_analyzer(self.system_infor_dir_4, [self.selected_species], [compt_i], self.PFAS_species_dir_4, input_file)
                keys_list=list(PFAS_DT_compt_i.keys())
                compt_i_label=self.system_infor_dir_4[compt_i]['label']
                input_mass= PFAS_DT_compt_i[keys_list[2]][0]
                if not isinstance(input_mass, str):
                    input_mass=round(input_mass, 3)
                output_mass=PFAS_DT_compt_i[keys_list[3]][0]
                if not isinstance(output_mass, str):
                    output_mass=round(output_mass, 3)
                mass_change=PFAS_DT_compt_i[keys_list[4]][0]
                reason=PFAS_DT_compt_i[keys_list[5]][0]
                
                DT_results={keys_list[2]: input_mass, compt_i: compt_i_label, keys_list[3]:output_mass, keys_list[4]: mass_change, keys_list[5]: reason}
                color_list=['blue', 'black', 'blue', 'red', 'red']
                for i, key_i in enumerate(list(DT_results.keys())):
                    Label_name=compt_i + '_' + key_i
                    Label_name_index=compt_i + '_'+key_i+'_index'
                    if key_i != 'reason':
                        Label_name=toga.Label(f'{key_i}: {DT_results[key_i]} >>', style=Pack(margin=5, color=color_list[i], font_weight="bold"))
                    else:
                        Label_name=toga.Label(f'{key_i}: {DT_results[key_i]}', style=Pack(margin=5, color=color_list[i], font_weight="bold"))
                    
                    compt_i_box.add(Label_name)
                    self.label_dict_4[Label_name_index]=Label_name
                    
                                         
                # Outer box acts as the border
                bordered_box = toga.Box(children=[compt_i_box], style=Pack(margin=0.5, background_color="grey", margin_bottom=40, alignment=CENTER))
                     
                PFAS_DT_inner_box.add(bordered_box)
                
            
            PFAS_DT_outer_box = toga.Box(children=[PFAS_DT_inner_box], style=Pack(margin=2, background_color="grey", direction=ROW, flex=0) )
            self.DT_analysis_box.add(PFAS_DT_outer_box)
              
              
              
            #Step4-4: put a button, allowing user to save all the DT analysis results into a csv file.
            save_file_box=toga.Box(style=Pack(direction=ROW, margin=10))
            self.DT_file_name=toga.TextInput(placeholder='Enter the csv file name ...', style=Pack(width=200, margin=10))
            save_all_DT_analysis = toga.Button('Save all the DT analysis results into a csv file', on_press=self.save_DT_analysis, style=Pack(margin=10, color= "black", font_weight="bold"))
            
            save_file_box.add(self.DT_file_name)
            save_file_box.add(save_all_DT_analysis)
            
            self.DT_analysis_box.add(save_file_box)
              
                        
        else:
            self.main_window.info_dialog("Error", "An invalid inputdata file path was assigned.")  
        
                      
        self.container_step4 = toga.Box(style=Pack(direction=COLUMN))
        
        self.container_step4.add(self.header)
        self.container_step4.add(self.DT_analysis_box)
                
        
        # self.container.clear()
        # self.container.add(self.header)
        # self.container.add(self.DT_analysis_box)
                
        self.main_window.content = self.container_step4
        self.main_window.show()
    

    

    def on_select_species(self, widget):
        selected = widget.value
        self.selected_species = selected
        self.selection_confrim_label.text = f"Selected species: {selected}"
        
        
        #update the DT results in the following box.
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        for compt_i in self.compt_list:                       
            PFAS_DT_compt_i=FunctionM_1.DT_analyzer(self.system_infor_dir_4, [self.selected_species], [compt_i], self.PFAS_species_dir_4, input_file)
            keys_list=list(PFAS_DT_compt_i.keys())
            
            compt_i_label=self.system_infor_dir_4[compt_i]['label']
            
            input_mass= PFAS_DT_compt_i[keys_list[2]][0]
            if not isinstance(input_mass, str):
                input_mass=round(input_mass, 3)
                
            output_mass=PFAS_DT_compt_i[keys_list[3]][0]
            if not isinstance(output_mass, str):
                output_mass=round(output_mass, 3)
                
            mass_change=PFAS_DT_compt_i[keys_list[4]][0]
            
            reason=PFAS_DT_compt_i[keys_list[5]][0]
            
            DT_results={keys_list[2]: input_mass, compt_i: compt_i_label, keys_list[3]:output_mass, keys_list[4]: mass_change, keys_list[5]: reason}
            for i, key_i in enumerate(list(DT_results.keys())):
                Label_name_index=compt_i + '_'+key_i+'_index'
                if Label_name_index in self.label_dict_4:
                    if key_i != 'reason':
                        self.label_dict_4[Label_name_index].text=f'{key_i}: {DT_results[key_i]} >>'
                    else:
                        self.label_dict_4[Label_name_index].text=f'{key_i}: {DT_results[key_i]}'
    
 
    

    def SD_consideration(self, widget):
        FunctionM_1.consider_SD_PFASmass_CalC
        FunctionM_1.consider_SD_PFASmass_CalC=widget.value
        
        self.SD_confrim_label.text=str(FunctionM_1.consider_SD_PFASmass_CalC)
        
        #update the DT results in the following box.
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        for compt_i in self.compt_list:                       
            PFAS_DT_compt_i=FunctionM_1.DT_analyzer(self.system_infor_dir_4, [self.selected_species], [compt_i], self.PFAS_species_dir_4, input_file)
            keys_list=list(PFAS_DT_compt_i.keys())
            
            compt_i_label=self.system_infor_dir_4[compt_i]['label']
            
            input_mass= PFAS_DT_compt_i[keys_list[2]][0]
            if not isinstance(input_mass, str):
                input_mass=round(input_mass, 3)
                
            output_mass=PFAS_DT_compt_i[keys_list[3]][0]
            if not isinstance(output_mass, str):
                output_mass=round(output_mass, 3)
                
            mass_change=PFAS_DT_compt_i[keys_list[4]][0]
            
            reason=PFAS_DT_compt_i[keys_list[5]][0]
            
            DT_results={keys_list[2]: input_mass, compt_i: compt_i_label, keys_list[3]:output_mass, keys_list[4]: mass_change, keys_list[5]: reason}
            for i, key_i in enumerate(list(DT_results.keys())):
                Label_name_index=compt_i + '_'+key_i+'_index'
                if Label_name_index in self.label_dict_4:
                    if key_i != 'reason':
                        self.label_dict_4[Label_name_index].text=f'{key_i}: {DT_results[key_i]} >>'
                    else:
                        self.label_dict_4[Label_name_index].text=f'{key_i}: {DT_results[key_i]}'



    def save_DT_analysis(self, widget):
        PFAS_dict_value=list(self.PFAS_species_dir_4.values())
        PFAS_list=[item for sublist in PFAS_dict_value for item in sublist]
        input_file= self.path_input.value.strip() #receive the input file path from step 1.
        
        all_PFAS_DT_result_dict=FunctionM_1.DT_analyzer(self.system_infor_dir_4, PFAS_list, self.compt_list, self.PFAS_species_dir_4, input_file)
        all_PFAS_DT_result_df=pd.DataFrame(all_PFAS_DT_result_dict, index=list(range(len(PFAS_list)*len(self.compt_list))) )
        
        
        output_folder_path=Path( self.output_path_input.value.strip() )
        output_filename=self.DT_file_name.value.strip()
        output_file_path=output_folder_path/output_filename
        
        all_PFAS_DT_result_df.to_csv(output_file_path, index=True)
        
        self.main_window.info_dialog("Successfully saved", f"Saved PFAS DT results into: {output_file_path}")
        
        





#---------------------------------------------------------------------------------------------
#define the actions for calibration:
    #Calibration: performe the parameter calibration based on observation for
    #the simplified model, in which these parameter is used to quantify the PFAS DT in WWTP. 
#--------------------------------------------------------------------------------------------
    def Para_Calibration(self, widget):
        #create the box to contain the content for PFAS DT analysis.
        self.Para_Calibration_box=toga.Box(style=Pack(direction=COLUMN, margin=5))
        
        input_file= self.path_input_1.value.strip() #receive the input file path from step 1.
        file_path=Path(input_file)
        if file_path.exists():
            #The windonw showing the calibration consists of three boxs in column: top, middle, and bottom boxes.
            
            #box 1: The top box, containing interface allowing user to select compartment, and run the calibration.
            self.heading_label=toga.Label("Select a compartment: ", style=Pack(margin=5, color="black", font_weight="bold"))
            self.selected_compt='Clarifier1' #set this as default.
            self.compt_selection = toga.Selection(items=['Clarifier1', 'Reactor_Clarifier2'], style=Pack(margin=10))
            self.compt_selection.on_change = self.on_select_compt    #a function to show the corresponding initial para value and obs in the middle box.       
                        
            self.Calibration_botton=toga.Button('Run Calibration', on_press=self.Run_Calibration, style=Pack(margin_left=15))
            
            compt_selection_box = toga.Box(children=[self.heading_label, self.compt_selection, self.Calibration_botton], style=Pack(direction=ROW, margin=5, align_items=CENTER))
            self.Para_Calibration_box.add(compt_selection_box)
            
            
            #box 2: The middle box, show the initial para value and obs.
            self.head_label=toga.Label("The initial value of parameters (theta_0) and number of observation", style=Pack(margin=10, color="black", font_weight="bold"))
            self.RLS_input = toga.TextInput(placeholder='Change the ratio solid-to-liquid (RSL) if needed. In default: it is 0.0001 for Clarifier1 and 0.001 for Reactor_Clarifier2.', style=Pack(flex=3, margin_bottom=10))
            self.opt_method_input = toga.TextInput(placeholder='Change the optimization method if needed. In default: it is dogbox, avaliable options: dogbox, trf, and lm.', style=Pack(flex=3, margin_bottom=10))
            self.opt_loss_input = toga.TextInput(placeholder='Change the loss function keyword if needed. In default: it is soft_l1, avaliable options: linear, soft_l1, huber, cauchy, and arctan.', style=Pack(flex=3, margin_bottom=10))
            
            self.theta_head_label_1=toga.Label("theta list: [fj, fk, PC_i, PC_j, PC_k]", style=Pack(margin=5, color="black"))
            self.theta_0_label=toga.Label("theta_0: ", style=Pack(margin=5, color="black"))
            self.theta_lb_label=toga.Label("theta_lb: ", style=Pack(margin=5, color="black"))
            self.theta_ub_label=toga.Label("theta_ub: ", style=Pack(margin=5, color="black"))
           
            
            self.X_obs_label=toga.Label("Number of input observation: ", style=Pack(margin=5, color="black"))
            self.Y_obs_label=toga.Label("Number of output observation: ", style=Pack(margin=5, color="black"))
                                       
            initial_para_obs_box=toga.Box(children=[self.head_label, self.RLS_input,self.opt_method_input, self.opt_loss_input, self.theta_head_label_1, self.theta_0_label, self.theta_lb_label, self.theta_ub_label, self.X_obs_label, self.Y_obs_label], style=Pack(direction=COLUMN, margin=5) )                           
            self.Para_Calibration_box.add(initial_para_obs_box)


            
            #box 3: The bottom box, show the optimized para value and basic infor of optimization performance.
            self.opt_head_label=toga.Label("The optimized value of parameters (theta_opt) and optimization performance", style=Pack(margin=10, color="black", font_weight="bold"))
            self.theta_head_label_2=toga.Label("theta list: [fj, fk, PC_i, PC_j, PC_k]", style=Pack(margin=5, color="green"))
            self.theta_opt_label=toga.Label("theta_opt: ", style=Pack(margin=5, color="green"))
                        
            
            self.nef_label=toga.Label("Number of model simulation (nvef): ", style=Pack(margin=5, color="green"))
            self.cost_label=toga.Label("Total residule square (cost): ", style=Pack(margin=5, color="green"))
            self.message_label=toga.Label("Human-readable stop reason: ", style=Pack(margin=5, color="green"))
            
            opt_para_box=toga.Box(children=[self.opt_head_label, self.theta_head_label_2, self.theta_opt_label, self.nef_label, self.cost_label, self.message_label], style=Pack(direction=COLUMN, margin=5) )
            self.Para_Calibration_box.add(opt_para_box)


                       
            
            
        else:
            self.main_window.info_dialog("Error", "An invalid inputdata file path was assigned.")
        
        
        self.container_Calibration = toga.Box(style=Pack(direction=COLUMN))
        
        self.container_Calibration.add(self.header)
        self.container_Calibration.add(self.Para_Calibration_box)
                
        
        # self.container.clear()
        # self.container.add(self.header)
        # self.container.add(self.DT_analysis_box)
                
        self.main_window.content = self.container_Calibration
        self.main_window.show()






    def on_select_compt(self, widget):
        selected = widget.value
        self.selected_compt = selected
        input_file= self.path_input_1.value.strip()
        
        if self.selected_compt == 'Reactor_Clarifier2':
            self.theta_head_label_1.text="theta list: [fj, fk, f_RAS, PC_i, PC_j, PC_k]"
            self.theta_head_label_2.text="theta list: [fj, fk, f_RAS, PC_i, PC_j, PC_k]"
        elif self.selected_compt == 'Clarifier1':
            self.theta_head_label_1.text="theta list: [fj, fk, PC_i, PC_j, PC_k]"
            self.theta_head_label_2.text="theta list: [fj, fk, PC_i, PC_j, PC_k]"
            
        initial_para_obs_dict=FunctionM_2.observation_reader(input_file, self.selected_compt)
        
        self.theta_0_label.text=f"theta_0: {initial_para_obs_dict['theta_0']}"
        self.theta_lb_label.text=f"theta_lb: {initial_para_obs_dict['lb']}"
        self.theta_ub_label.text=f"theta_ub: {initial_para_obs_dict['ub']}"
        
        X_obs_array=initial_para_obs_dict['X']
        Y_obs_array=initial_para_obs_dict['Y_obs']
        
        self.X_obs_label.text=f"Number of input observation: {X_obs_array.shape[0]}"
        self.Y_obs_label.text=f"Number of output observation: {Y_obs_array.shape[0]}"
        
        
        
    def Run_Calibration(self, widget):
        
        input_file= self.path_input_1.value.strip()
                            
        initial_para_obs_dict=FunctionM_2.observation_reader(input_file, self.selected_compt)
        theta_0=initial_para_obs_dict['theta_0']
        lb=initial_para_obs_dict['lb']
        ub=initial_para_obs_dict['ub']
        X_obs_array=initial_para_obs_dict['X']
        X=X_obs_array.tolist()
        Y_obs_array=initial_para_obs_dict['Y_obs']
        Y_obs=Y_obs_array.tolist()
        
        
        #updated the RLS if new value is inputed.
        RLS=self.RLS_input.value.strip()
        method=self.opt_method_input.value.strip()
        loss=self.opt_loss_input.value.strip()
        
        if RLS !="":
            RLS=float(RLS)
            if self.selected_compt == 'Clarifier1':
                FunctionM_2.RSL_Clarifier1=RLS
            elif self.selected_compt == 'Reactor_Clarifier2':
                FunctionM_2.RSL_Reactor_Clarifier2=RLS
        
        if method !="":
            FunctionM_2.optimization_method=method
        if loss !="":
            FunctionM_2.optimization_loss=loss
        
        
        #performe the parameter optimization
        if self.selected_compt == 'Clarifier1':
            res=FunctionM_2.residules_Clarifier1_optimization(theta_0, lb, ub, X, Y_obs)
        elif self.selected_compt == 'Reactor_Clarifier2':
            res=FunctionM_2.residules_Reactor_Clarifier2_optimization(theta_0, lb, ub, X, Y_obs)
        
        self.theta_opt_dict[self.selected_compt] = res['theta_hat'] #store the optimized theta for DT_Simulation.
        
        self.theta_opt_label.text=f"theta_opt: {res['theta_hat']}"
        self.nef_label.text=f"Number of model simulation (nvef): {res['nfev']}"
        self.cost_label.text=f"Total residule square (cost): {2*res['cost']}"
        self.message_label.text=f"Human-readable stop reason: {res['message']}"





#---------------------------------------------------------------------------------------------------------------
#define the actions for DT_simulation:
    #DT_simulation in WWTP: performe the PFAS DT simulation based on the optimized parameters and obsered input data. 
#------------------------------------------------------------------------------------------------------------
    def DT_Simulation(self, widget):
        #create the box to contain the content for PFAS DT analysis.
        self.DT_Simulation_box=toga.Box(style=Pack(direction=COLUMN, margin=5))
        
        #The windonw showing the DT_simulation consists of three boxs in column: top, middle, and bottom boxes.
        
        #box 1: The top box, containing interface allowing user to select compartment, and run the calibration.
        heading_label_1=toga.Label("Select a compartment: ", style=Pack(margin=5, color="black", font_weight="bold"))
        self.selected_compt_simulation='Clarifier1' #set this as default.
        self.compt_selection_simulation = toga.Selection(items=['Clarifier1', 'Reactor_Clarifier2'], style=Pack(margin=10))
        self.compt_selection_simulation.on_change = self.select_compt_simulation    #a function to show the corresponding opt para value and input obs in the middle box.       
                    
        self.simulation_botton=toga.Button('Run Simulation', on_press=self.Run_Simulation, style=Pack(margin_left=15))
        
        compt_selection_box = toga.Box(children=[heading_label_1, self.compt_selection_simulation, self.simulation_botton], style=Pack(direction=ROW, margin=5, align_items=CENTER))
        self.DT_Simulation_box.add(compt_selection_box)
        
        
        #box 2: The middle box contain the interface allow user to enter the input data and theta_opt.
        heading_label_2=toga.Label("Input the parameters (theta_opt) and input observation", style=Pack(margin=10, color="black", font_weight="bold"))
        self.input_paralist_lable=toga.Label("theta list: fj, fk, PC_i, PC_j, PC_k", style=Pack(margin=5, color="black"))
        self.para_input = toga.TextInput(placeholder='Typing the parameters list for simulation. In default: it is the optimized theta obtained from Calibration.', style=Pack(flex=3, margin_bottom=10))
        self.input_Xlist_label=toga.Label("Input X list: V, Cli, Clj, Clk", style=Pack(margin=5, color="black"))
        self.X_input = toga.TextInput(placeholder='Typing the X list for simulation.', style=Pack(flex=3, margin_bottom=10))
        
        data_input_box = toga.Box(children=[heading_label_2, self.input_paralist_lable, self.para_input, self.input_Xlist_label, self.X_input], style=Pack(direction=COLUMN, margin=5))
        self.DT_Simulation_box.add(data_input_box)
        
        
        #box 3: The bottom box show the calculated results.
        self.output_Ylist_label=toga.Label("Output Y list: [V, Cli_out, Clj_out, Clk_out] = ...", style=Pack(margin=5, color="green"))
        self.DT_Simulation_box.add(self.output_Ylist_label)
        
        
        
        self.container_Simulation = toga.Box(style=Pack(direction=COLUMN))
            
        self.container_Simulation.add(self.header)
        self.container_Simulation.add(self.DT_Simulation_box)
                
        
        # self.container.clear()
        # self.container.add(self.header)
        # self.container.add(self.DT_analysis_box)
                
        self.main_window.content = self.container_Simulation
        self.main_window.show()
    
    
  
        
  
    
    def select_compt_simulation(self, widget):
        selected = widget.value
        self.selected_compt_simulation=selected
        
        if self.selected_compt_simulation == 'Reactor_Clarifier2':
            self.input_paralist_lable.text= "theta list: fj, fk, f_RAS, PC_i, PC_j, PC_k"
            self.input_Xlist_label.text= "Input X list: V, Cli, Clj, Clk, Csi, Csj, Csk"
        elif self.selected_compt_simulation == 'Clarifier1':
            self.input_paralist_lable.text= "theta list: fj, fk, PC_i, PC_j, PC_k"
            self.input_Xlist_label.text= "Input X list: V, Cli, Clj, Clk"
        
        theta_list_1=self.para_input.value.strip()
        if theta_list_1 != "":
            try:
                theta_list_2 = [float(x) for x in theta_list_1.split(",")]
            except ValueError:
                self.main_window.info_dialog('Error', "A invalid para list is entered. \nNotice: each element is seperated with commas and no Square Brackets.")
        else:
            if self.selected_compt_simulation in list(self.theta_opt_dict.keys()):
                theta_list_2 = self.theta_opt_dict[self.selected_compt_simulation]
            else:
                self.main_window.info_dialog('Error', 'Required parameter list is neither entered nor generated in Parameter Calibration.')
                  
                
        X_list_1=self.X_input.value.strip()
        try:
            X_list_2=[float(x) for x in X_list_1.split(",")]
        except ValueError:
            self.main_window.info_dialog('Error', "A invalid X list is entered. \nNotice: each element is seperated with commas and no Square Brackets.")
            
        
        #assigning the received parameter list and X list into a global dict variable, that will be used in simulation.
        self.received_para_X_dict={}
        self.received_para_X_dict['theta']=theta_list_2
        self.received_para_X_dict['X']=X_list_2
        
        
        
    def Run_Simulation(self, widget):
        #extract the parameter and X for simualtion.
        if self.selected_compt_simulation == 'Reactor_Clarifier2':
            fj_0, fk_0, f_RAS_0, PC_i_0, PC_j_0, PC_k_0=self.received_para_X_dict['theta']
            V_I, C_i, C_j, C_k, Ci_RAS, Cj_RAS, Ck_RAS = self.received_para_X_dict['X']
            Y_hat_dict=FunctionM_2.PFAS_DT_Ractor_Clarifier2(V_I = V_I, Cl_i= C_i, Cl_j=C_j, Cl_k=C_k, Ci_RAS=Ci_RAS, Cj_RAS=Cj_RAS, Ck_RAS=Ck_RAS, fj = fj_0, fk = fk_0, f_RAS= f_RAS_0, PC_i = PC_i_0, PC_j = PC_j_0, PC_k = PC_k_0)
            Y_hat_list=[Y_hat_dict[i] for i in list(Y_hat_dict.keys()) ]
            self.output_Ylist_label.text=f"Output Y list: [V, Cli_out, Clj_out, Clk_out, Csi_out, Csj_out, Csk_out] = {Y_hat_list}"
            
        elif self.selected_compt_simulation == 'Clarifier1':
            fj_0, fk_0, PC_i_0, PC_j_0, PC_k_0=self.received_para_X_dict['theta']
            V_0, C_i, C_j, C_k = self.received_para_X_dict['X']
            Y_hat_dict=FunctionM_2.PFAS_DT_Clarifier1(V0 = V_0, Ci = C_i, Cj = C_j, Ck = C_k, fj = fj_0, fk = fk_0, PC_i = PC_i_0, PC_j = PC_j_0, PC_k = PC_k_0)
            Y_hat_list=[Y_hat_dict[i] for i in list(Y_hat_dict.keys()) ]
            self.output_Ylist_label.text=f"Output Y list: [V, Cli_out, Clj_out, Clk_out] = {Y_hat_list}"
        
        
        









def main():
    return MyApp("PFAS_DT_Investigator", "org.example.pfas_dt_investigator")

