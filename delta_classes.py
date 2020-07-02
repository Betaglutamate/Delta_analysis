import numpy as np
from more_itertools import split_after, split_before

class Chamber:
    def __init__(self, chamber_nb, cells):
        self.chamber_nb = chamber_nb
        self.cells = cells
    
    def __repr__(self):
        return f"Chamber no: {self.chamber_nb} contains {len(self.cells)} cell objects"

class Cell:
    def __init__(self, experiment_cell_number, cell_number, mother_nb,
                 chamber_nb, framenbs, daughters,
                 length, width, area, fluorescence,
                 family_tree):
        self.experiment_cell_number = experiment_cell_number
        self.cell_number = cell_number
        self.mother_number = mother_nb
        self.chamber_nb = chamber_nb
        self.framenbs = np.asarray(framenbs)
        self.daughters = np.asarray(daughters)
        self.length = np.asarray(length)
        self.width = np.asarray(width)
        self.area = np.asarray(area)
        self.fluorescence = np.asarray(fluorescence)
        self.family_tree = family_tree
    
    def __repr__(self):
        return(f"object number: {self.experiment_cell_number},\
            cell number: {self.cell_number},\
            chamber number: {self.chamber_nb}")
    
    def __len__(self):
        return len(self.framenbs)
    
    time_interval_per_frame = 5

    def make_ancestry_map(self):
        ancestry_map = {}
        for ancestry in self.family_tree:
            ancestry_map[ancestry['cell_nb']] = ancestry['mother_nb']
        self.ancestry_map = ancestry_map
        return self.ancestry_map


    def find_lineage(self):
        cell_nb = self.cell_number
        ancestry = []
        while cell_nb:
            cell_nb = self.ancestry_map[cell_nb]
            ancestry.append(cell_nb)
        self.ancestry = ancestry
        return self.ancestry

    def calculate_fluorescence_by_area(self):
        fluorescence_by_area = [x/y for x, y in zip(self.fluorescence, self.area)]
        self.fluorescence_by_area = np.asarray(fluorescence_by_area)
        return self.fluorescence_by_area


    def calculate_length_growth(self, time_interval_per_frame=time_interval_per_frame):
        length = self.length
        delta_length = [j-i for i, j in zip(length[:-1], length[1:])] 
        dL_dT = [x/time_interval_per_frame for x in delta_length]
        self.length_growth = np.asarray(dL_dT)
        return self.length_growth
    
    def calculate_width_growth(self, time_interval_per_frame=time_interval_per_frame):
        width = self.width
        delta_width = [j-i for i, j in zip(width[:-1], width[1:])] 
        dW_dT = [x/time_interval_per_frame for x in delta_width]
        self.width_growth = np.asarray(dW_dT)
        return self.width_growth
        
    def calculate_area_growth(self, time_interval_per_frame=time_interval_per_frame):
        area = self.area
        delta_area = [j-i for i, j in zip(area[:-1], area[1:])] 
        dA_dT = [x/time_interval_per_frame for x in delta_area]
        self.area_growth = np.asarray(dA_dT)
        return self.area_growth
    
    def calculate_fluorescence_growth(self, time_interval_per_frame=time_interval_per_frame):
        fluorescence = self.fluorescence
        delta_fluorescence = [j-i for i, j in zip(fluorescence[:-1], fluorescence[1:])] 
        dF_dT = [x/time_interval_per_frame for x in delta_fluorescence]
        self.fluorescence_growth = np.asarray(dF_dT)
        return self.fluorescence_growth

    def calculate_fluorescence_growth_by_area(self):
        try:
            dF_dT_byArea = [x/y for x, y in zip(self.fluorescence_growth, self.area)]
            self.fluorescence_growth_by_area = np.asarray(dF_dT_byArea)
            return self.fluorescence_growth_by_area
        except:
             NameError("unable to perform operation did you run function calculate_fluorescence_growth?")


    def create_subCells(self):
        '''This function splits all of the cell data into subcell objects
        to do this I find the length of the list before a fission even occurs
        I then take the length of the list till a fission event occurs and use
        it to slice the original data'''
        subcell_list = list(split_before(self.daughters, lambda x: x != 0))
        subcell_list_lengths = [len(length_subcell) for length_subcell in subcell_list]
        all_subcells = []
        current_position = 0
        subcell_number_tracker = 1

        #Here I slice the original cell data and put it into subcell objects
        for number_of_frames_in_subcell in subcell_list_lengths:
            end_position = current_position+number_of_frames_in_subcell
            
            subcell_framenbs = np.asarray(self.framenbs[current_position:end_position])
            subcell_length = np.asarray(self.length[current_position:end_position])
            subcell_width = np.asarray(self.width[current_position:end_position])
            subcell_area = np.asarray(self.area[current_position:end_position])
            subcell_fluorescence = np.asarray(self.fluorescence[current_position:end_position])
            subcell_fluorescence_by_area = np.asarray(self.fluorescence_by_area[current_position:end_position])
            subcell_length_growth = np.asarray(self.length_growth[current_position:end_position-1])
            subcell_width_growth = np.asarray(self.width_growth[current_position:end_position-1])
            subcell_area_growth = np.asarray(self.area_growth[current_position:end_position-1])
            subcell_fluorescence_growth = np.asarray(self.fluorescence_growth[current_position:end_position-1])
            subcell_fluorescence_growth_by_area = np.asarray(self.fluorescence_growth_by_area[current_position:end_position-1])


            current_subcell = SubCell(
            subcell_number_tracker,
            subcell_framenbs,
            subcell_length,
            subcell_width,
            subcell_area,
            subcell_fluorescence,
            subcell_fluorescence_by_area,
            subcell_length_growth,
            subcell_width_growth,
            subcell_area_growth,
            subcell_fluorescence_growth,
            subcell_fluorescence_growth_by_area
            )

            current_position = end_position
            all_subcells.append(current_subcell)
            subcell_number_tracker += 1
        self.subcells= all_subcells
        
        return self.subcells

class SubCell:
    def __init__(self, 
                subcell_nb, 
                framenbs,
                length,
                width,
                area,
                fluorescence, 
                fluorescence_by_area,
                length_growth,
                width_growth,
                area_growth,
                fluorescence_growth,
                fluorescence_growth_by_area):
        self.subcell_nb = subcell_nb
        self.framenbs = framenbs
        self.length = length
        self.width = width
        self.area = area
        self.fluorescence = fluorescence
        self.fluorescence_by_area = fluorescence_by_area
        self.length_growth = length_growth
        self.width_growth = width_growth
        self.area_growth = area_growth
        self.fluorescence_growth = fluorescence_growth
        self.fluorescence_growth_by_area = fluorescence_growth_by_area
    
    def __repr__(self):
        return f"subcell: {self.subcell_nb}"

    def __len__(self):
        return len(self.framenbs)


