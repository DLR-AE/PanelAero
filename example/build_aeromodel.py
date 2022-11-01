import numpy as np

class AeroModel():
    
    def __init__(self, filename):
        self.filename = filename

    def build_aerogrid(self):
        caero_grid, caero_panels = self.read_CAERO(self.filename, 0)
        ID = []
        l = [] # length of panel
        A = [] # area of one panel
        N = [] # unit normal vector 
        offset_l = [] # 25% point l
        offset_k = [] # 50% point k
        offset_j = [] # 75% downwash control point j
        offset_P1 = [] # Vortex point at 25% chord, 0% span
        offset_P3 = [] # Vortex point at 25% chord, 100% span
        r = [] # vector P1 to P3, span of panel
        
        
        for i_panel in range(len(caero_panels['ID'])):
    
            #
            #                   l_2
            #             4 o---------o 3
            #               |         |
            #  u -->    b_1 | l  k  j | b_2
            #               |         |
            #             1 o---------o 2
            #         y         l_1
            #         |
            #        z.--- x
    
            index_1 = np.where(caero_panels['cornerpoints'][i_panel][0]==caero_grid['ID'])[0][0]
            index_2 = np.where(caero_panels['cornerpoints'][i_panel][1]==caero_grid['ID'])[0][0]
            index_3 = np.where(caero_panels['cornerpoints'][i_panel][2]==caero_grid['ID'])[0][0]
            index_4 = np.where(caero_panels['cornerpoints'][i_panel][3]==caero_grid['ID'])[0][0]
            
            l_1 = caero_grid['offset'][index_2] - caero_grid['offset'][index_1]
            l_2 = caero_grid['offset'][index_3] - caero_grid['offset'][index_4]
            b_1 = caero_grid['offset'][index_4] - caero_grid['offset'][index_1]
            b_2 = caero_grid['offset'][index_3] - caero_grid['offset'][index_2]
            l_m = (l_1 + l_2) / 2.0
            b_m = (b_1 + b_2) / 2.0
            
            ID.append(caero_panels['ID'][i_panel])    
            l.append(l_m[0])
            # A.append(l_m[0]*b_m[1])
            A.append(np.linalg.norm(np.cross(l_m, b_m)))
            N.append(np.cross(l_1, b_1)/np.linalg.norm(np.cross(l_1, b_1)))
            offset_l.append(caero_grid['offset'][index_1] + 0.25*l_m + 0.50*b_1)
            offset_k.append(caero_grid['offset'][index_1] + 0.50*l_m + 0.50*b_1)
            offset_j.append(caero_grid['offset'][index_1] + 0.75*l_m + 0.50*b_1)
            offset_P1.append(caero_grid['offset'][index_1] + 0.25*l_1)
            offset_P3.append(caero_grid['offset'][index_4] + 0.25*l_2)
            r.append((caero_grid['offset'][index_4] + 0.25*l_2) - (caero_grid['offset'][index_1] + 0.25*l_1))
       
        n = len(ID)
        set_l = np.arange(n*6).reshape((n,6))
        set_k = np.arange(n*6).reshape((n,6))
        set_j = np.arange(n*6).reshape((n,6))
        aerogrid = {'ID': np.array(ID),
                    'l': np.array(l),
                    'A': np.array(A),
                    'N': np.array(N),
                    'offset_l': np.array(offset_l),
                    'offset_k': np.array(offset_k),
                    'offset_j': np.array(offset_j),
                    'offset_P1': np.array(offset_P1),
                    'offset_P3': np.array(offset_P3),
                    'r': np.array(r),
                    'set_l': set_l,
                    'set_k': set_k,
                    'set_j': set_j,
                    'CD': caero_panels['CD'],
                    'CP': caero_panels['CP'],
                    'n': n,
                    'coord_desc': 'bodyfixed',
                    'cornerpoint_panels': caero_panels['cornerpoints'],
                    'cornerpoint_grids': np.hstack((caero_grid['ID'][:,None],caero_grid['offset']))
                   }
        self.aerogrid = aerogrid
    
    def read_CAERO(self, filename, i_file):
        print('Read CAERO1 and/or CAERO7 cards from Nastran/ZAERO bdf: %s' %filename)
        caerocards = []
        with open(filename, 'r') as fid:
            while True:
                read_string = fid.readline()
                if str.find(read_string, 'CAERO1') !=-1 and read_string[0] != '$':
                    # read first line of CAERO card
                    caerocard = {'EID': nastran_number_converter(read_string[8:16], 'ID'),
                                 'CP': nastran_number_converter(read_string[24:32], 'ID'),
                                 'n_span': nastran_number_converter(read_string[32:40], 'ID'), # n_boxes
                                 'n_chord': nastran_number_converter(read_string[40:48], 'ID'), # n_boxes
                                 'l_span': nastran_number_converter(read_string[48:56], 'ID'),
                                 'l_chord': nastran_number_converter(read_string[56:64], 'ID'),
                                }
                    # read second line of CAERO card
                    read_string = fid.readline()  
                    caerocard['X1'] = np.array([nastran_number_converter(read_string[ 8:16], 'float'), nastran_number_converter(read_string[16:24], 'float'), nastran_number_converter(read_string[24:32], 'float')])
                    caerocard['length12'] = nastran_number_converter(read_string[32:40], 'float')
                    caerocard['X2'] = caerocard['X1'] + np.array([caerocard['length12'], 0.0, 0.0])
                    caerocard['X4'] =np.array([nastran_number_converter(read_string[40:48], 'float'), nastran_number_converter(read_string[48:56], 'float'), nastran_number_converter(read_string[56:64], 'float')])
                    caerocard['length43'] = nastran_number_converter(read_string[64:72], 'float')
                    caerocard['X3'] = caerocard['X4'] + np.array([caerocard['length43'], 0.0, 0.0])
                    caerocards.append(caerocard)
                if str.find(read_string, 'CAERO7') !=-1 and read_string[0] != '$':
                    # The CAERO7 cards of ZAERO is nearly identical to Nastran'S CAERO1 card. 
                    # However, it uses 3 lines, which makes the card more readable to the human eye.
                    # Also, not the number of boxes but the number of divisions is given (n_boxes = n_division-1)
                    # read first line of CAERO card
                    caerocard = {'EID': nastran_number_converter(read_string[8:16], 'ID'),
                                 'CP': nastran_number_converter(read_string[24:32], 'ID'),
                                 'n_span': nastran_number_converter(read_string[32:40], 'ID') - 1,
                                 'n_chord': nastran_number_converter(read_string[40:48], 'ID') - 1,
                                }
                    if np.any([caerocard['n_span'] == 0, caerocard['n_chord'] == 0]):
                        print('Assumption of equal spaced CAERO7 panels is violated!')
                    # read second line of CAERO card
                    read_string = fid.readline()  
                    caerocard['X1'] = np.array([nastran_number_converter(read_string[ 8:16], 'float'), nastran_number_converter(read_string[16:24], 'float'), nastran_number_converter(read_string[24:32], 'float')])
                    caerocard['length12'] = nastran_number_converter(read_string[32:40], 'float')
                    caerocard['X2'] = caerocard['X1'] + np.array([caerocard['length12'], 0.0, 0.0])
                    # read third line of CAERO card
                    read_string = fid.readline()  
                    caerocard['X4'] =np.array([nastran_number_converter(read_string[ 8:16], 'float'), nastran_number_converter(read_string[16:24], 'float'), nastran_number_converter(read_string[24:32], 'float')])
                    caerocard['length43'] = nastran_number_converter(read_string[32:40], 'float')
                    caerocard['X3'] = caerocard['X4'] + np.array([caerocard['length43'], 0.0, 0.0])
                    caerocards.append(caerocard)
                elif read_string == '':
                    break
        
        # from CAERO cards, construct corner points... '
        # then, combine four corner points to one panel
        grid_ID = i_file * 100000 # the file number is used to set a range of grid IDs 
        grids = {'ID':[], 'offset':[]}
        panels = {"ID": [], 'CP':[], 'CD':[], "cornerpoints": []}
        for caerocard in caerocards:
            # calculate LE, Root and Tip vectors [x,y,z]^T
            LE   = caerocard['X4'] - caerocard['X1']
            Root = caerocard['X2'] - caerocard['X1']
            Tip  = caerocard['X3'] - caerocard['X4']
            
            if caerocard['n_chord'] == 0:
                print('AEFACT cards are not supported by this reader.')
            else:
                # assume equidistant spacing
                d_chord = np.linspace(0.0, 1.0, caerocard['n_chord']+1 ) 
                 
            if caerocard['n_span'] == 0:
                print('AEFACT cards are not supported by this reader.')
            else:
                # assume equidistant spacing
                d_span = np.linspace(0.0, 1.0, caerocard['n_span']+1 ) 
            
            # build matrix of corner points
            # index based on n_divisions
            grids_map = np.zeros((caerocard['n_chord']+1,caerocard['n_span']+1), dtype='int')
            for i_strip in range(caerocard['n_span']+1):
                for i_row in range(caerocard['n_chord']+1):
                    offset = caerocard['X1'] \
                           + LE * d_span[i_strip] \
                           + (Root*(1.0-d_span[i_strip]) + Tip*d_span[i_strip]) * d_chord[i_row]
                    grids['ID'].append(grid_ID)
                    grids['offset'].append(offset)
                    grids_map[i_row,i_strip ] = grid_ID
                    grid_ID += 1
            # build panels from cornerpoints
            # index based on n_boxes
            panel_ID =  caerocard['EID']                  
            for i_strip in range(caerocard['n_span']):
                for i_row in range(caerocard['n_chord']):
                    panels['ID'].append(panel_ID)
                    panels['CP'].append(caerocard['CP']) # applying CP of CAERO card to all grids
                    panels['CD'].append(caerocard['CP'])
                    panels['cornerpoints'].append([ grids_map[i_row, i_strip], grids_map[i_row+1, i_strip], grids_map[i_row+1, i_strip+1], grids_map[i_row, i_strip+1] ])
                    panel_ID += 1 
        panels['ID'] = np.array(panels['ID'])
        panels['CP'] = np.array(panels['CP'])
        panels['CD'] = np.array(panels['CD'])
        panels['cornerpoints'] = np.array(panels['cornerpoints'])
        grids['ID'] = np.array(grids['ID'])
        grids['offset'] = np.array(grids['offset'])
        return grids, panels
    
def nastran_number_converter(string_in, type, default=0):
    if type in ['float']:
        try:
            out = float(string_in)
        except:
            
            string_in = string_in.replace(' ', '') # remove leading spaces
            for c in ['\n', '\r']: string_in = string_in.strip(c) # remove end of line
            if '-' in string_in[1:]:
                if string_in[0] in ['-', '+']:
                    sign = string_in[0]
                    out = float(sign + string_in[1:].replace('-', 'E-'))
                else:
                    out = float(string_in.replace('-', 'E-'))
            elif '+' in string_in[1:]:
                if string_in[0] in ['-', '+']:
                    sign = string_in[0]
                    out = float(sign + string_in[1:].replace('+', 'E+'))
                else:
                    out = float(string_in.replace('+', 'E+'))
            elif string_in == '':
                print("Could not interpret the following number: '" + string_in + "' -> setting value to "+str(default))
                out = float(default)
            else: 
                print("Could not interpret the following number: " + string_in)
                return
    elif type in ['int', 'ID', 'CD', 'CP']:
        try:
            out = int(string_in)
        except:
            out = int(default)  
    return out


