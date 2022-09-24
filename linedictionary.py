# based on splatalogue and https://help.almascience.org/kb/articles/what-spectral-resolution-will-i-get-for-a-given-channel-spacing and original proposal for spectral resolution and channel width
# created by Y. Yamato in Sep. 4th, 2022

line_dict = {
    "C18O":     {"qn": "2--1",         "freq": 219.56035410, "logA": -6.22103, "E_u": 15.80580, "g_u": 5.0,  "spin": None, "dchan": 61.035,  "res": 122.070}, # CDMS
    "13CO":     {"qn": "2--1",         "freq": 220.39868420, "logA": -6.51752, "E_u": 15.86618, "g_u": 10.0, "spin": None, "dchan": 61.035,  "res": 122.070}, # CDMS
    "12CO":     {"qn": "2--1",         "freq": 230.53800000, "logA": -6.16050, "E_u": 16.59608, "g_u": 5.0,  "spin": None, "dchan": 244.141, "res": 488.282}, # CDMS
    "SO":       {"qn": "6$_5$--5$_4$", "freq": 219.94944200, "logA": -3.87446, "E_u": 34.98470, "g_u": 13.0, "spin": None, "dchan": 61.035,  "res": 122.070}, # CDMS
    "H2CO":     {"t1": {"qn": "3$_{2,1}$--2$_{2,0}$", "freq": 218.76006600, "logA": -3.80205, "E_u": 68.11081, "g_u": 7.0, "spin": "para", "dchan": 61.035,  "res": 122.070}, # CDMS
                 "t2": {"qn": "3$_{0,3}$--2$_{0,2}$", "freq": 218.22219200, "logA": -3.55007, "E_u": 20.95640, "g_u": 7.0, "spin": "para", "dchan": 488.281, "res": 976.562}, # CDMS
                 "t3": {"qn": "3$_{2,2}$--2$_{2,1}$", "freq": 218.47563200, "logA": -3.80373, "E_u": 68.09370, "g_u": 7.0, "spin": "para", "dchan": 488.281, "res": 976.562}, # CDMS
                },
    "c-C3H2":   {"t1": {"qn": "6$_{1,6}$--5$_{0,5}$", "freq": 217.82215, "logA": -3.26791, "E_u": 38.60744, "g_u": 39.0, "spin": "ortho", "dchan": 488.281, "res": 976.562}, # JPL
                 "t2": {"qn": "6$_{0,6}$--5$_{1,5}$", "freq": 217.82215, "logA": -3.26789, "E_u": 38.60744, "g_u": 13.0, "spin": "para",  "dchan": 488.281, "res": 976.562}, # JPL
                 "t3": {"qn": "5$_{1,4}$--4$_{2,3}$", "freq": 217.94005, "logA": -3.35391, "E_u": 35.41701, "g_u": 33.0, "spin": "ortho", "dchan": 488.281, "res": 976.562}, # JPL
                 "t4": {"qn": "5$_{2,4}$--4$_{1,3}$", "freq": 218.16044, "logA": -3.35225, "E_u": 35.41781, "g_u": 11.0, "spin": "ortho", "dchan": 488.281, "res": 976.562}, # JPL
                },
    "DCN":      {"qn": "3--2",         "freq": 217.23853780, "logA": -3.33964, "E_u": 20.85164, "g_u": 21.0, "spin": None, "dchan": 488.281, "res": 976.562}, # CDMS
    "CH3OH":    {"qn": "4$_2$--3$_1$", "freq": 218.44006300, "logA": -4.32917, "E_u": 45.45988, "g_u": 9.0,  "spin": None, "dchan": 488.281, "res": 976.562}, # JPL
    "SiO":      {"qn": "5--4",         "freq": 217.10498000, "logA": -3.28288, "E_u": 31.25889, "g_u": 11.0, "spin": None, "dchan": 488.281, "res": 976.562}, # JPL
    }

# transition list following nomenclature in imaging script by John
transition_list = [
    "C18O", # main line
    "13CO", # main line
    "12CO",
    "SO", # main line
    "H2CO_3_21-2_20_218.76GHz",
    "H2CO_3_03-2_02_218.22GHz",
    "H2CO_3_22-2_21_218.47GHz",
    "c-C3H2_217.82",
    "cC3H2_217.94",
    "cC3H2_218.16",
    "DCN",
    "CH3OH",
    "SiO",
]