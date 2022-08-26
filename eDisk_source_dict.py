# ver0; based on SourceList_noCam_sai.pdf
source_dict = {
    "L1489IRS": {
        # "radec": "04h04m43.08s 26d18m56.10s",
        "radec": '4h04m43.07997408s 26d18m56.11868681s', # updated 2022.06.08 based on 2D gaussian fit on the image plane
        "rep_robust": 1.0,
        "distance": 146, # Roccatagliata et al. 2020
        "PA": 67.2, # from visibility fit
        "incl": 70.6, # from visibility fit
        "vsys": 7.22, # LSR systemic velocity in km/s; Sai et al. 2020
        "emission_extent": {"12CO": (-13.0, 27), # updated 2022.06.19 based on inspection on the channel map
                            "13CO": (-2.0, 16), # updated 2022.06.19 based on inspection on the channel map
                            "C18O": (-1.0, 15), # updated 2022.06.19 based on inspection on the channel map
                            "SO": (-8, 22)} # inspected on casaviewer by eye in robust=0.5 images
    },
    "IRAS04169": {
        "radec": "04h19m58.449s  27d09m56.936s", # Takakuwa et al. 2018 (2D Gaussian fit)
    },
    "IRAS04302": {
        "radec": "04h33m16.49977s +22d53m20.225224s",
    },
    "Ced110IRS4": {
        "radec": "11h06m46.37687s -77d22m32.881218s",
    },
    "GSS30IRS3": {
        "radec": "16h26m21.72s -24d22m50.7s",
    },
    "OphIRS43": {
        "radec": "16h27m26.905457s -24d40m50.83194s",
    },
    "OphIRS63": {
        "radec": "16h31m35.70s -24d01m29.6s",
    },
    "IRS5N": {
        "radec": "19h01m48.479616s -36d57m15.38531s",
    },
    "IRS7B": {
        "radec": "19h01m56.419063s -36d57m28.67292s",
    },
    "IRAS32": {
        "radec": "19h02m58.72279s -37d07m37.38115s",
    },
    "IRAS04166+2706": {
        "radec": "04h19m42.50s  27d13m36.0s",
    },
    "L1527IRS": {
        "radec": "04h39m53.91s  26d03m09.8s",
    },
    "BHR71_IRS1": {
        "radec": "12h01m36.474422s -65d08m49.35978s", 
    },
    "BHR71_IRS2": {
        "radec": "12h01m34.01015s -65d08m48.0695s",
    },
    "IRAS15398": {
        "radec": "15h43m02.23327s -034d09m06.943163s",
    },
    "IRAS16253": {
        "radec": "16h28m21.616s -24d36m24.33s",
        "PA": 113,
        "incl": 70.6,
        "rep_robust": 0.5
    },
    "CB68": {
        "radec": "16h57m19.642s -16d09m24.03s",
        "PA": 45, #?
        "incl": 74,
        "rep_robust": 1.0
    },
}
