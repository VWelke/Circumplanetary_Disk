data = {
    'Sky coordinates': ['14h08m10s'] * 9,
    'Polarisation': ['Dual'] * 9,
    'Bandwidth per Polarization': ['7.5 GHz'] * 6 + ['15 GHz'] * 3,
    'Column Density': ['0.913mm (3rd Octile)'] * 9,
    'Resolution_C10 (arcsecond)': [0.1100, 0.0420, 0.0280, 0.0230, 0.0180, 0.0120, 0.0091, 0.0065, 0.0048],
    'Resolution_C9 (arcsecond)': [0.1400, 0.0570, 0.0380, 0.0310, 0.0250, 0.0170, 0.0120, 0.0088, 0.0066],
    'Resolution_C8 (arcsecond)': [0.240, 0.096, 0.064, 0.052, 0.042, 0.028, 0.021, 0.015, 0.011],
    'Resolution_C7 (arcsecond)': [0.530, 0.210, 0.140, 0.110, 0.092, 0.061, 0.046, 0.033, 0.024],
    'Band': [1, 3, 4, 5, 6, 7, 8, 9, 10],
    "Wavelength (micron)": [7494.81, 2997.92,1998.62, 1620.50,1304.45, 868.96, 651.72, 461.22, 344.59],
    'Wavelength (GHz)': [40, 100, 150, 185, 230, 345, 460, 650, 870],
    'beam_ALMA_C8 (arcsec^2)': [0.065266, 0.010443, 0.004641, 0.003064, 0.001999, 0.000888, 0.000500, 0.000255, 0.000137],
    'beam_ALMA_C7 (arcsec^2)': [0.318285, 0.049969, 0.022209, 0.013710, 0.009590, 0.004216, 0.002398, 0.001234, 0.000653],
    'Noise_C10 (Jy/pixel)': [7.703e-07, 3.689e-06, 6.300e-06, 9.125e-06, 9.573e-06, 1.926e-05, 2.999e-05, 3.040e-05, 2.864e-05],
    'Noise_C9 (Jy/pixel)': [1.248e-06, 6.795e-06, 1.160e-05, 1.658e-05, 1.847e-05, 3.866e-05, 5.214e-05, 5.571e-05, 5.414e-05],
    'Noise_C8 (Jy/pixel)': [3.667e-06, 1.927e-05, 3.291e-05, 4.664e-05, 5.212e-05, 1.049e-04, 1.597e-04, 1.619e-04, 1.504e-04],
    'Noise_C7 (Jy/pixel)': [1.788e-05, 9.223e-05, 1.575e-04, 2.087e-04, 2.501e-04, 4.978e-04, 7.662e-04, 7.835e-04, 7.159e-04],
    'Time_C10 (s)': [265880, 17000, 6800, 339900, 5140, 5100, 25320, 348350, 2274680],
    'Time_C9 (s)': [101300, 5010, 2010, 102960, 1380, 1270, 8380, 103730, 636550],
    'Time_C8 (s)': [11740, 630, 250, 13020, 180, 180, 900, 12290, 82490],
    'Time_C7 (s)': [500, 30, 20, 650, 10, 10, 40, 530, 3650]
}


def get_alma_data():
    return data


