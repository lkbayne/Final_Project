import numpy as np

def calc_Si_conc_and_uncert(avg_abs, Si_yint, Si_slope, Sxo_Si, Si_slope_err, Si_yint_err):
    '''
    Converts measured average absorbance to Si concentration (μmol/L) and calculates combined uncertainty.

    Parameters
    ----------
    avg_abs : float
        Average absorbance measured from sample
    Si_yint : float
        Y-intercept of Si calibration curve
    Si_slope : float
        Slope of Si calibration curve
    Sxo_Si : float
        Standard deviation of absorbance blank (instrument scatter)
    Si_slope_err : float
        Standard error of the slope from calibration
    Si_yint_err : float
        Standard error of the intercept from calibration

    Returns
    -------
    conc : float
        Calculated Si concentration (μmol/L)
    u_combined : float
        Combined uncertainty (μmol/L)
    '''
    conc = (avg_abs - Si_yint) / Si_slope
    
    #Calculate uncertainty 
    u_abs_term = (1 / Si_slope) * Sxo_Si
    u_slope_term = ((avg_abs - Si_yint) / (Si_slope**2)) * Si_slope_err
    u_yint_term = (1 / Si_slope) * Si_yint_err
    
    #Combine uncertainties
    u_combined = np.sqrt(u_abs_term**2 + u_slope_term**2 + u_yint_term**2)
    
    return conc, u_combined


def calc_PO4_conc_and_uncert(avg_abs, a, b, c, Sxo_PO4):
    '''
    Converts measured average absorbance to PO4 concentration (μmol/L) using a quadratic calibration, and calculates combined uncertainty.

    Parameters
    ----------
    avg_abs : float
        Average absorbance measured from sample
    a, b, c : floats
        Quadratic calibration coefficients
    Sxo_PO4 : float
        Standard deviation of absorbance blank (instrument scatter)

    Returns
    -------
    conc : float
        Calculated PO4 concentration (μmol/L)
    u_combined : float
        Combined uncertainty (μmol/L)
    '''
    #Calculate concentration
    conc = poly_abs_to_conc(a, b, c, avg_abs)
    
    # Absorbance by ±Sxo_PO4
    conc_plus = poly_abs_to_conc(a, b, c, avg_abs + Sxo_PO4)
    conc_minus = poly_abs_to_conc(a, b, c, avg_abs - Sxo_PO4)
    
    # Calculate sample uncertainty from half the difference
    u_sample = abs(conc_plus - conc_minus) / 2
    
    # Combine sample and blank uncertainties
    u_combined = np.sqrt(Sxo_PO4**2 + u_sample**2)
    
    return conc, u_combined


    
def conc_to_dil(c1=None, v1=None, c2=None, v2=None):
    '''
    stock calculations to get from concentrated to diluted concentrations and volumes

    c1 * v1 = c2 * v2

    Parameters
    ----------
    c: concentration 
    c1 and c2 in the same units
    v: volume
    v1 and v2 in the same units

    Return
    ------
    concentration or volume not listed

    example and helpful print statement:
    c1 = 4000.65
    c2 = 5
    v1 = None
    v2 = 50
    calc_v1 = cc.conc_to_dil(c1, v1, c2, v2)
    print(f'Spike volume= {calc_v1:.2f} mL')
    '''
    if not c1:
        output = (c2 * v2) / v1
    if not c2:
        output = (c1 * v1) / v2
    if not v1:
        output = (c2 * v2) / c1
    if not v2:
        output = (c1 * v1) / c2
    return (output)


import math
def poly_abs_to_conc(a, b, c, ABS):
    '''
    absorbance measured to concentration with polynomial fit using the 
    quadratic equation in uM

    Parameters
    ----------
    ABS: absorbance measured
    a: a from polynomial calibration curve
    b: b from polynomial calibration curve
    c: c from polynomial calibration curve

    '''
    conc = (-b + math.sqrt(b**2 - 4*a*(c - ABS))) / (2*a)
    return(conc)

