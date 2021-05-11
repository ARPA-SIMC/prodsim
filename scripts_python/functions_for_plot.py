import argparse, numpy as np
from math import *

# -----------------------------------------------------------------------------------------
# FUNZIONI
# -----------------------------------------------------------------------------------------
# Input da linea di comando
def command_line():
    parser = argparse.ArgumentParser(description = 'PLOT GRIB CON PYTHON. Lo script può   \
        plottare: un campo scalare, un campo vettoriale (vento) o la sovrapposizione dei  \
        due. Ogni file GRIB dato in input deve contenere UN SOLO messaggio. Maggiori      \
        informazioni si possono trovare nello script.')

    parser.add_argument('--plot_scalar',default = False, action = 'store_true',
                        help = 'Attiva il plot di un campo scalare, il cui path è definito\
                                da --fname. Può essere usato in combinazione con          \
                                --plot_wind')
    parser.add_argument('--plot_wind',  default = False, action = 'store_true',
                        help = 'Attiva il plot del campo vettoriale di vento, definito dai\
                                campi U e V specificati da --fname_u e --fname_v. Può     \
                                essere usato in combinazione con --plot_scalar')
    parser.add_argument('--fname',      default = None,  type = str,
                        help = "Path del file da plottare (campo scalare)") 
    parser.add_argument('--fname_u',    default = None,  type = str,
                        help = "Path del file da plottare (componente U del campo         \
                                vettoriale)")
    parser.add_argument('--fname_v',    default = None,  type = str,
                        help = "Path del file da plottare (componente V del campo         \
                                vettoriale)")
    parser.add_argument('--fold_out',   default = "./",  type = str,
                        help = "Cartella in cui salvare il plot. Default: path corrente")
    parser.add_argument('--outname',    default = "",    type = str,
                        help = "Nome del file in output (se non specificato viene creato  \
                                automaticamente")
    parser.add_argument('--aver_wind',  default = False, action = 'store_true',
                        help = "Applica un filtro gaussiano al campo di vento 2D per      \
                                renderlo più smooth")
    parser.add_argument('--sigma',      default = 2,     type = float,
                        help = "Sigma della gaussiana per mediare il campo di vento (deve \
                                essere attivato --aver_wind). Default: 2")
    parser.add_argument('--scal_wind',  default = 80,    type = float,
                        help = "Fattore di scala che regola la lunghezza dei vettori del  \
                                vento.  Default: 80")
    parser.add_argument('--width_arrow',default = 0.002, type = float,
                        help = "spessore vettori vento. Default: 0.002")
    parser.add_argument('--subarea',    default = False, action = 'store_true',
                        help = "Attiva il ritaglio di una sottoarea del plot, secondo i   \
                                valori definiti da latmin, latmax, lonmin, lonmax")
    parser.add_argument('--latmin',     default = None,   type = float,
                        help = "Latitudine minima per l'opzione --subarea")
    parser.add_argument('--latmax',     default = None,   type = float,
                        help = "Latitudine massima per l'opzione --subarea")
    parser.add_argument('--lonmin',     default = None,   type = float,
                        help = "Longitudine minima per l'opzione --subarea")
    parser.add_argument('--lonmax',     default = None,   type = float,
                        help = "Longitudine massima per l'opzione --subarea")
    parser.add_argument('--grid',       default = False,  action = 'store_true',
                        help = "Attiva il plot della griglia con labels di latitudine e   \
                                 longitudine")
    parser.add_argument('--dpi',        default = 150,   type = int,
                        help = "Risoluzione immagine in dpi (default: 150)")
    parser.add_argument('--title',      default = "",     type = str,
                        help = "Titolo del plot. Se vuoto (default) viene generato        \
                                automaticamente. Se None, nessun titolo")
    return parser.parse_args()


# Ricavo alcune variabili dal grib (tutte le keys sono in grb.key())
def get_var_grib(grib_data):
    variabile = grib_data.shortName
    data      = grib_data.dataDate
    ora       = "%04d"  %grib_data.dataTime
    Ni        = int(grib_data['Ni'])
    Nj        = int(grib_data['Nj'])
    startstep = grib_data.startStep
    endstep   = "%02d" %grib_data.endStep
    hcum      = int(endstep) - int(startstep)
    
    print("\nVariabile: ", variabile,   "\nData:      ",   data,
          "\nOra:       ", ora,         "\nScadenza:  +", endstep, sep="")
    if hcum != 0: print("Cumulata su: ", hcum, "h", sep="")

    return variabile, data, ora, Ni, Nj, startstep, endstep, hcum


# Trasformazione di coordinate da "regular latlon" a "rotated latlon" (option=1) e 
# viceversa (option=0)
def rotated_grid_transform(lon_arr, lat_arr, option, SP_coor):
    # Conversione da gradi a radianti
    lon = (lon_arr*pi)/180
    lat = (lat_arr*pi)/180

    # Salvo coordinate del Southern Pole
    SP_lon, SP_lat = SP_coor[0], SP_coor[1];

    # Definisco l'angolo (in radianti) di rotazione attorno all'asse y (theta) e z (phi)
    theta = (90+SP_lat)*pi/180
    phi   = (SP_lon)*pi/180

    # Converto da coordinate sferiche a cartesiane
    x = np.multiply(np.cos(lon), np.cos(lat)) 
    y = np.multiply(np.sin(lon), np.cos(lat))
    z = np.sin(lat);

    if option == 1: # Regular -> Rotated 
        cosphi_x = np.multiply(np.cos(phi), x)
        senphi_y = np.multiply(np.sin(phi), y)
        x_new =  np.multiply(cosphi_x, np.cos(theta)) + \
                 np.multiply(senphi_y, np.cos(theta)) + np.multiply(np.sin(theta), z)
        y_new = -np.multiply(np.sin(phi), x) + np.multiply(np.cos(phi), y)
        z_new = -np.multiply(cosphi_x, np.sin(theta)) - \
                 np.multiply(senphi_y, np.sin(theta)) + np.multiply(np.cos(theta), z)

    else:  # Rotated -> Regular
        phi = -phi;
        theta = -theta;

        costheta_x = np.multiply(np.cos(theta), x)
        sentheta_z = np.multiply(np.sin(theta), z)

        x_new =  np.multiply(costheta_x, np.cos(phi)) + np.multiply(np.sin(phi), y) + \
                 np.multiply(sentheta_z, np.cos(phi))
        y_new = -np.multiply(costheta_x, np.sin(phi)) + np.multiply(np.cos(phi), y) - \
                 np.multiply(sentheta_z, np.sin(phi))
        z_new = -np.multiply(np.sin(theta), x)        + np.multiply(np.cos(theta), z)

    # Ri-converto da coordinate cartesiane a sferiche e trasformo in radianti
    lon_new = (np.arctan2(y_new,x_new))*180/pi
    lat_new = (np.arcsin(z_new))*180/pi

    return lon_new, lat_new
