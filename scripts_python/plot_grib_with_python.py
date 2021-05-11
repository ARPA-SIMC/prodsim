import os, sys, pygrib, numpy as np
import cartopy, cartopy.crs as ccrs, cartopy.feature as cfeature
import matplotlib, matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors
from variables import *
from functions_for_plot import *

# -----------------------------------------------------------------------------------------
# PLOT GRIB CON PYTHON - AGGIORNATO AL 03/06/2020
# -----------------------------------------------------------------------------------------
# DESCRIZIONE
# Lo script può plottare:
#   - un campo scalare
#   - un campo vettoriale (vento)
#   - un campo vettoriale (vento) sovrapposto ad un campo scalare.
# Ogni file GRIB dato in input deve contenere UN SOLO messaggio.
#
# I parametri che definiscono cosa plottare (path del file, area, ecc.) possono essere 
# definiti da linea di comando oppure nella sezione INPUT all'interno di questo script
# (vedi DEFINIZIONE INPUT per maggiori dettagli).
#
# Le specifiche relative alla variabile plottata (scala colori, livelli...) sono definite
# mediante una una classe (identificata dallo shortName) nel file esterno 'variables.py'.
# Se la variabile che vuoi plottare non è presente devi:
#   - aggiungere la classe in 'variables.py';
#   - associare la classe alla variabile modificando le righe di codice sotto il commento
#     "Associo la variabile letta nel grib alla classe che la definisce"
# Altrimenti ti accontenti del default.
#
#
# DEFINIZIONE INPUT
# La descrizione delle variabili si ottiene:
#   - per input da linea di comando       (LC): digitando --help
#   - per input dall'interno dello script (IS): leggendo la sezione INPUT
#
# Va sempre definito il tipo di campo da plottare. Quindi:
#   - per un campo scalare:     --plot_scalar (LC) o plot_scalar=True (IS)
#   - per un campo vettoriale:  --plot_wind   (LC) o plot_wind=True   (IS)
#   - per entrambi:             --plot_scalar      e --plot_wind      (LC) o 
#                               plot_scalar=True   e plot_wind=True   (IS)
#
#
# PRIMO UTILIZZO
# 1) Assicurati di aver installato i pacchetti "pygrib" e "cartopy" per python3.
# 2) Vai su "trentadue" (o analogo): serve per lo scaricamento dei files relativi alla
#    geografia. Dal secondo utilizzo potrai lanciarlo dal tuo PC.
# 3) Lancialo: python3 plot_grib_with_python.py [--OPTIONAL_ARGUMENTS]
#
# Nota che il punto 2) va effettuato anche se decidi di aggiungere dei files relativi
# alla geografia o se decidi di usare una differente risoluzione o origine.


# -----------------------------------------------------------------------------------------
# INPUT
# -----------------------------------------------------------------------------------------
# SPECIFICHE I/O
# plot_scalar = se True, plotta un campo scalare; va devinito 'fname'
# plot_wind   = se True, plotta un campo vettoriale (vento): vanno definiti 'fname_u' e 
#               'fname_v'.
# fname       = file da plottare; necessario se "plot_scalar = True"
# fname_u     = file da plottare contenente la componente U del vento; necessario se 
#               "plot_wind = True"
# fname_v     = file da plottare contenente la componente V del vento; necessario se
#               "plot_wind = True"
# fold_out    = cartella in cui verranno salvati i plot
# outname     = nome assegnato al file salvato (senza .png). Se non specificato ("") viene
#               creato automaticamente

plot_wind     = False
plot_scalar   = True
fname         = ""
fname_u       = ""
fname_v       = ""
fold_out      = "./"
outname       = ""


# SPECIFICHE PER IL VENTO
# aver_wind   = se True, applica un filtro gaussiano al campo di vento 2D per renderlo
#               più smooth
# sigma       = sigma della gaussiana; necessario se "aver_wind = True"
# scal_wind   = fattore di scala che regola la lunghezza dei vettori del vento
# width_arrow = spessore frecce

aver_wind     = False
sigma         = 2 
scal_wind     = 80
width_arrow   = 0.002


# SPECIFICHE PLOT
# subarea     = se True, ritaglia il plot, secondo i valori autoesplicativi definiti da
#               latmin, latmax, lonmin, lonmax.
# grid        = se True, plotta una griglia con labels di latitudine e longitudine
# title       = se ""   : titolo generato automaticamente
#               se "str": la stringa 'str' sarà il titolo del plot
#               se None : nessun titolo
# dpi         = risoluzione in dpi

subarea       = False
latmin        = 43
latmax        = 46
lonmin        = 11
lonmax        = 14
grid          = False
title         = ""
dpi           = 150


# -----------------------------------------------------------------------------------------
# LETTURA PARAMETRI DA LINEA DI COMANDO E CHECK SULLE VARIABILI DEFINITE
# -----------------------------------------------------------------------------------------
# Se lo script è stato lanciato con almeno un argomento, sovrascrivo i valori delle 
# variabili definiti ad inizio programma nella sezione "INPUT". Per --title, visto che
# è un str, lo pongo None se è 'None' o 'none'
if len(sys.argv) > 1:
    args = command_line()
    plot_scalar, plot_wind  = args.plot_scalar, args.plot_wind
    fname, fname_u, fname_v = args.fname,       args.fname_u, args.fname_v
    fold_out,  outname      = args.fold_out,    args.outname
    aver_wind, sigma        = args.aver_wind,   args.sigma
    scal_wind, width_arrow  = args.scal_wind,   args.width_arrow
    subarea                 = args.subarea
    latmin, latmax          = args.latmin,      args.latmax
    lonmin, lonmax          = args.lonmin,      args.lonmax
    grid                    = args.grid
    dpi                     = args.dpi
    title                   = args.title
    if title == 'None' or title == 'none': title = None 

# Check: almeno uno tra 'plot_scalar' e 'plot_wind' deve essere attivo (True)
if not (plot_scalar or plot_wind):
    sys.exit("\nERRORE! Devi selezionare almeno uno tra 'plot_scalar' e 'plot_wind'\n")

# Check: se 'plot_scalar', allora 'fname' deve essere definito
if plot_scalar and fname is None:
     sys.exit("\nERRORE! Hai selezionato 'plot_scalar' quindi devi definire 'fname'\n")

# Check: se 'plot_scalar', allora 'fname_u' e 'fname_v' devono essere definiti
if plot_wind  and (fname_u is None or fname_v is None):
     sys.exit("\nERRORE! Hai selezionato 'plot_wind' quindi devi definire sia         \
               \n'fname_u' che 'fname_v'\n")

# Check: se è attivata 'subarea', devono esserci i valori di 'latmin', 'latmax', 'lonmin',
# 'lonmax' e devono essere sensati
if subarea:
    if not (isinstance(latmin, (int, float)) and isinstance(latmax, (int, float)) and \
            isinstance(lonmin, (int, float)) and isinstance(lonmax, (int, float))):
        sys.exit("\nERRORE! Se attivi 'subarea' devi anche specificare correttamente  \
                  \nlatmin, latmax, lonmin, lonmax\n")
    if latmax <= latmin:
        sys.exit("\nERRORE! 'latmin' deve essere minore di 'latmax'\n")
    if lonmax <= lonmin:
        sys.exit("\nERRORE! 'lonmin' deve essere minore di 'lonmax'\n")


# -----------------------------------------------------------------------------------------
# LETTURA FILE GRIB
# -----------------------------------------------------------------------------------------
# Definisco cosa plottare (serve per sovrappore il campo di vento ad un campo scalare)
if plot_wind and plot_scalar:
    plot_seq = ['scalar', 'wind']
elif plot_wind and (not plot_scalar):
    plot_seq = ['wind']
elif (not plot_wind) and plot_scalar:
    plot_seq = ['scalar']
else:
    sys.exit("ERRORE! Devi decidere cosa plottare!")

# Loop sul tipo di file da plottare
for pseq in plot_seq:
    # Apro il file grib e seleziono il primo messaggio (al momento il programma non è
    # pensato per leggerne più di uno)
    print("\n\n---------------------------------------------------------------------")
    if pseq == 'wind':
        print("Lettura file:", fname_u)
        print("             ", fname_v)
        grbs_u = pygrib.open(fname_u)
        grbs_v = pygrib.open(fname_v)
        grb_u  = grbs_u.select()[0]
        grb_v  = grbs_v.select()[0]
    else:
        print("Lettura file:", fname)
        grbs   = pygrib.open(fname)
        grb    = grbs.select()[0]

    # Ricavo alcune variabili del grib. Per il vento controllo che i due campi siano 
    # consistenti (si dovrebbe controllare che corrispandano anche ad una U e V).
    if pseq == 'wind':
        grib_info_u = get_var_grib(grb_u)
        grib_info_v = get_var_grib(grb_v)
        if grib_info_u[1:] != grib_info_v[1:]:
            print("\nERRORE! Almeno uno dei parametri letti nei due file",
                  "\n(data, ora, Ni, Nj, startstep, endstep, hcum) non coincide!",
                  "\nfile contenete U: ", grib_info_u, 
                  "\nfile contenete V: ", grib_info_v)
            sys.exit()
        variabile, data, ora, Ni, Nj, startstep, endstep, hcum = grib_info_u
    else:
        variabile, data, ora, Ni, Nj, startstep, endstep, hcum = get_var_grib(grb)

    # Associo la variabile letta nel grib alla classe che la definisce (il vento è
    # gestito separatamente).
    # Devi associare la variabile ad una nuova classe? Aggiungi:
    # elif:
    #   var = nuovo_shortName
    if   variabile == 'pmsl':
        var = pmsl
    elif variabile == 'tp':
        var = tp
    elif variabile == 't':
        var = t
    elif variabile == '2t':
        var = t2m
    elif variabile == 'w_so':
        var = w_so
    else:
        var = default
    if pseq == 'wind': var = wind

    # Se sovrappongo il vento al campo scalare, salvo la variabile del campo scalare per
    # tenerne traccia per il titolo e il nome del file
    if plot_wind and  pseq == 'scalar': var_scalar = var

    # Salvo i valori della variabile letta e li modifico in base all'offset e al fattore di
    # scala definiti nella classe. Per il vento, pongo "grb = grb_u" per comodità nel
    # calcolo della griglia
    if pseq == 'wind':
        if aver_wind:
            dati_u = grb_u.values
            dati_v = grb_v.values
            sig    = [sigma, sigma]
            tmp_u  = sp.ndimage.filters.gaussian_filter(dati_u, sig, mode='constant')
            tmp_v  = sp.ndimage.filters.gaussian_filter(dati_v, sig, mode='constant')
            dati_u = tmp_u[0::var.thin_fact, 0::var.thin_fact]
            dati_v = tmp_v[0::var.thin_fact, 0::var.thin_fact]
        else:
            dati_u = grb_u.values[0::var.thin_fact, 0::var.thin_fact]
            dati_v = grb_v.values[0::var.thin_fact, 0::var.thin_fact]
        dati_u = dati_u*var.scal_fact + var.scal_offset
        dati_v = dati_v*var.scal_fact + var.scal_offset
        grb    = grb_u
    else:
        dati   = grb.values[0::var.thin_fact, 0::var.thin_fact]
        dati   = dati*var.scal_fact   + var.scal_offset


    # -----------------------------------------------------------------------------------------
    # CREAZIONE GRIGLIA E RITAGLIO
    # -----------------------------------------------------------------------------------------
    # Se la griglia non è "rotated_ll" o "regular_ll", non so gestirla quindi esco 
    if not (grb.gridType == "rotated_ll"  or  grb.gridType == "regular_ll"):
        sys.exit("\nERRORE! La griglia %s non è gestita\n" %grb.gridType)

    # Ricavo latitudine e longirudine del primo e dell'ultimo punto della griglia e controllo
    # la consistenza dei valori di longitudine
    lat_first_point = float(grb['latitudeOfFirstGridPointInDegrees'])
    lon_first_point = float(grb['longitudeOfFirstGridPointInDegrees'])
    lat_last_point  = float(grb['latitudeOfLastGridPointInDegrees'])
    lon_last_point  = float(grb['longitudeOfLastGridPointInDegrees'])
    if lon_first_point > lon_last_point: lon_first_point = lon_first_point - 360

    # Ricavo longitudini e latitudini
    lons = np.linspace(lon_first_point, lon_last_point, Ni)[0::var.thin_fact]
    lats = np.linspace(lat_first_point, lat_last_point, Nj)[0::var.thin_fact]
    grid_lon, grid_lat = np.meshgrid(lons, lats)

    # Se la griglia è "rotated_ll", la ruoto
    if grb.gridType == "rotated_ll":
        SP = [float(grb['longitudeOfSouthernPoleInDegrees']),
              float(grb['latitudeOfSouthernPoleInDegrees'])]
        grid_lon, grid_lat = rotated_grid_transform(grid_lon, grid_lat, 0, SP)

    # Controllo che la griglia sia a posto
    if pseq == 'wind' and ((np.shape(dati_u) != np.shape(grid_lon)) or \
                           (np.shape(dati_u) != np.shape(grid_lat))):
        sys.exit("ERRORE! La griglia spaziale non coincide con quella di U e V!")
    elif (pseq == 'scalar') and ((np.shape(dati) != np.shape(grid_lon)) or \
                                 (np.shape(dati) != np.shape(grid_lat))):
        sys.exit("ERRORE! La griglia spaziale non coincide con quella della variabile!")
    
    
    # -----------------------------------------------------------------------------------------
    # PLOT
    # -----------------------------------------------------------------------------------------
    # Set-up della mappa
    if not (pseq == 'wind' and len(plot_seq) == 2):
        fig = plt.figure(figsize=(11,11), dpi=dpi)
        ax = plt.axes(projection=ccrs.PlateCarree())
        font = {'size': 14}
        plt.rc('font', **font)

    # Plot dati vento/campo
    if pseq == 'wind':
        pcm = ax.quiver(grid_lon, grid_lat, dati_u, dati_v, scale=scal_wind,
                        scale_units='inches', width=width_arrow)
                        # headwidth=3, headlength=2, headaxislength=2, width=0.0012)
        # Legenda
        qk  = ax.quiverkey(pcm, 0.9, 1.02, 10, '10 m/s', labelpos='E', coordinates='axes')
    else:
        # Definisco la colormap
        if var.livelli == []:
            levmin, levmax  = np.amin(dati), np.amax(dati)
            levstep = (levmax - levmin)/len(var.colori)
            livelli = np.arange(levmin, 1.00001*levmax, levstep)
            ext     = 'neither'
        else:
            livelli = var.livelli
            ext     = var.ext_lim
        cmap, norm = from_levels_and_colors(livelli, var.colori, extend=ext)

        # Plot dati
        pcm     = ax.pcolormesh(grid_lon, grid_lat, dati, cmap=cmap, norm=norm)

        # Colorbar
        cax, kw = matplotlib.colorbar.make_axes(ax,location='right',pad=0.05,shrink=0.7)
        out     = fig.colorbar(pcm,cax=cax,extend=ext, ticks=livelli, **kw)

# Aggiungo la linea di costa, i confini e i laghi. Note:
# - Si possono usare i file di GSHHS o NaturalEarth (qui si è usato il secondo). In 
#   ogni caso, per scaricarli la PRIMA volta, lo script va lanciato su 'trenta'; poi 
#   rimngono salvati 
# - Possibili risoluzioni con NaturalEarth: 10m, 50m, 110m
# - Esempio con GSHHS:
#     ax.add_feature(cfeature.GSHHSFeature(scale='i', levels=[1,2,3]));
ax.coastlines(resolution='10m') 
ax.add_feature(cfeature.NaturalEarthFeature(category='cultural', 
               name='admin_0_boundary_lines_land', scale='10m'),
               edgecolor='k', facecolor='none')
ax.add_feature(cartopy.feature.LAKES.with_scale('10m'), edgecolor='k', facecolor='none')

# Ritaglio sottoarea
if subarea: 
    ax.set_extent([lonmin, lonmax, latmin, latmax])

# Plot griglia con labels di latitudine e longitudine a sinistra e in basso
if grid:
    gl = ax.gridlines(draw_labels=True, linestyle='--')
    gl.top_labels   = False
    gl.right_labels = False

# Titolo
if title is not None:
    # Titolo creato automaticmente
    if title == "":
        # Nome variabile e unità di misura
        if plot_scalar and plot_wind:
            var_tit = "%s [%s] and %s [%s]" %(var_scalar.varln, var_scalar.units, 
                      var.varln, var.units)
        else:
            if var.varln != '' and var.units != '':
                var_tit = "%s [%s]" %(var.varln, var.units)
            else:
                var_tit = "%s [%s]" %(grb.name, grb.units)

        # In base alla lunghezza di 'var_tit' decido se separare le due parti del titolo
        # con un - o se spezzarlo in due righe
        if len(var_tit) < 30:
            sep = '-'
        else:
            sep = '\n'

        # Separazione ore e minuti
        ora_str = "%s:%s" %(str(ora)[:2], str(ora)[2:4])

        # Completamento con data e ora
        if   hcum != 0 and int(endstep) != 0:
            title = "%sh %s %s %s %s at +%sh" %(hcum, var_tit, sep, data, ora_str, endstep)
        elif hcum != 0 and int(endstep) == 0:
            title = "%sh %s %s %s %s at +%sh" %(hcum, var_tit, sep, data, ora_str, hcum)
        elif hcum == 0 and int(endstep) != 0:
            title = "%s %s %s %s at +%sh" %(var_tit, sep, data, ora_str, endstep)
        else:
            title = "%s %s %s %s" %(var_tit, sep, data, ora_str)
        
    # Titolo manuale
    else:
        tit = title

    # Creazione titolo
    ax.set_title(title)


# -----------------------------------------------------------------------------------------
# SALVATAGGIO
# -----------------------------------------------------------------------------------------
# Creo la directory "fold_out" se non esiste
if not os.path.exists("%s" %fold_out):
    os.makedirs("%s" %fold_out)

# Nome del file automatico
if outname == "":
    # Identificativo per il nome della variabile
    if plot_scalar and  plot_wind:
        varname = "%s_wind" %var_scalar.varsn
    else:
        varname =  var.varsn
    
    # Completamento con informazioni su data, ora ecc.
    if   hcum != 0 and int(endstep) != 0:
        fileout = "%s/%s%sh_%s%s_+%sh.png" %(fold_out, varname, hcum, data, ora, endstep)
    elif hcum != 0 and int(endstep) == 0:
        fileout = "%s/%s%sh_%s%s_+%sh.png" %(fold_out, varname, hcum, data, ora, hcum)
    elif hcum == 0 and int(endstep) != 0:
        fileout = "%s/%s_%s%s_+%sh.png" %(fold_out, varname, data, ora, endstep)
    else:
        fileout = "%s/%s_%s%s.png" %(fold_out, varname, data, ora)

# Nome del file manuale
else:
    fileout = "%s/%s.png" %(fold_out, outname)

# Salvataggio
fig.savefig(fileout, bbox_inches='tight')

