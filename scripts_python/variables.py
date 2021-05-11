# -----------------------------------------------------------------------------------------
# VARIABILI
# -----------------------------------------------------------------------------------------
# Ogni variabile è definita da una classe il cui nome corrisponde con il suo 'shortName'.
# Vanno specificati:
# varsn       = nome abbreviato (per nome file in output automatico)  
# varln       = nome esteso     (per titolo del plot automatico)
# units       = unità di misura (per titolo del plot automatico)
# scal_fact   = fattore per cui moltiplicare i dati
# scal_offset = offset da aggiungere ai dati
# thin_fact   = fattore di thinning dei dati in input (se 1, prende tutti i dati)
# livelli     = lista dei livelli per la colorbar. Se [] viene calcolato automaticamente
#               a partire dal minimo e dal massimo.
# ext_lim     = ['neither', 'both', 'min, 'max']. Determina la colorazione del contorno dei
#               valori che sono fuori dall'intervallo dei livelli (definisce la variabile
#               'extend' di contourf o analoghi). Se 'neither', i valori fuori 
#               dall'intervallo dei livelli non sono colorati. Se 'min', 'max' o 'both', 
#               colora i valori sotto, sopra o sotto e sopra l'intervallo dei livelli.
#               Se 'livelli = []' viene settato automaticamente a 'neiher'
# colori      = lista dei colori da usare nella colorbar. Se i livelli sono stati definiti,
#               cioè livelli != [], allora il numero di colori deve coincidere con quello
#               dei livelli - 1.

# Precipitazione
class tp:
    varsn       = 'tp'
    varln       = 'Accumulated precipitation'
    units       = 'mm'
    scal_fact   = 1
    scal_offset = 0.0
    thin_fact   = 1
    livelli     = [0.,1.,2.,5., 10., 15., 20., 40.]
    ext_lim     = 'neither'
    colori      = [[255/255,255/255,217/255], [237/255,248/255,177/255],
                   [ 90/255,195/255,190/255], [ 29/255,145/255,192/255],
                   [ 34/255, 94/255,168/255], [ 12/255, 44/255,132/255],
                   [  8/255, 29/255, 88/255]]

# Pressione media alla superficie
class pmsl:
    varsn       = 'pmsl'
    varln       = 'Mean sea level pressure'
    units       = 'hPa'
    scal_fact   = 0.01
    scal_offset = 0.0
    thin_fact   = 1
    livelli     = []
    ext_lim     = 'neither'
    colori      = [[ 31/255, 27/255,120/255], [ 49/255, 54/255,149/255],
                   [ 69/255,117/255,180/255], [116/255,173/255,209/255],
                   [171/255,217/255,233/255], [255/255,255/255,191/255],
                   [254/255,224/255,144/255], [253/255,174/255, 97/255],
                   [244/255,109/255, 67/255], [215/255, 48/255, 39/255],
                   [165/255,  0/255, 38/255]]

# Temperatura
class t:
    varsn       = 't'
    varln       = 'Temperature'
    units       = '°C'
    scal_fact   = 1
    scal_offset = -273.15
    thin_fact   = 1
    livelli     = []
    ext_lim     = 'neither'
    colori      =  [[ 26/255,152/255, 80/255], [145/255,207/255, 96/255],
                    [217/255,239/255,139/255], [255/255,255/255,191/255],
                    [254/255,224/255,139/255], [252/255,141/255, 89/255],
                    [215/255, 48/255, 39/255]]

# Temperatura a 2 metri
class t2m:
    varsn       = '2t'
    varln       = '2m Temperature'
    units       = '°C'
    scal_fact   = 1
    scal_offset = -273.15
    thin_fact   = 1
    livelli     = [5,12,15,18,21,24,27,40]
    ext_lim     = 'neither'
    colori      =  [[ 26/255,152/255, 80/255], [145/255,207/255, 96/255],
                    [217/255,239/255,139/255], [255/255,255/255,191/255],
                    [254/255,224/255,139/255], [252/255,141/255, 89/255],
                    [215/255, 48/255, 39/255]]

# Umidità del suolo
class w_so:
    varsn       = 'w_so'
    varln       = 'Soil moisture'
    units       = r'$kg/m^2$'
    scal_fact   = 1
    scal_offset = 0
    thin_fact   = 1
    livelli     = []
    livelli     = [ -40,  -20,  -10,   -2,   2,  10,  20,  40]      # lv 54
                 #[ -20,  -10,   -5,   -1,   1,   5,  10,  20]      # lv 18
                 #[-1.2, -0.8, -0.4, -0.1, 0.1, 0.4, 0.8, 1.2]      # lv  1
    ext_lim     = 'both'
    colori      = [[165/255,  0/255, 38/255], [215/255, 48/255, 39/255],
                   [244/255,109/255, 67/255], [253/255,174/255, 97/255],
                    'white',
                   [116/255,173/255,209/255], [ 69/255,117/255,180/255],
                   [ 49/255, 54/255,149/255], [ 24/255, 19/255,120/255]]

# Vento
# Essendo un campo vettoriale, non è necessario definire livelli, ext_lim e colori
class wind:
    varsn       = 'wind'
    varln       = 'Wind'
    units       = 'm/s'
    scal_fact   = 1
    scal_offset = 0.0
    thin_fact   = 10

# Tutte le altre variabili non definite dalla precedenti classi
class default:
    varsn       = 'plot'
    varln       = ''
    units       = ''
    scal_fact   = 1
    scal_offset = 0.0
    thin_fact   = 1
    livelli     = []
    ext_lim     = 'neither'
    colori      =  [[ 26/255,152/255, 80/255], [145/255,207/255, 96/255],
                    [217/255,239/255,139/255], [255/255,255/255,191/255],
                    [254/255,224/255,139/255], [252/255,141/255, 89/255],
                    [215/255, 48/255, 39/255]]
