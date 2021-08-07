#!/usr/bin/env python
"""
Implementation of the Daisyworld model described in:

    Watson, A.J.; Lovelock, J.E (1983). "Biological homeostasis of
    the global environment: the parable of Daisyworld". Tellus.35B:
    286–9.

    Bibcode:1983TellB..35..284W. doi:10.1111/j.1600-0889.1983.tb00031.x.

Copyright (c) 2021 Fausto Zanasi
All rights reserved.

Redistribution and use in source and binary forms are permitted
provided that the above copyright notice and this paragraph are
duplicated in all such forms and that any documentation,
advertising materials, and other materials related to such
distribution and use acknowledge that the software was developed
by the authors. The name of the authors may not be used to endorse
or promote products derived from this software without specific prior
written permission.

THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
"""

# Imports
import numpy as np
import matplotlib.pyplot as plt

REVERSE = False
# Temperature
KELVIN_OFFSET = 273.15
temperaturaMargheriteMinima = 5 + KELVIN_OFFSET
temperaturaMargheriteMassima = 40 + KELVIN_OFFSET
temperaturaIdealeMargheriteNere = 22.5 + KELVIN_OFFSET
temperaturaIdealeMargheriteBianche = 22.5 + KELVIN_OFFSET

# Albedo
albedoBianche = 0.75
areaBianche = 0.01 # 0.01
albedoNere = 0.25
areaNere = 0.01 # 0.01
areaTotale = 1
albedoVuoto = 0.5
isolamento = 1 # 20 # isolamento fattore di proporzionalita'  tra temperatura e albedo

tassoDecrescita = 0.3 # le margherite hanno una vita media di 3.3 anni. 
# Calano di un 30% (0.3) all'anno

# Criteri di convergenza
numeroMassimoIterazioni = 1000
tolleranza = 0.000001

# livelloEnergiaSolare terms
costanteIrraggiamento = 917 # (W/mq)
sigma = 5.67032e-8 # costante di Stefan-Boltzmann

# EnergiaSolare (fattore di luminosità) limiti e step
energiaSolare_min = 0.5
energiaSolare_max = 1.3 # 1.6
energiaSolare_step = 0.002

# If run from command line, do the whole thing
if __name__ == '__main__':
    """Attiva la simulazione di daisyworld"""
    livelliEnergiaSolare = np.arange(energiaSolare_min, energiaSolare_max, energiaSolare_step)
    if REVERSE:
        livelliEnergiaSolare = livelliEnergiaSolare[::-1]

    areaNereStoria = np.zeros_like(livelliEnergiaSolare)
    areaBiancheStoria = np.zeros_like(livelliEnergiaSolare)
    areaVuotoStoria = np.zeros_like(livelliEnergiaSolare)
    temperaturaPianetaStoria = np.zeros_like(livelliEnergiaSolare)
    temperaturaSenzaVita = np.zeros_like(livelliEnergiaSolare)

    # Loop over livelliEnergiaSolare
    for j, livelloEnergiaSolare in enumerate(livelliEnergiaSolare):
        # Minimum daisy coverage
        if areaNere < 0.01:
            areaNere = 0.01
        if areaBianche < 0.01:
            areaBianche = 0.01
        areaVuoto = 1 - (areaNere + areaBianche )

        # Ripristina valori iniziali metrica iterazione
        it = 0
        cambiamentoPrecedenteIterazioneNere = 2*tolleranza
        cambiamentoPrecedenteIterazioneBianche = 2*tolleranza
        cambiamentoAreaNerePrecedente = 0
        cambiamentoAreaBianchePrecedente = 0
        temperaturaMediaPianeta = 0

        while it <= numeroMassimoIterazioni and cambiamentoPrecedenteIterazioneNere > tolleranza and cambiamentoPrecedenteIterazioneBianche > tolleranza:
            # albedo del pianeta
            albedoPianeta = (areaNere * albedoNere
                     + areaBianche * albedoBianche
                     + areaVuoto * albedoVuoto) / areaTotale
            # temperatura media del pianeta
            temperaturaMediaPianeta = np.power(livelloEnergiaSolare*(1-albedoPianeta)*costanteIrraggiamento/sigma, 0.25)
            # temperatura margherite
            temperaturaMargheriteNere = temperaturaMediaPianeta + isolamento*(albedoPianeta-albedoNere)
            temperaturaMargheriteBianche = temperaturaMediaPianeta + isolamento*(albedoPianeta-albedoBianche) 

            # Determina tassi di crescita
            if (temperaturaMargheriteNere >= temperaturaMargheriteMinima
                    and temperaturaMargheriteNere <= temperaturaMargheriteMassima
                    and areaNere >= 0.01):
                tassoCrescitaMargheriteNere = 1 - 0.003265*(temperaturaIdealeMargheriteNere-temperaturaMargheriteNere)**2
            else:
                tassoCrescitaMargheriteNere = 0.0

            if (temperaturaMargheriteBianche >= temperaturaMargheriteMinima
                    and temperaturaMargheriteBianche <= temperaturaMargheriteMassima
                    and areaBianche >= 0.01):
                tassoCrescitaMargheriteBianche = 1 - 0.003265*(temperaturaIdealeMargheriteBianche-temperaturaMargheriteBianche)**2
            else:
                tassoCrescitaMargheriteBianche = 0.0

            # Cambiamenti negli areali
            cambiamentoAreaNere = areaNere*(tassoCrescitaMargheriteNere*areaVuoto-tassoDecrescita)
            cambiamentoAreaBianche = areaBianche*(tassoCrescitaMargheriteBianche*areaVuoto-tassoDecrescita)

            # Cambiamenti dalla precedente iterazione
            cambiamentoPrecedenteIterazioneNere = abs(cambiamentoAreaNere-cambiamentoAreaNerePrecedente)
            cambiamentoPrecedenteIterazioneBianche = abs(cambiamentoAreaBianche-cambiamentoAreaBianchePrecedente)

            # Aggiorna aree, stati econteggio iterazioni
            cambiamentoAreaNerePrecedente = cambiamentoAreaNere
            cambiamentoAreaBianchePrecedente = cambiamentoAreaBianche
            areaNere = areaNere+cambiamentoAreaNere
            areaBianche = areaBianche+cambiamentoAreaBianche 
            areaVuoto = 1-(areaNere+areaBianche)
            it += 1

        # Vettore stati (storia)
        areaNereStoria[j] = areaNere
        areaBiancheStoria[j] = areaBianche 
        areaVuotoStoria[j] = areaVuoto
        temperaturaPianetaStoria[j] = temperaturaMediaPianeta
        temperaturaSenzaVita[j] = np.power(((costanteIrraggiamento*livelloEnergiaSolare/sigma)*(1-albedoVuoto)), (1/4)) 

    fig, ax = plt.subplots(2, 1)
    ax[0].plot(livelliEnergiaSolare, 100*areaNereStoria, color='black', label='nere')
    ax[0].plot(livelliEnergiaSolare, 100*areaBiancheStoria, color='red', label='bianche')
    #ax[0].plot(livelliEnergiaSolare, 100*areaNereStoria +100*areaBiancheStoria, label='area totale')
    # ax[0].set_xlabel('luminosità solare')
    ax[0].set_ylabel('area (%)')
    ax[0].legend()

    ax[1].plot(livelliEnergiaSolare, temperaturaPianetaStoria-KELVIN_OFFSET, color='black', label="con vita")
    ax[1].plot(livelliEnergiaSolare, temperaturaSenzaVita -KELVIN_OFFSET, color='red', label="senza vita")
    ax[1].set_xlabel('luminosità solare')
    ax[1].set_ylabel('temperatura globale (C)')
    ax[1].legend()
    plt.show()