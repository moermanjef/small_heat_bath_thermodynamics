import random
import matplotlib.pyplot as plt
import numpy as np

vrijheidsgraden = 20
volumegrootte = 1    # gemiddeld aantal deeltjes in systeem
T0 = 1
simulatietijd = 200000

U = T0*vrijheidsgraden

energieniveaus = np.array([random.uniform(0, 1) for i in range (vrijheidsgraden)])
energieniveaus *= U/sum(energieniveaus)

tracker = np.zeros(simulatietijd)
nummers = np.array([i for i in range (vrijheidsgraden)])
aantallen = np.array([i for i in range (vrijheidsgraden+1)])

def combinatoriek(a):
    getallen = np.zeros(a+1)+1
    index = 0
    for c in range (a+1):
        getal = 1
        for i in range (a-c+1, a+1):
            getal *= i
        for i in range (1, c+1):
            getal /= i
        getallen[index] = getal
        index += 1
    
    return getallen


# het aantal deeltjes in het systeem hebben een gewogen voorkomen
gewichten = combinatoriek(vrijheidsgraden)* (volumegrootte/vrijheidsgraden)**(aantallen) * (1 - volumegrootte/vrijheidsgraden)**(vrijheidsgraden-aantallen)
print(gewichten)


for i in range (simulatietijd):
    fluctuaties = np.zeros(vrijheidsgraden)
    for j in range (vrijheidsgraden):
        warmtebadtemperatuur = (U-energieniveaus[j])/(vrijheidsgraden-1)
        fluctuatie = random.gauss(0, warmtebadtemperatuur)
        # minimale energie = 0
        if fluctuatie + energieniveaus[j] < 0:
            fluctuatie *= energieniveaus[j]/fluctuatie
        energieniveaus[j] += fluctuatie
    # respecteer energiebehoud
    energieniveaus *= U/sum(energieniveaus)
    
    # een systeem waar er vrij deeltjes in en uit kunnen fluctueren
    aantal = random.choices(aantallen, weights=gewichten)
    aantal = aantal[0]
    deeltjes = np.random.choice(nummers, aantal, replace=False)
    for deeltje in deeltjes:
        tracker[i] += energieniveaus[deeltje]


def tsallis(x, vrijheidsgraden):
    return (vrijheidsgraden-1)/U * (1-x/U)**(vrijheidsgraden-2)

binaantal = 300

mogelijkheden = np.linspace(0, U, binaantal)
voorspelling = (vrijheidsgraden/volumegrootte)*simulatietijd*tsallis(mogelijkheden, vrijheidsgraden)/binaantal
exponentieel = (vrijheidsgraden/volumegrootte)*simulatietijd*np.exp(-mogelijkheden/(volumegrootte*T0))/(T0*(1-np.exp(-vrijheidsgraden))*binaantal)
    
plt.figure()
plt.title('Warmtebad met {} vrijheidsgraden \n T0: {}'.format(vrijheidsgraden, U/vrijheidsgraden), size = 20)
plt.xlabel('energietoestand van systeem met gemiddeld {} deeltjes'.format(volumegrootte), size = 14)
plt.ylabel('dichtheid', size = 14)
plt.grid()
plt.hist(tracker, bins=mogelijkheden, label="gesampelde energieen")
plt.plot(mogelijkheden, voorspelling, label="Tsallisdistributie")
plt.plot(mogelijkheden, exponentieel, label="Exponentiele")
plt.legend()
plt.show()

