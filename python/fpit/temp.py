

a = {13: -0.037146231792860515, 15: -0.037111590757275444, 16: 0.0, 36: 0.0, 38: 0.0, 41: 0.0, 62: 0.0, 63: 0.0, 87: 0.0, 88: -0.03719595413234706, 112: -0.036848710612285604}

a.values()
a.keys()
Aarr = np.zeros(max(a.keys())+1)

Aarr[np.array(list(a.keys()))] = np.array(list(a.values()))